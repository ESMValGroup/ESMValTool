#!/usr/bin/env python
"""Generate a HTML summary report from a Cylc SQLite database."""

import os
import shutil
import sqlite3
import subprocess
import traceback
from datetime import datetime
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape
from requests.exceptions import ConnectionError, HTTPError, Timeout

# Import from the ESMValTool package for testing.
try:
    from esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.fetch_commit_info import (
        fetch_commit_details_from_github_api,
    )
    from esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.shas_via_git import (
        get_shas_from_git,
    )
    from esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.shas_via_singularity import (
        get_shas_from_singularity,
    )
# Import locally for running in Cylc.
except ImportError:
    from fetch_commit_info import fetch_commit_details_from_github_api
    from shas_via_git import get_shas_from_git
    from shas_via_singularity import get_shas_from_singularity


# Load environment variables required at all sites.
CYLC_DB_PATH = os.environ.get("CYLC_DB_PATH")
CYLC_SHARE_DIR = os.environ.get("CYLC_WORKFLOW_SHARE_DIR")
CYLC_TASK_CYCLE_POINT = os.environ.get("CYLC_TASK_CYCLE_POINT")
CYLC_TASK_CYCLE_YESTERDAY = os.environ.get("CYLC_TASK_CYCLE_YESTERDAY")
CYLC_WORKFLOW_RUN_DIR = os.environ.get("CYLC_WORKFLOW_RUN_DIR")
OUTPUT_DIR = os.environ.get("OUTPUT_DIR")
PRODUCTION = os.environ.get("PRODUCTION")
REPORT_PATH = os.environ.get("REPORT_PATH")
SITE = os.environ.get("SITE")

# TODO: Remove/move to main/DKRZ after development. PRODUCTION too?
VM_PATH = os.environ.get("VM_PATH")
VM_PATH = Path(CYLC_WORKFLOW_RUN_DIR) / "mock_vm"
MOCK_PRODUCTION = True
if VM_PATH:
    VM_DEBUG_LOG_DIR = Path(VM_PATH) / "debug_logs"
    if VM_DEBUG_LOG_DIR.exists():
        shutil.rmtree(VM_DEBUG_LOG_DIR)
    VM_DEBUG_LOG_DIR.mkdir(parents=True, exist_ok=True)

ESMVAL_VERSIONS_TODAY = None
ESMVAL_VERSIONS_YESTERDAY = None
REPOS = None

if SITE == "dkrz":
    ESMVAL_VERSIONS_TODAY = os.environ.get("ESMVAL_VERSIONS_CURRENT")
    ESMVAL_VERSIONS_YESTERDAY = os.environ.get("ESMVAL_VERSIONS_PREVIOUS")

if SITE == "metoffice":
    REPOS = {
        "ESMValCore_today": os.environ.get("ESMVALCORE_DIR"),
        "ESMValTool_today": os.environ.get("ESMVALTOOL_DIR"),
    }
    if CYLC_TASK_CYCLE_YESTERDAY:
        path_to_yesterdays_cycle = Path(CYLC_TASK_CYCLE_YESTERDAY)
        REPOS["ESMValCore_yesterday"] = path_to_yesterdays_cycle / "ESMValCore"
        REPOS["ESMValTool_yesterday"] = path_to_yesterdays_cycle / "ESMValTool"


def main(
    db_file_path=CYLC_DB_PATH,
    site=SITE,
    report_path=REPORT_PATH,
    cylc_task_cycle_point=CYLC_TASK_CYCLE_POINT,
    esmval_versions_today=ESMVAL_VERSIONS_TODAY,
    esmval_versions_yesterday=ESMVAL_VERSIONS_YESTERDAY,
    repos=REPOS,
):
    """
    Main function to generate the HTML report.

    Parameters
    ----------
    db_file_path : str, default CYLC_DB_FILE_PATH
        The path to the SQLite database file.
    site : str
        The site the Recipe Test Workflow is being run at.
    report_path : str
        The path to output the HTML report.
    cylc_task_cycle_point : str
        The cycle point of the task as a string in ISO8601 format.
    esmval_versions_today : str | None
        The path to today's singularity container, if the site uses a
        singularity container, or None.
    esmval_versions_yesterday : str | None
        The path to yesterday's singularity container, if the site uses a
        singularity container and it exists, or None.
    repos : dict[str, str] | None
        A dictionary of git repos if the site uses git repos, or None.
    """
    commit_info = None
    raw_db_data = fetch_report_data(db_file_path, cylc_task_cycle_point)
    processed_db_data = process_db_output(raw_db_data)

    # Commits/SHAs will only be included for these sites. The report will run
    # at other sites without commit/SHA information.
    try:
        if site == "dkrz":
            sha_info = get_shas_from_singularity(
                esmval_versions_today, esmval_versions_yesterday
            )
        elif site == "metoffice":
            sha_info = get_shas_from_git(repos)
    # Catch the following errors so the report generates without commit/SHA
    # information. These errors are either propagated on purpose or
    # indicate a probable minor issue.
    except (
        ValueError,
        KeyError,
        IndexError,
        HTTPError,
        Timeout,
        ConnectionError,
    ):
        print(
            "Report generating with results only. Error while fetching commit "
            "data. See std.err log for details."
        )
        traceback.print_exc()

    if sha_info:
        commit_info = fetch_commit_details_from_github_api(sha_info)
        add_report_messages_to_commits(commit_info)

    # TODO: move to dkrz section after development and use PRODUCTION.
    if VM_DEBUG_LOG_DIR and MOCK_PRODUCTION:
        debug_log_processor(processed_db_data)

    reinstate_backslashes_to_recipe_names(processed_db_data)
    subheader = create_subheader(cylc_task_cycle_point)

    rendered_html = render_html_report(
        subheader=subheader,
        report_data=processed_db_data,
        commit_info=commit_info,
    )

    write_report_to_file(rendered_html, report_path)


def add_report_messages_to_commits(commit_info):
    """
    Add report messages to a commit info dictionary in-place.

    Parameters
    ----------
    commit_info : dict[str, list[dict]]
        A dictionary where keys are package names and values are lists of
        commit details. E.g.
        {
            "ESMValCore": [
                {"sha": "abcd123", "message": "Fix bug", ...},
                {"sha": "efgh456", "message": "Add feature", ...}
            ],
            "ESMValTool": [
                {"sha": "ijkl789", "message": "Update docs", ...}
            ]
        }
    """
    if commit_info:
        for package_commits in commit_info.values():
            package_commits[0]["report_flag"] = "Commit tested this cycle >>>"
            if len(package_commits) > 1:
                package_commits[-1]["report_flag"] = (
                    "Commit tested last cycle >>>"
                )


def fetch_report_data(db_file_path, target_cycle_point):
    """
    Fetch report data for a single cycle from the Cylc SQLite database.

    Parameters
    ----------
    db_file_path : str
        The path to the SQLite database file.
    target_cycle_point : str
        The cycle point to collect data for.

    Returns
    -------
    list
        A list of tuples containing rows of a single cycle's data from the
        cylc database.
    """
    sql_query_task_states_target_cycle = (
        "SELECT name, status FROM task_states WHERE cycle = ?"
    )
    connection = sqlite3.connect(db_file_path)
    cursor = connection.cursor()
    cursor.execute(sql_query_task_states_target_cycle, (target_cycle_point,))
    fetched_data = cursor.fetchall()
    connection.close()
    return fetched_data


def process_db_task(task_name, status):
    """
    Process db output data for a single task.

    Create a dictionary in the format:
    ```
    "process_task": {
        "status": "succeeded",
        "style": "color: green"
    },
    ```

    Parameters
    ----------
    task_name: str
        The name of the cylc task.
    status: str
        The task status.

    Returns
    -------
    tuple[str, dict]
        A tuple containing the name of the recipe as a string and the task data
        as a dictionary.
    """
    styles = {
        "succeeded": "color: green",
        "failed": "color: red",
    }
    task_name_parts = task_name.split("_", 1)
    recipe_name = task_name_parts[1]
    processed_task_name = task_name_parts[0] + "_task"
    # Restore directories to a "/"
    # recipe_name = task_name_parts[1].replace("--", "/")
    style = styles.get(status, "color: black")
    task_data = (
        recipe_name,
        {processed_task_name: {"status": status, "style": style}},
    )
    return task_data


def process_db_output(report_data):
    """
    Process the database output into a dictionary sorted by recipe name.

    Group target tasks by recipe. Filter out other tasks. Remove task prefix
    and add style information. Sort by recipe name. E.g.
    ```
    {
        "recipe_name": {
            "process_task": {
                "status": "succeeded",
                "style": "color: green"
            },
            "compare_task": {
                "status": "failed",
                "style": "color: red"
            }
        }
    }
    ```

    Parameters
    ----------
    report_data : list
        The report data fetched from the database. Each item is a tuple
        containing row data.

    Returns
    -------
    dict
        A dictionary with recipe names as keys and tasks/task data as values.
    """
    processed_db_data = {}
    # A tuple is required for the `startswith` func.
    tasks_to_include_in_report = ("process", "compare")
    for task_name, status in report_data:
        if task_name.startswith(tasks_to_include_in_report):
            recipe, task_data = process_db_task(task_name, status)
            if not processed_db_data.get(recipe):
                processed_db_data[recipe] = task_data
            else:
                processed_db_data[recipe].update(task_data)
    sorted_processed_db_data = dict(sorted(processed_db_data.items()))
    return sorted_processed_db_data


def copy_a_debug_log_to_vm(target_debug_log, recipe, task):
    """
    Copy a debug log file to a VM directory ``debug_logs/<recipe>/<task>``.

    Parameters
    ----------
    target_debug_log : Path
        The target log file. E.g. From the Cylc run directory.
    recipe : str
        The name of the recipe.
    task : str
        The name of the task e.g. "process_task", "compare_task".

    Returns
    -------
    Path
        The path to the debug log on the VM.
    """
    if not VM_DEBUG_LOG_DIR or not MOCK_PRODUCTION:
        return None
    file_name = target_debug_log.name
    vm_debug_log_dir_for_recipe_task = VM_DEBUG_LOG_DIR / recipe / task

    if not vm_debug_log_dir_for_recipe_task.exists():
        vm_debug_log_dir_for_recipe_task.mkdir(parents=True, exist_ok=True)

    command = f"rsync -a {target_debug_log} {vm_debug_log_dir_for_recipe_task}"
    subprocess.run(command, shell=True)
    return vm_debug_log_dir_for_recipe_task / file_name


def esmvaltool_debug_log_processor(recipe):
    """
    Copy the ESMValTool ``main_log_debug.txt`` to the VM, if it exists.

    Parameters
    ----------
    recipe: str
        The recipe.

    Returns
    -------
    dict
        A dict containing the path to the debug log file on the VM, or an empty
        dict.
    """
    if OUTPUT_DIR:
        output_dir = Path(OUTPUT_DIR)
        if output_dir.is_dir():
            # ESMValTool only uses the last part of a recipe name when creating
            # it's directory. E.g. ``examples--recipe_python`` will be in
            # ``recipe_python_<date>_<time>/run/``
            recipe_name_without_dir = recipe.split("--")[-1]
            for path in output_dir.iterdir():
                if path.name.startswith(recipe_name_without_dir):
                    debug_file_path = path / "run" / "main_log_debug.txt"
                    if debug_file_path.exists():
                        vm_debug_file_path = copy_a_debug_log_to_vm(
                            debug_file_path, recipe, "process_task"
                        )
                        return {"debug": vm_debug_file_path}
    return {}


def cylc_debug_log_processor(recipe, task):
    """
    Copy the Cylc stderr and stdout log files to the VM, if they exist.

    Parameters
    ----------
    recipe : str
        The name of the recipe.
    task : str
        The name of the task e.g. "process_task", "compare_task".

    Returns
    -------
    dict
        Dict of debug logs. The key is the HTML display name (e.g. ``stderr``)
        and the value is the path to the debug log on the VM. If no Cylc debug
        logs exist, the returned dict will be empty.
    """
    cylc_debug_logs = {}
    task_prefix = task.split("_")[0]
    for display_name, file in [("stderr", "job.err"), ("stdout", "job.out")]:
        path_to_run_cylc_log = (
            Path(CYLC_WORKFLOW_RUN_DIR)
            / "log"
            / "job"
            / CYLC_TASK_CYCLE_POINT
            / f"{task_prefix}_{recipe}"
            / "01"
            / file
        )
        if path_to_run_cylc_log.exists():
            vm_debug_file_path = copy_a_debug_log_to_vm(
                path_to_run_cylc_log, recipe, task
            )
            cylc_debug_logs[display_name] = vm_debug_file_path
    return cylc_debug_logs


def debug_log_processor(processed_db_data):
    """
    Copy debug logs to the VM and add the VM paths to the database task data.

    Parameters
    ----------
    processed_db_data : dict
        A dictionary with recipe names as keys and tasks/task data as values.
        Debug logs are added to the dict as task data e.g.
        ``{"<recipe>" : "<task>" {"status": ...  "debug_logs" { ...}}}``
    """
    for target_task in ("process_task", "compare_task"):
        for recipe, task_data in processed_db_data.items():
            target_task_data = task_data.get(target_task)
            if target_task_data:
                if target_task_data.get("status") == "failed":
                    target_task_data["debug_logs"] = {}
                    cylc_debug_log_paths = cylc_debug_log_processor(
                        recipe, target_task
                    )
                    target_task_data["debug_logs"].update(cylc_debug_log_paths)
                    # Only process tasks have ESMValTool debug logs.
                    if target_task == "process_task":
                        esmvaltool_debug_log_path = (
                            esmvaltool_debug_log_processor(recipe)
                        )
                        target_task_data["debug_logs"].update(
                            esmvaltool_debug_log_path
                        )


def reinstate_backslashes_to_recipe_names(processed_db_data):
    """
    Reinstate backslashes in recipe names for display in the report.

    Recipes nested in directories have backslashes converted to ``--`` for Cylc
    compatibility. This function reverts to the backslash for display in the
    report.

    Parameters
    ----------
    processed_db_data : dict
        A dictionary with recipe names as keys and tasks/task data as values.
        Debug logs are added to the dict as task data e.g.
        ``{"<recipe>" : "<task>" {"status": ...  "debug_logs" { ...}}}``
    """
    for recipe_name in processed_db_data.keys():
        if "--" in recipe_name:
            revised_recipe_name = recipe_name.replace("--", "/")
            processed_db_data[revised_recipe_name] = processed_db_data.pop(
                recipe_name
            )


def create_subheader(cylc_task_cycle_point):
    """
    Create the subheader for the HTML report.

    Parameters
    ----------
    cylc_task_cycle_point : str, default CYLC_TASK_CYCLE_POINT
        The cycle point of the task as a string in ISO8601 format.

    Returns
    -------
    str
        The formatted subheader string.
    """
    parsed_datetime = datetime.strptime(cylc_task_cycle_point, "%Y%m%dT%H%MZ")
    formatted_datetime = parsed_datetime.strftime("%Y-%m-%d %H:%M")
    subheader = f"Cycle start: {formatted_datetime} UTC"
    return subheader


def render_html_report(report_data, subheader, commit_info):
    """
    Render the HTML report using Jinja2.

    Parameters
    ----------
    report_data : dict
        The report data to be rendered in the HTML template.
    subheader : str
        The subheader for the HTML report.
    commit_info : dict
        The commit information for ESMValCore and ESMValTool.

    Returns
    -------
    str
        The rendered HTML content.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    env = Environment(
        loader=FileSystemLoader(script_dir),
        autoescape=select_autoescape(),
    )
    template = env.get_template("report_template.jinja")
    rendered_html = template.render(
        subheader=subheader,
        report_data=report_data,
        esmval_core_commits=commit_info.get("ESMValCore", []),
        esmval_tool_commits=commit_info.get("ESMValTool", []),
    )
    return rendered_html


def write_report_to_file(rendered_html, output_file_path):
    """
    Write the report data to an HTML file.

    Parameters
    ----------
    rendered_html : str
        The rendered HTML content.
    output_file_path : str
        The path to the output HTML file.
    """
    with open(output_file_path, "w") as file:
        file.write(rendered_html)
    print(f"HTML report written to: {output_file_path}")


if __name__ == "__main__":
    main()
