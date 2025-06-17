#!/usr/bin/env python
"""Generate a HTML summary report from a Cylc SQLite database."""

import os
import sqlite3
import traceback
from datetime import datetime
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape
from requests.exceptions import HTTPError

# Import from the ESMValTool package for testing.
try:
    from esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.fetch_commit_info import (
        fetch_commit_details_from_github_api,
    )
    from esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.shas_via_git import (
        CommitInfo,
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
CYLC_TASK_CYCLE_POINT = os.environ.get("CYLC_TASK_CYCLE_POINT")
CYLC_TASK_CYCLE_YESTERDAY = os.environ.get("CYLC_TASK_CYCLE_YESTERDAY")
REPORT_PATH = os.environ.get("REPORT_PATH")
SITE = os.environ.get("SITE")

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
    except (ValueError, KeyError, IndexError, HTTPError):
        print(
            "Report generating with results only. Error while fetching commit "
            "data. See std.err log for details."
        )
        traceback.print_exc()
    if sha_info:
        commit_info = fetch_commit_details_from_github_api(sha_info)
        print(commit_info)

    # raw_db_data = fetch_report_data(db_file_path, cylc_task_cycle_point)
    # processed_db_data = process_db_output(raw_db_data)
    # subheader = create_subheader(cylc_task_cycle_point)

    # rendered_html = render_html_report(
    #     subheader=subheader,
    #     report_data=processed_db_data,
    #     commit_info=commit_info,
    # )

    # write_report_to_file(rendered_html, report_path)


def add_report_messages_to_commits(commit_info):
    """
    Add report messages to a CommitInfo dataclass.

    Parameters
    ----------
    commit_info : CommitInfo
        A CommitInfo dataclass containing two lists of package commits.
    """
    for package_commits in [commit_info.core, commit_info.tool]:
        package_commits[0]["report_flag"] = "Version tested this cycle >>>"
        if len(package_commits) > 1:
            package_commits[-1]["report_flag"] = (
                "Version tested last cycle >>>"
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
    recipe_name = task_name_parts[1].replace("--", "/")
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


def copy_debug_log_to_vm(failed_task):
    """ """
    # TODO: Current work is in python-playground/lumberjack.py
    pass


def add_debug_log_for_failed_tasks(processed_db_data):
    """
    {
        "recipe_1": {
            "process_task": {
                "status": "succeeded",
                "style": "color: green",
                "debug_log": "path/to/debug.log",
            },
            "compare_task": {
                "status": "failed",
                "style": "color: red",
                "debug_log": "path/to/debug.log",
            },
        }
    },
    """

    for recipe, tasks in processed_db_data.items():
        for task_name, task_data in tasks.items():
            if task_data["status"] == "failed":
                # Assuming the debug log path is constructed from the recipe name
                # and task name. Adjust as necessary.
                debug_log_path = (
                    f"/path/to/debug/logs/{recipe}/{task_name}.log"
                )
                task_data["debug_log"] = debug_log_path


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


def render_html_report(report_data, subheader, commit_info, sha_info):
    """
    Render the HTML report using Jinja2.

    Parameters
    ----------
    report_data : dict
        The report data to be rendered in the HTML template.
    subheader : str
        The subheader for the HTML report.
    commit_info : CommitInfo | None
        The commit information for ESMValCore and ESMValTool, if the site uses
        git repos, or None.
    sha_info : dict | None
        The SHA information for ESMValCore and ESMValTool, if the site uses
        singularity containers, or None.

    Returns
    -------
    str
        The rendered HTML content.
    """
    commit_info = commit_info or CommitInfo([], [])
    script_dir = os.path.dirname(os.path.abspath(__file__))
    env = Environment(
        loader=FileSystemLoader(script_dir),
        autoescape=select_autoescape(),
    )
    template = env.get_template("report_template.jinja")
    rendered_html = template.render(
        subheader=subheader,
        report_data=report_data,
        esmval_core_commits=commit_info.core,
        esmval_tool_commits=commit_info.tool,
        sha_info=sha_info,
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
