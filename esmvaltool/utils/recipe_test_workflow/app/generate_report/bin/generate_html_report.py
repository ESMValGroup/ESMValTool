#!/usr/bin/env python
import os
import sqlite3
from datetime import datetime
from pathlib import Path
import subprocess

from jinja2 import Environment, FileSystemLoader, select_autoescape

# UNCOMMENT FOR LOCAL
# CYLC_DB_PATH = "hello"
# CYLC_TASK_CYCLE_POINT = "20250516T1053Z"
# REPORT_PATH="/home/users/christopher.billows/Code/ESMValTool/esmvaltool/utils/recipe_test_workflow/app/generate_report/bin/report.html"
# cached_raw_db_data = [
#     ('get_esmval', 'waiting'),
#     ('install_env_file','succeeded'),
#     ('get_esmval','succeeded'),
#     ('configure','succeeded'),
#     ('compare_recipe_radiation_budget','succeeded'),
#     ('process_recipe_radiation_budget','succeeded'),
#     ('process_recipe_albedolandcover','succeeded'),
#     ('process_recipe_ocean_amoc','succeeded'),
#     ('process_recipe_autoassess_landsurface_soilmoisture','succeeded'),
#     ('process_recipe_heatwaves_coldwaves','succeeded'),
#     ('process_recipe_ocean_multimap','succeeded'),
#     ('process_recipe_ensclus','succeeded'),
#     ('process_recipe_consecdrydays','succeeded'),
#     ('compare_recipe_albedolandcover','succeeded'),
#     ('compare_recipe_consecdrydays','succeeded'),
#     ('generate_report','running'),
#     ('compare_recipe_autoassess_landsurface_soilmoisture','succeeded'),
#     ('compare_recipe_heatwaves_coldwaves','succeeded'),
#     ('compare_recipe_ensclus','succeeded'),
#     ('compare_recipe_ocean_multimap','succeeded'),
#     ('compare_recipe_ocean_amoc','succeeded')
# ]


# UNCOMMENT FOR MET OFFICE / DKRZ
CYLC_DB_PATH = os.environ.get("CYLC_DB_PATH")
CYLC_TASK_CYCLE_POINT = os.environ.get("CYLC_TASK_CYCLE_POINT")
CYLC_TASK_PREVIOUS_CYCLE = os.environ.get("ROSE_DATACP1D")
REPORT_PATH = os.environ.get("REPORT_PATH")

# UNCOMMENT FOR MET OFFICE
# ESMVAL_CORE_CURRENT = os.environ.get("ESMVALCORE_DIR")
# ESMVAL_TOOL_CURRENT = os.environ.get("ESMVALTOOL_DIR")
# ESMVAL_CORE_PREVIOUS = Path(CYLC_TASK_PREVIOUS_CYCLE) / "ESMValCore"
# ESMVAL_TOOL_PREVIOUS = Path(CYLC_TASK_PREVIOUS_CYCLE) / "ESMValTool"

# UNCOMMENT FOR DKRZ
CONTAINER_DIR = os.environ.get("CONTAINER_DIR")
CONTAINER = "esmvaltool.sif"
CONTAINER_PATH = os.environ.get("CONTAINER_PATH")

SHARE_BIN = os.environ.get("SHARE_BIN")
ENV_FILE = os.environ.get("ENV_FILE")
SING_ENV_FILE = os.environ.get("SING_ENV_FILE")

SQL_QUERY_TASK_STATES = "SELECT name, status FROM task_states"


def main(db_file_path=CYLC_DB_PATH):
    """
    Main function to generate the HTML report.

    Parameters
    ----------
    db_file_path : str, default CYLC_DB_FILE_PATH
        The path to the SQLite database file.
    """
    # UNCOMMENT FOR LOCAL
    # raw_db_data = cached_raw_db_data
    # processed_db_data = process_db_output(raw_db_data)
    # local_esmvaltool = Path("/home/users/christopher.billows/Code/ESMValTool/")
    # local_esmvalcore = Path("/home/users/christopher.billows/Code/ESMValCore/")
    # esmval_core_previous_commit_sha = "170a93893"
    # esmval_tool_previous_commit_sha = "4515a2b92"
    # esmval_core_all_commits = fetch_git_commits(
    #     local_esmvalcore, esmval_core_previous_commit_sha
    # )
    # esmval_tool_all_commits = fetch_git_commits(
    #     local_esmvaltool, esmval_tool_previous_commit_sha
    # )

    # UNCOMMENT FOR MET OFFICE
    # raw_db_data = fetch_report_data(db_file_path)
    # processed_db_data = process_db_output(raw_db_data)

    # esmval_core_previous_commit_sha = None
    # if ESMVAL_CORE_PREVIOUS.exists():
    #     esmval_core_previous_commit_sha = (
    #         fetch_git_commits(ESMVAL_CORE_PREVIOUS)['sha']
    #     )

    # esmval_tool_previous_commit_sha = None
    # if ESMVAL_TOOL_PREVIOUS.exists():
    #     esmval_tool_previous_commit_sha = (
    #         fetch_git_commits(ESMVAL_TOOL_PREVIOUS)['sha']
    #     )

    # esmval_core_all_commits = fetch_git_commits(
    #     ESMVAL_CORE_CURRENT, esmval_core_previous_commit_sha
    # )
    # esmval_tool_all_commits = fetch_git_commits(
    #     ESMVAL_TOOL_CURRENT, esmval_tool_previous_commit_sha
    # )

    print("Container dir: ", CONTAINER_DIR)
    print("Container path: ", CONTAINER_PATH)
    print("Previous cycle point: ", CYLC_TASK_PREVIOUS_CYCLE)
    print("Share bin: ", SHARE_BIN)
    print("Env file: ", ENV_FILE)
    print("Singularity env: ", SING_ENV_FILE)

    print("Share bin exists?: ", Path(SHARE_BIN).exists())
    print("Env file exists?: ", Path(ENV_FILE).exists())
    print("Sing env exists?", Path(SING_ENV_FILE).exists())

    print("Env file content", Path(ENV_FILE).read_text())

    raw_db_data = fetch_report_data(db_file_path)
    processed_db_data = process_db_output(raw_db_data)

    print("Fetching package versions")
    current_package_versions = fetch_package_versions_from_container()
    print("Package versions fetched")

    subheader = create_subheader()
    rendered_html = render_html_report(
        subheader=subheader,
        report_data=processed_db_data,
        # esmval_core_commits=esmval_core_all_commits,
        # esmval_tool_commits=esmval_tool_all_commits,
    )
    write_report_to_file(rendered_html)


def fetch_report_data(db_file_path):
    """
    Fetch report data from the Cylc SQLite database.

    Parameters
    ----------
    db_file_path : str
        The path to the SQLite database file.

    Returns
    -------
    list
        A list of tuples containing rows of data from the database.
    """
    connection = sqlite3.connect(db_file_path)
    cursor = connection.cursor()
    cursor.execute(SQL_QUERY_TASK_STATES)
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


def fetch_package_versions_from_container():
    print("Cwd", CONTAINER_DIR)
    command = [SING_ENV_FILE, "singularity", "run", "esmvaltool.sif", "version"]
    print("Command: ", command)
    raw_version_info = subprocess.run(
        command,
        cwd=CONTAINER_DIR,
        capture_output=True,
        check=True,
        text=True
    )
    print(raw_version_info.stdout)
    print(raw_version_info.stderr)
    return raw_version_info.stdout


def add_report_message_to_git_commits(git_commits_info):
    """
    Add report messages to a git commit information dictionary.

    Parameters
    ----------
    list[dict]
        A list of git commits.
    """
    git_commits_info[0]['report_flag'] = "Version tested this cycle >>>"
    if len(git_commits_info) > 1:
        git_commits_info[-1]['report_flag'] = "Version tested last cycle >>>"


def fetch_git_commits(package_path, sha=None):
    """
    Fetch git commit information for an installed package.

    Parameters
    ----------
    package_path : str
        Path to a package's git repo.
    sha: str | None
        Optional. The sha of a previously tested commit. If provided, commits
        from HEAD back to the passed sha (inclusive) will be retrieved.

    Returns
    -------
    list[dict]
        A list of dicts where each dict represents one commit. If ``sha`` is
        passed, multiple commits/dicts may be returned.
    """
    command = [
        "git", "log", "-1", "--date=iso-strict", "--pretty=%cd^_^%h^_^%an^_^%s"
        ]

    if sha:
        command[2] = f"{sha}^..HEAD"

    raw_commit_info = subprocess.run(
        command, cwd=package_path, capture_output=True, check=True, text=True
    )

    processed_commit_info = []
    raw_commits = raw_commit_info.stdout.splitlines()
    for commit in raw_commits:
        split_fields = commit.split("^_^")
        processed_commit_info.append(
            {
                "report_flag": "",
                "date": split_fields[0],
                "sha": split_fields[1],
                "author": split_fields[2],
                "message": split_fields[3],
            }
        )
    add_report_message_to_git_commits(processed_commit_info)
    return processed_commit_info


def create_subheader(cylc_task_cycle_point=CYLC_TASK_CYCLE_POINT):
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
    formated_datetime = parsed_datetime.strftime("%Y-%m-%d %H:%M")
    subheader = f"Cycle start: {formated_datetime} UTC"
    return subheader


def render_html_report(
        report_data, subheader, esmval_core_commits, esmval_tool_commits,
    ):
    """
    Render the HTML report using Jinja2.

    Parameters
    ----------
    report_data : dict
        The report data to be rendered in the HTML template.
    subheader : str
        The subheader for the HTML report.
    esmval_core_commits : dict
        The ESMValCore commits information.
    esmval_tool_commits : dict
        The ESMValTool commits information.

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
        esmval_core_commits=esmval_core_commits,
        esmval_tool_commits=esmval_tool_commits,
    )
    return rendered_html


def write_report_to_file(rendered_html, output_file_path=REPORT_PATH):
    """
    Write the report data to an HTML file.

    Parameters
    ----------
    rendered_html : str
        The rendered HTML content.
    output_file_path : str, default OUTPUT_FILE_PATH
        The path to the output HTML file.
    """
    with open(output_file_path, "w") as file:
        file.write(rendered_html)
    print(f"HTML report written to: {output_file_path}")


if __name__ == "__main__":
    main()
