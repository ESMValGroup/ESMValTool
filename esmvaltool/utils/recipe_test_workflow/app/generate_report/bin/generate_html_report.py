#!/usr/bin/env python
import os
import sqlite3
from datetime import datetime
from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape

from bin.commits_via_git import get_commits_from_git
from bin.sha_via_singularity import get_shas_from_singularity

# Load environment variables required at all sites.
CYLC_DB_PATH = os.environ.get("CYLC_DB_PATH")
CYLC_TASK_CYCLE_POINT = os.environ.get("CYLC_TASK_CYCLE_POINT")
CYLC_TASK_CYCLE_YESTERDAY = os.environ.get("ROSE_DATACP1D")
REPORT_PATH = os.environ.get("REPORT_PATH")
SITE = os.environ.get("SITE")

if SITE == "dkrz":
    ESMVAL_VERSIONS_TODAY = os.environ.get("ESMVAL_VERSIONS_CURRENT")
    ESMVAL_VERSIONS_YESTERDAY = os.environ.get("ESMVAL_VERSIONS_PREVIOUS")

elif SITE == "metoffice":
    REPOS = {
        "core_today": os.environ.get("ESMVALCORE_DIR"),
        "tool_today": Path(CYLC_TASK_CYCLE_YESTERDAY) / "ESMValCore",
        "core_yesterday": os.environ.get("ESMVALTOOL_DIR"),
        "tool_yesterday": Path(CYLC_TASK_CYCLE_YESTERDAY) / "ESMValTool",
    }
    print("Repos", REPOS)

SQL_QUERY_TASK_STATES = "SELECT name, status FROM task_states"


def main(db_file_path=CYLC_DB_PATH, site=SITE):
    """
    Main function to generate the HTML report.

    Parameters
    ----------
    db_file_path : str, default CYLC_DB_FILE_PATH
        The path to the SQLite database file.
    """
    raw_db_data = fetch_report_data(db_file_path)
    processed_db_data = process_db_output(raw_db_data)
    subheader = create_subheader()

    try:
        if site == "dkrz":
            # A single commit SHA for each package is expected.
            commit_info = get_shas_from_singularity(
                ESMVAL_VERSIONS_TODAY, ESMVAL_VERSIONS_YESTERDAY
            )
        elif site == "metoffice":
            # At least a single commit for each package is expected. There may
            # be multiple commits per package, and info for each commit has
            # multiple fields.
            commit_info = get_commits_from_git(REPOS)
            print("commit_info_messages", commit_info)
            add_report_message_to_git_commits(commit_info)
            print("commit_info_messages", commit_info)
        else:
            # No commit information for either package.
            commit_info = None
    # Catch as likely indicate a minor issue e.g. unexpected data content at
    # some point in the pipeline.
    except (ValueError, KeyError) as err:
        "Report generating without commit data. Error while fetching "
        f"commit data: {err}"
        commit_info = None

    rendered_html = render_html_report(
        subheader=subheader,
        report_data=processed_db_data,
        commit_info=commit_info,
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


def add_report_message_to_git_commits(git_commits_info):
    """
    Add report messages to a git commit information dictionary.

    Parameters
    ----------
    list[dict]
        A list of git commits.
    """
    git_commits_info[0]["report_flag"] = "Version tested this cycle >>>"
    if len(git_commits_info) > 1:
        git_commits_info[-1]["report_flag"] = "Version tested last cycle >>>"


def render_html_report(
    report_data,
    subheader,
    commit_info,
):
    """
    Render the HTML report using Jinja2.

    Parameters
    ----------
    report_data : dict
        The report data to be rendered in the HTML template.
    subheader : str
        The subheader for the HTML report.
    package_info : dict

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
        esmval_core_commits=commit_info["ESMValCore"]["commits"],
        esmval_tool_commits=commit_info["ESMValTool"]["commits"],
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
