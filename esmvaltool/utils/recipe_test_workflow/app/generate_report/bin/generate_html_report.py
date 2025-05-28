#!/usr/bin/env python
import os
import sqlite3
from datetime import datetime

from jinja2 import Environment, FileSystemLoader, select_autoescape

CYLC_DB_PATH = os.environ.get("CYLC_DB_PATH")
CYLC_TASK_CYCLE_POINT = os.environ.get("CYLC_TASK_CYCLE_POINT")
CYLC_WORKFLOW_SHARE_DIR = os.environ.get("CYLC_WORKFLOW_SHARE_DIR")
REPORT_PATH = os.environ.get("REPORT_PATH")


def main(db_file_path=CYLC_DB_PATH):
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
    rendered_html = render_html_report(
        subheader=subheader,
        report_data=processed_db_data,
    )
    write_report_to_file(rendered_html)


def fetch_report_data(db_file_path, target_cycle_point=CYLC_TASK_CYCLE_POINT):
    """
    Fetch report data for a single cycle from the Cylc SQLite database.

    Parameters
    ----------
    db_file_path : str
        The path to the SQLite database file.
    target_cycle_point : str, default CYLC_TASK_CYCLE_POINT
        The cycle point to collect data for. Defaults to the current cylc
        cycle.

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


def render_html_report(report_data, subheader):
    """
    Render the HTML report using Jinja2.

    Parameters
    ----------
    report_data : dict
        The report data to be rendered in the HTML template.
    subheader : str
        The subheader for the HTML report.

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
