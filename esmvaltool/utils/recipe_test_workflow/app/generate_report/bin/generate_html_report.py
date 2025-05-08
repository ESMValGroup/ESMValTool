#!/usr/bin/env python
from collections import defaultdict
from datetime import datetime
import os
import sqlite3

from jinja2 import Environment, FileSystemLoader, select_autoescape


CYLC_DB_FILE_PATH = os.environ.get("CYLC_DB_FILE_PATH")
CYLC_TASK_CYCLE_POINT = os.environ.get("CYLC_TASK_CYCLE_POINT")
CYLC_WORKFLOW_SHARE_DIR = os.environ.get("CYLC_WORKFLOW_SHARE_DIR")
OUTPUT_FILE_PATH = os.path.join(
    CYLC_WORKFLOW_SHARE_DIR, "recipe_test_workflow_status_report.html"
)

SQL_QUERY_TASK_STATES = "SELECT name, status FROM task_states"


def main():
    """
    Main function to generate the HTML report.
    """
    raw_db_data = fetch_report_data()
    procesed_db_data = process_db_output(raw_db_data)
    subheader = create_subheader()
    rendered_html = render_html_report(
        subheader=subheader,
        report_data=procesed_db_data,
    )
    write_report_to_file(rendered_html)


def fetch_report_data(db_file_path = CYLC_DB_FILE_PATH):
    """
    Fetch report data from the Cylc SQLite database.

    Parameters
    ----------
    db_file_path : str, default CYLC_DB_FILE_PATH
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


def process_db_output(report_data):
    """
    Process the database output.

    Group process and compare tasks by recipe. Filter out other tasks. Remove
    task prefix and add style information. E.g.
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
    styles = {
        "succeeded": "color: green",
        "failed": "color: red",
    }

    processed_db_data = defaultdict(dict)
    for task_name, status in report_data:
        if task_name.startswith("process_"):
            recipe_name = task_name.removeprefix("process_")
            processed_db_data[recipe_name]["process_task"] = {
                "status": status,
                "style": styles[status],
            }
        elif task_name.startswith("compare_"):            
            recipe_name = task_name.removeprefix("compare_")
            processed_db_data[recipe_name]["compare_task"] = {
                "status": status,
                "style": styles[status],
            }
    return processed_db_data


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
        autoescape=select_autoescape()
    )
    template = env.get_template("report_template.jinja")
    rendered_html = template.render(
        subheader=subheader,
        report_data=report_data,
    ) 
    return rendered_html


def create_subheader(cylc_task_cycle_point = CYLC_TASK_CYCLE_POINT):
    """
    Create the subheader for the HTML report.

    Parameters
    ----------
    cylc_task_cycle_point : str, default CYLC_TASK_CYCLE_POINT
        The cycle point of the task as str in ISO8601 format.

    Returns
    -------
    str
        The formatted subheader string.
    """
    parsed_datetime = datetime.strptime(CYLC_TASK_CYCLE_POINT, "%Y%m%dT%H%MZ")
    formated_datetime = parsed_datetime.strftime("%Y-%m-%d %H:%M")
    subheader = (
        f"Cycle start: {formated_datetime} UTC"
    )
    return subheader


def write_report_to_file(rendered_html, output_file_path = OUTPUT_FILE_PATH):
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


if __name__ == "__main__":
    main()
