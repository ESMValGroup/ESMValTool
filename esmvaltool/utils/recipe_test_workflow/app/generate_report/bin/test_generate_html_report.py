import sqlite3
import tempfile
from contextlib import contextmanager
from pathlib import Path

import pytest

from esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.generate_html_report import (
    create_subheader,
    fetch_report_data,
    process_db_output,
    render_html_report,
)


@contextmanager
def mock_db_with_passed_values(row_data):
    """
    Context manager yielding a path to a temp db populated with passed values.

    Parameters
    ----------
    row_data : list[tuple]
        The mock task state data in the form [('recipe', 'status', 'cycle')]

    Yields
    ------
    Path
        A pathlib.Path to the mock db.
    """
    sql_string = ",".join(
        [f"('{row[0]}', '{row[1]}', '{row[2]}')" for row in row_data],
    )
    # This nested context manager cleans up the resource, so no
    # "try/except/finally" clean up is required.
    with tempfile.TemporaryDirectory() as temp_directory:
        path_to_synthetic_db = Path(temp_directory) / "synthetic_cylc_db"
        connection = sqlite3.connect(path_to_synthetic_db)
        cursor = connection.cursor()
        cursor.execute("CREATE TABLE task_states(name, status, cycle)")
        cursor.execute(f"INSERT INTO task_states VALUES {sql_string}")
        connection.commit()
        yield path_to_synthetic_db


def test_fetch_report_data_single_cycle():
    mock_cycle = "20250521T0100Z"
    mock_data = [
        ("process_recipe_1", "succeeded", mock_cycle),
        ("compare_recipe_1", "succeeded", mock_cycle),
        ("process_recipe_2", "succeeded", mock_cycle),
        ("compare_recipe_2", "failed", mock_cycle),
    ]
    expected = [
        ("process_recipe_1", "succeeded"),
        ("compare_recipe_1", "succeeded"),
        ("process_recipe_2", "succeeded"),
        ("compare_recipe_2", "failed"),
    ]
    with mock_db_with_passed_values(mock_data) as mock_cylc_db:
        actual = fetch_report_data(mock_cylc_db, mock_cycle)
    assert actual == expected


def test_fetch_report_data_multi_cycle():
    mock_cycle = "20250521T0100Z"
    mock_cycle_minus_1d = "20250520T0100Z"
    mock_cycle_minus_2d = "20250519T1700Z"

    # The RTW is run at 5pm on May, 19th. It runs automatically at 1am on May,
    # 20th and May 21st. On May 21st, the process task fails and no compare
    # task is run.
    mock_data = [
        ("process_recipe", "succeeded", mock_cycle_minus_2d),
        ("compare_recipe", "succeeded", mock_cycle_minus_2d),
        ("process_recipe", "succeeded", mock_cycle_minus_1d),
        ("compare_recipe", "succeeded", mock_cycle_minus_1d),
        ("process_recipe", "succeeded", mock_cycle),
    ]
    expected = [
        ("process_recipe", "succeeded"),
    ]
    with mock_db_with_passed_values(mock_data) as mock_cylc_db:
        actual = fetch_report_data(mock_cylc_db, mock_cycle)
    assert actual == expected


@pytest.mark.parametrize(
    "mock_db_output, expected",
    [
        (
            [("process_recipe_1", "succeeded")],
            {
                "recipe_1": {
                    "process_task": {
                        "status": "succeeded",
                        "style": "color: green",
                    },
                },
            },
        ),
        (
            [("compare_recipe_1", "failed")],
            {
                "recipe_1": {
                    "compare_task": {
                        "status": "failed",
                        "style": "color: red",
                    },
                },
            },
        ),
        (
            [
                ("process_recipe_1", "succeeded"),
                ("compare_recipe_1", "failed"),
            ],
            {
                "recipe_1": {
                    "process_task": {
                        "status": "succeeded",
                        "style": "color: green",
                    },
                    "compare_task": {
                        "status": "failed",
                        "style": "color: red",
                    },
                },
            },
        ),
        ([("other_task", "succeeded")], {}),
        (
            [("process_recipe_1", "running")],
            {
                "recipe_1": {
                    "process_task": {
                        "status": "running",
                        "style": "color: black",
                    },
                },
            },
        ),
        (
            [("process_examples--recipe_python", "succeeded")],
            {
                "examples/recipe_python": {
                    "process_task": {
                        "status": "succeeded",
                        "style": "color: green",
                    },
                },
            },
        ),
    ],
)
def test_process_db_output(mock_db_output, expected):
    actual = process_db_output(mock_db_output)
    assert actual == expected


def test_process_db_output_sorting():
    mock_db_output = [
        ("process_c-ecipe", "succeeded"),
        ("compare_a-ecipe", "failed"),
        ("process_a-ecipe", "succeeded"),
        ("process_b-ecipe", "failed"),
        ("comapre_a-ecipe", "running"),
    ]
    actual = process_db_output(mock_db_output)
    assert list(actual.keys()) == ["a-ecipe", "b-ecipe", "c-ecipe"]


def test_create_subheader():
    mock_cylc_task_cycle_point = "20250101T0001Z"
    actual = create_subheader(mock_cylc_task_cycle_point)
    assert actual == "Cycle start: 2025-01-01 00:01 UTC"


def test_render_html_report():
    mock_subheader = "Cycle start: 2025-01-01 00:01 UTC"
    mock_report_data = {
        "recipe_1": {
            "process_task": {"status": "failed", "style": "color: red"},
        },
        "recipe_2": {
            "process_task": {"status": "succeeded", "style": "color: green"},
            "compare_task": {"status": "succeeded", "style": "color: green"},
        },
        "recipe_3": {
            "process_task": {"status": "succeeded", "style": "color: green"},
            "compare_task": {"status": "failed", "style": "color: red"},
        },
    }
    actual = render_html_report(mock_report_data, mock_subheader)
    expected = '<!doctype html>\n<html>\n    <head>\n        <title>Recipe Test Workflow</title>\n        <link rel="icon" href="https://esmvaltool.org/favicon.ico">\n    </head>\n    <style>\n        #recipes {\n            font-family: Arial, Helvetica, sans-serif;\n            border-collapse: collapse;\n            width: 100%;\n        }\n\n        #recipes td, #recipes th {\n            border: 1px solid #ddd;\n            padding: 8px;\n        }\n\n        #recipes tr:nth-child(even){background-color: #f2f2f2;}\n\n        #recipes tr:hover {background-color: #ddd;}\n\n        #recipes th {\n            padding-top: 12px;\n            padding-bottom: 12px;\n            text-align: left;\n            background-color: hsl(200, 50%, 50%);\n            color: white;\n        }\n    </style>\n    <body>\n        <table id="recipes">\n            <h1>Recipe Test Workflow - Last Run Status</h1>\n            <h2>Cycle start: 2025-01-01 00:01 UTC</h2>\n            <tr>\n                <th>Recipe</th>\n                <th>Recipe Run</th>\n                <th>Compare KGOs</th>\n            </tr>\n            \n            <tr>\n                <td>recipe_1</td>\n                <td style="color: red;">failed</td>\n                <td style="color: black;">-</td>\n            </tr>\n            \n            <tr>\n                <td>recipe_2</td>\n                <td style="color: green;">succeeded</td>\n                <td style="color: green;">succeeded</td>\n            </tr>\n            \n            <tr>\n                <td>recipe_3</td>\n                <td style="color: green;">succeeded</td>\n                <td style="color: red;">failed</td>\n            </tr>\n            \n        </table>\n    </body>\n    <footer>\n        <div style="text-align: center;">\n            <a\n                href="https://www.dkrz.de/en/about-en/contact/impressum"\n                target="_blank"\n            >Imprint</a> and\n            <a\n                href="https://www.dkrz.de/en/about-en/contact/en-datenschutzhinweise"\n                target="_blank"\n            >Privacy Policy</a>\n        </div>\n    </footer>\n</html>'
    assert actual == expected
