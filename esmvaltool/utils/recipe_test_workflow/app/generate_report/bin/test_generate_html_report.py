import os
import sqlite3
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.generate_html_report import (
    SQL_QUERY_TASK_STATES,
    create_subheader,
    fetch_report_data,
    main,
    process_db_output,
)
from esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.test_shas_via_singularity import (
    mock_scm_version_output,
)


@pytest.fixture()
def mock_cylc_db():
    with tempfile.TemporaryDirectory() as temp_directory:
        # Create a basic SQLite database file in a temporary directory.
        path_to_synthetic_db = Path(temp_directory) / "synthetic_cylc_db"
        connection = sqlite3.connect(path_to_synthetic_db)
        cursor = connection.cursor()
        cursor.execute("CREATE TABLE task_states(name, status)")
        cursor.execute("""
            INSERT INTO task_states VALUES
                ('process_recipe_1', 'succeeded'),
                ('compare_recipe_1', 'succeeded'),
                ('process_recipe_2', 'succeeded'),
                ('compare_recipe_2', 'failed')
        """)
        connection.commit()
        yield path_to_synthetic_db


def test_mock_cylc_db_and_sql_query(mock_cylc_db):
    connection = sqlite3.connect(mock_cylc_db)
    cursor = connection.cursor()
    actual = cursor.execute(SQL_QUERY_TASK_STATES)
    expected = [
        ("process_recipe_1", "succeeded"),
        ("compare_recipe_1", "succeeded"),
        ("process_recipe_2", "succeeded"),
        ("compare_recipe_2", "failed"),
    ]
    assert actual.fetchall() == expected


def test_main_for_site_dkrz(mock_cylc_db):
    # Duplicate input doesn't impact test.
    esmval_versions_today = mock_scm_version_output()
    esmval_version_yesterday = mock_scm_version_output()
    with patch.dict(os.environ, {"MY_ENV_VAR": "mocked_value"}):
        with tempfile.TemporaryDirectory() as temp_directory:
            actual = main(
                mock_cylc_db,
                "dkrz",
                temp_directory,
                "20250601T2015Z",
                esmval_versions_today,
                esmval_version_yesterday,
                None,
            )
            assert actual is None
            expected_report = Path(temp_directory) / "status_report.html"
            assert expected_report.exists()
            # report_content = expected_report.read_text()


def test_fetch_report_data(mock_cylc_db):
    actual = fetch_report_data(mock_cylc_db)
    expected = [
        ("process_recipe_1", "succeeded"),
        ("compare_recipe_1", "succeeded"),
        ("process_recipe_2", "succeeded"),
        ("compare_recipe_2", "failed"),
    ]
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
                    }
                }
            },
        ),
        (
            [("compare_recipe_1", "failed")],
            {
                "recipe_1": {
                    "compare_task": {"status": "failed", "style": "color: red"}
                }
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
                }
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
                }
            },
        ),
        (
            [("process_examples--recipe_python", "succeeded")],
            {
                "examples/recipe_python": {
                    "process_task": {
                        "status": "succeeded",
                        "style": "color: green",
                    }
                }
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


# def test_render_html_report_no_commits_no_shas():
#     mock_subheader = "Cycle start: 2025-01-01 00:01 UTC"
#     mock_report_data = {
#         "recipe_1": {
#             "process_task": {"status": "failed", "style": "color: red"},
#         },
#         "recipe_2": {
#             "process_task": {"status": "succeeded", "style": "color: green"},
#             "compare_task": {"status": "succeeded", "style": "color: green"},
#         },
#         "recipe_3": {
#             "process_task": {"status": "succeeded", "style": "color: green"},
#             "compare_task": {"status": "failed", "style": "color: red"},
#         },
#     }
#     actual = render_html_report(mock_report_data, mock_subheader)
#     expected = '<!doctype html>\n<html>\n    <head>\n        <title>Recipe Test Workflow</title>\n        <link rel="icon" href="https://esmvaltool.org/favicon.ico">\n    </head>\n    <style>\n        #recipes {\n            font-family: Arial, Helvetica, sans-serif;\n            border-collapse: collapse;\n            width: 100%;\n        }\n\n        #recipes td, #recipes th {\n            border: 1px solid #ddd;\n            padding: 8px;\n        }\n\n        #recipes tr:nth-child(even){background-color: #f2f2f2;}\n\n        #recipes tr:hover {background-color: #ddd;}\n\n        #recipes th {\n            padding-top: 12px;\n            padding-bottom: 12px;\n            text-align: left;\n            background-color: hsl(200, 50%, 50%);\n            color: white;\n        }\n    </style>\n    <body>\n        <table id="recipes">\n            <h1>Recipe Test Workflow - Last Run Status</h1>\n            <h2>Cycle start: 2025-01-01 00:01 UTC</h2>\n            <tr>\n                <th>Recipe</th>\n                <th>Recipe Run</th>\n                <th>Compare KGOs</th>\n            </tr>\n            \n            <tr>\n                <td>recipe_1</td>\n                <td style="color: red;">failed</td>\n                <td style="color: black;">-</td>\n            </tr>\n            \n            <tr>\n                <td>recipe_2</td>\n                <td style="color: green;">succeeded</td>\n                <td style="color: green;">succeeded</td>\n            </tr>\n            \n            <tr>\n                <td>recipe_3</td>\n                <td style="color: green;">succeeded</td>\n                <td style="color: red;">failed</td>\n            </tr>\n            \n        </table>\n    </body>\n    <footer>\n        <div style="text-align: center;">\n            <a\n                href="https://www.dkrz.de/en/about-en/contact/impressum"\n                target="_blank"\n            >Imprint</a> and\n            <a\n                href="https://www.dkrz.de/en/about-en/contact/en-datenschutzhinweise"\n                target="_blank"\n            >Privacy Policy</a>\n        </div>\n    </footer>\n</html>'
#     assert actual == expected


# def test_render_html_report_partial():
#     mock_subheader = "Cycle start: 2025-01-01 00:01 UTC"
#     mock_report_data = {
#         "recipe_1": {
#             "process_task": {"status": "failed", "style": "color: red"},
#         },
#         "recipe_2": {
#             "process_task": {"status": "succeeded", "style": "color: green"},
#             "compare_task": {"status": "succeeded", "style": "color: green"},
#         },
#         "recipe_3": {
#             "process_task": {"status": "succeeded", "style": "color: green"},
#             "compare_task": {"status": "failed", "style": "color: red"},
#         },
#     }
#     commit_info = None
#     actual = render_html_report(mock_report_data, mock_subheader, commit_info)
#     assert expected in actual
