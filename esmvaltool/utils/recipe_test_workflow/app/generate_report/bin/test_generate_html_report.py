"""
Tests for generating a HTML report from a Cylc database.

NOTE: The imports used mean these tests can only be run from within ESMValTool.
Ensure your ESMValTool working copy is installed in your environment before
running these tests.
"""

import sqlite3
import tempfile
from collections import namedtuple
from contextlib import contextmanager
from pathlib import Path

import pytest

from esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.generate_html_report import (
    create_subheader,
    fetch_report_data,
    process_db_output,
)

MockDbData = namedtuple("MockDbdata", ["cycle", "row_data"])


@pytest.fixture
def mock_db_data_single_cycle():
    """Fixture providing a mock Cylc database with a single cycle of data."""
    cycle = "20250521T0100Z"
    row_data = [
        ("process_recipe_1", "succeeded", cycle),
        ("compare_recipe_1", "succeeded", cycle),
        ("process_recipe_2", "succeeded", cycle),
        ("compare_recipe_2", "failed", cycle),
    ]
    return MockDbData(cycle, row_data)


@contextmanager
def mock_db_with_passed_values(row_data):
    """
    Context manager yielding a path to a temp db populated with passed values.

    Parameters
    ----------
    row_data : list[tuple]
        The mock task state data in the form [('recipe', 'status', 'cycle')]

    Yields
    -------
    Path
        A pathlib.Path to the mock db.
    """
    sql_string = ",".join(
        [f"('{row[0]}', '{row[1]}', '{row[2]}')" for row in row_data]
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


def test_fetch_report_data_single_cycle(mock_db_data_single_cycle):
    """Test fetching report data from a DB containing a single cycle's data."""
    expected = [
        ("process_recipe_1", "succeeded"),
        ("compare_recipe_1", "succeeded"),
        ("process_recipe_2", "succeeded"),
        ("compare_recipe_2", "failed"),
    ]
    with mock_db_with_passed_values(
        mock_db_data_single_cycle.row_data
    ) as mock_cylc_db:
        actual = fetch_report_data(
            mock_cylc_db, mock_db_data_single_cycle.cycle
        )
    assert actual == expected


def test_fetch_report_data_multi_cycle():
    """Test fetching report data from a DB containing multiple cycles' data."""
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
                    "recipe_display_name": "recipe_1",
                    "process_task": {
                        "status": "succeeded",
                        "style": "color: green",
                    },
                }
            },
        ),
        (
            [("compare_recipe_1", "failed")],
            {
                "recipe_1": {
                    "recipe_display_name": "recipe_1",
                    "compare_task": {
                        "status": "failed",
                        "style": "color: red",
                    },
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
                    "recipe_display_name": "recipe_1",
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
                    "recipe_display_name": "recipe_1",
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
                "examples--recipe_python": {
                    "recipe_display_name": "examples/recipe_python",
                    "process_task": {
                        "status": "succeeded",
                        "style": "color: green",
                    },
                }
            },
        ),
    ],
)
def test_process_db_output(mock_db_output, expected):
    """Test processing the DB output into a structured report data format."""
    actual = process_db_output(mock_db_output)
    assert actual == expected


def test_process_db_output_sorting():
    """Test that the DB output is sorted by recipe name."""
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
    """Test creating a subheader with a formatted cycle point."""
    mock_cylc_task_cycle_point = "20250101T0001Z"
    actual = create_subheader(mock_cylc_task_cycle_point)
    assert actual == "Cycle start: 2025-01-01 00:01 UTC"
