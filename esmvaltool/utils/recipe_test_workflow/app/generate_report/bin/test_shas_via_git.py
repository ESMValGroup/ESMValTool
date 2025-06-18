"""
Tests for getting git commit SHAs via locally cloned git repositories.

NOTE: The imports used mean these tests can only be run from within ESMValTool.
Ensure your ESMValTool working copy is installed in your environment before
running these tests.
"""

import pytest

from esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.generate_html_report import (
    add_report_messages_to_commits,
)


@pytest.mark.skip(reason="To be reimplemented")
def test_add_report_message_to_git_commits_today_only():
    """Test adding report messages to git commits for today only."""
    mock_commit_info = (
        [{"date": "a", "sha": "abc", "author": "x", "message": "core"}],
        [{"date": "a", "sha": "xyz", "author": "y", "message": "tool"}],
    )
    add_report_messages_to_commits(mock_commit_info)
    expected = "Version tested this cycle >>>"
    assert mock_commit_info[0][0]["report_flag"] == expected
    assert mock_commit_info[1][0]["report_flag"] == expected


@pytest.mark.skip(reason="To be reimplemented")
def test_add_report_message_to_git_commits_both_days():
    """Test adding report messages to git commits for today and yesterday."""
    mock_commit_info = (
        [{"date": "a", "sha": "123", "author": "x", "message": "core"}],
        [{"date": "a", "sha": "xyz", "author": "y", "message": "tool"}],
    )
    add_report_messages_to_commits(mock_commit_info)
    expected = "Version tested this cycle >>>"
    assert mock_commit_info[0][0]["report_flag"] == expected
    assert mock_commit_info[1][0]["report_flag"] == expected


@pytest.mark.skip(reason="To be reimplemented")
def test_add_report_message_git_commits_both_days_multiple_commits():
    """Test adding report messages for both days with multiple commits."""
    mock_commit_info = (
        [
            {"date": "a", "sha": "123", "author": "x", "message": "core_1"},
            {"date": "a", "sha": "456", "author": "z", "message": "core_2"},
        ],
        [
            {"date": "b", "sha": "xyz", "author": "y", "message": "tool_1"},
            {"date": "b", "sha": "abc", "author": "z", "message": "tool_2"},
        ],
    )
    add_report_messages_to_commits(mock_commit_info)
    expected_1 = "Version tested this cycle >>>"
    expected_2 = "Version tested last cycle >>>"
    assert (
        mock_commit_info[0][0]["report_flag"]
        and mock_commit_info[1][0]["report_flag"]
    ) == expected_1
    assert (
        mock_commit_info[0][1]["report_flag"]
        and mock_commit_info[1][1]["report_flag"]
    ) == expected_2
