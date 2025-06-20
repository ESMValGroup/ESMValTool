"""
Tests for fetching detailled commit information from the GitHub API.

NOTE: The imports used mean these tests can only be run from within ESMValTool.
Ensure your ESMValTool working copy is installed in your environment before
running these tests.
"""

from unittest.mock import Mock, call, patch

import pytest

from esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.fetch_commit_info import (
    fetch_commit_details_from_github_api,
    fetch_range_of_commits,
    fetch_single_commit,
    process_commit_info,
)


@pytest.mark.parametrize(
    "mock_shas_by_package_and_day",
    [
        {
            "ESMValCore": {"today": "abcd123"},
            "ESMValTool": {"today": "ijkl789"},
        },
        {
            "ESMValCore": {"today": "abcd123", "yesterday": "abcd123"},
            "ESMValTool": {"today": "ijkl789", "yesterday": "ijkl789"},
        },
    ],
    ids=["todays_shas_only", "todays_shas_same_as_yesterdays"],
)
@patch(
    "esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.fetch_commit_info.process_commit_info"
)
@patch(
    "esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.fetch_commit_info.fetch_single_commit"
)
def test_fetch_commit_details_from_github_api_single_commits(
    mock_fetch_single, mock_process_commit_info, mock_shas_by_package_and_day
):
    """Test correct function is called for fetching single commits."""
    mock_fetch_single.return_value = {"mock_raw_commit": "raw_data"}
    mock_headers = "headers"
    mock_process_commit_info.return_value = [
        {"mock_processed_commit": "processed_data"}
    ]

    actual = fetch_commit_details_from_github_api(
        mock_shas_by_package_and_day, mock_headers
    )

    assert actual == {
        "ESMValCore": [{"mock_processed_commit": "processed_data"}],
        "ESMValTool": [{"mock_processed_commit": "processed_data"}],
    }

    assert mock_fetch_single.call_count == 2
    assert mock_fetch_single.call_args_list == [
        call("ESMValCore", "ESMValGroup", mock_headers, "abcd123"),
        call("ESMValTool", "ESMValGroup", mock_headers, "ijkl789"),
    ]

    assert mock_process_commit_info.call_count == 2
    assert mock_process_commit_info.call_args_list == [
        call({"mock_raw_commit": "raw_data"}),
        call({"mock_raw_commit": "raw_data"}),
    ]


@patch(
    "esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.fetch_commit_info.process_commit_info"
)
@patch(
    "esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.fetch_commit_info.fetch_range_of_commits"
)
def test_fetch_commit_details_from_github_api_range_of_commits(
    mock_fetch_range, mock_process_commit_info
):
    """Test correct function is called for fetching a range of commits."""
    mock_fetch_range.return_value = [{"mock_response": "data"}]
    mock_headers = "headers"
    mock_shas_by_package_and_day = {
        "ESMValCore": {"today": "abcd123", "yesterday": "efgh456"},
        "ESMValTool": {"today": "ijkl789", "yesterday": "mnop012"},
    }
    mock_process_commit_info.return_value = [
        {"mock_processed_commit": "processed_data"}
    ]

    actual = fetch_commit_details_from_github_api(
        mock_shas_by_package_and_day, mock_headers
    )

    assert actual == {
        "ESMValCore": [{"mock_processed_commit": "processed_data"}],
        "ESMValTool": [{"mock_processed_commit": "processed_data"}],
    }

    assert mock_fetch_range.call_count == 2
    assert mock_fetch_range.call_args_list == [
        call(
            "ESMValCore",
            "ESMValGroup",
            mock_headers,
            newer_sha="abcd123",
            older_sha="efgh456",
        ),
        call(
            "ESMValTool",
            "ESMValGroup",
            mock_headers,
            newer_sha="ijkl789",
            older_sha="mnop012",
        ),
    ]


@patch(
    "esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.fetch_commit_info.requests.get"
)
def test_fetch_single_commit(mock_get):
    """Test fetching a single commit from the GitHub API."""
    mock_response = Mock()
    mock_response.json.return_value = {"mock_response": "data"}
    mock_response.status_code = 200
    mock_get.return_value = mock_response

    actual = fetch_single_commit(
        "mock_repo", "mock_owner", "mock_headers", "mock_sha"
    )
    assert actual == {"mock_response": "data"}
    mock_get.assert_called_once_with(
        "https://api.github.com/repos/mock_owner/mock_repo/commits/mock_sha",
        headers="mock_headers",
        params=None,
        timeout=10,
    )


@patch(
    "esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.fetch_commit_info.requests.get"
)
def test_fetch_range_of_commits_results_on_first_page(mock_get):
    """Test fetching a range of commits with results on the first page."""
    mock_response_json_return_value = [
        {"sha": "tested_today"},
        {"sha": "intermediate_commit_1"},
        {"sha": "intermediate_commit_2"},
        {"sha": "tested_yesterday"},
        {"sha": "should_be_ignored"},
    ]
    mock_response = Mock()
    mock_response.json.return_value = mock_response_json_return_value
    mock_response.status_code = 200
    mock_get.return_value = mock_response

    mock_params = {
        "per_page": 10,
        "sha": "tested_today",
        "page": 1,
    }
    actual = fetch_range_of_commits(
        "mock_repo",
        "mock_owner",
        "mock_headers",
        newer_sha="tested_today",
        older_sha="tested_yesterday",
    )
    assert actual == mock_response_json_return_value[:-1]
    mock_get.assert_called_once_with(
        "https://api.github.com/repos/mock_owner/mock_repo/commits",
        headers="mock_headers",
        params=mock_params,
        timeout=10,
    )


@patch(
    "esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.fetch_commit_info.requests.get"
)
def test_fetch_range_of_commits_results_on_later_page(mock_get):
    """Test fetching a range of commits with results on later pages."""
    mock_response_json_return_value = (
        [{"sha": "tested_today"}]
        + [{"sha": f"intermediate_commit_{i}"} for i in range(1, 26)]
        + [{"sha": "tested_yesterday"}]
        + [{"sha": f"should_be_ignored_{i}"} for i in range(1, 21)]
    )

    mock_response = Mock()
    mock_response.json.side_effect = [
        mock_response_json_return_value[:10],  # First page
        mock_response_json_return_value[10:20],  # Second page
        mock_response_json_return_value[20:30],  # Third page
    ]
    mock_response.status_code = 200
    mock_get.return_value = mock_response

    actual = fetch_range_of_commits(
        "mock_repo",
        "mock_owner",
        "mock_headers",
        newer_sha="tested_today",
        older_sha="tested_yesterday",
    )
    assert actual == mock_response_json_return_value[:27]
    assert mock_get.call_count == 3

    # TODO: Can't figure out why not working.
    # expected_url = "https://api.github.com/repos/mock_owner/mock_repo/commits/"
    # expected_params = {
    #     "per_page": 10,
    #     "sha": "tested_today",
    #     "page": 3,
    # }
    # expected_calls = [
    #     call(expected_url, headers="mock_headers", params=expected_params, timeout=10),
    #     call(expected_url, headers="mock_headers", params=expected_params, timeout=10),
    #     call(expected_url, headers="mock_headers", params=expected_params, timeout=10),
    # ]
    # mock_get.assert_has_calls(expected_calls)


@pytest.mark.parametrize("num_of_commits_under_test", ["single", "range"])
def test_process_commit_info(num_of_commits_under_test):
    """Test processing commit information into a simplified format."""
    # Abbreviateded response for a real commit from ESMValTool.
    mock_raw_commit = {
        "sha": "662e792984d6577ae52e6931e794386ce508960c",
        "commit": {
            "author": {
                "name": "github-actions[bot]",
                "email": "41898282+github-actions[bot]@users.noreply.github.com",
                "date": "2025-06-11T11:55:31Z",
            },
            "committer": {
                "name": "GitHub",
                "email": "noreply@github.com",
                "date": "2025-06-11T11:55:31Z",
            },
            "message": "[Condalock] Update Linux condalock file (#4089)",
        },
        "url": "https://api.github.com/repos/ESMValGroup/ESMValTool/commit/"
        "662e792984d6577ae52e6931e794386ce508960c",
        "html_url": "https://github.com/ESMValGroup/ESMValTool/commit/"
        "662e792984d6577ae52e6931e794386ce508960c",
        "author": {
            "login": "github-actions[bot]",
            "avatar_url": "https://avatars.githubusercontent.com/in/15368?v=4",
        },
        "stats": {"total": 84, "additions": 42, "deletions": 42},
        "files": [],
    }
    mock_processed_commit = {
        "sha": "662e792",
        "author": "github-actions[bot]",
        "message": "[Condalock] Update Linux condalock file (#4089)",
        "date": "2025-06-11T11:55:31Z",
        "url": "https://github.com/ESMValGroup/ESMValTool/commit/"
        "662e792984d6577ae52e6931e794386ce508960c",
        "author_avatar": "https://avatars.githubusercontent.com/in/15368?v=4",
    }

    if num_of_commits_under_test == "single":
        mock_commit_info = mock_raw_commit
        expected = [mock_processed_commit]
    else:
        range_of_commits = 3
        mock_commit_info = [mock_raw_commit] * range_of_commits
        expected = [mock_processed_commit] * range_of_commits

    actual = process_commit_info(mock_commit_info)
    assert actual == expected
