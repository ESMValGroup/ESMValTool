import pytest

from esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.shas_via_singularity import (
    extract_scm_shas,
    get_shas_from_singularity,
)


@pytest.fixture()
def mock_scm_version_output():
    return (
        "ESMValCore: 2.13.0.dev54+g82d795ec\n"
        "ESMValTool: 2.13.0.dev66+g53c339c5c.d20250523"
    )


def test_get_shas_from_singularity(mock_scm_version_output):
    dev_version_today = mock_scm_version_output
    dev_version_yesterday = mock_scm_version_output
    actual = get_shas_from_singularity(
        dev_version_today, dev_version_yesterday
    )
    expected = {
        "core_today": "82d795ec",
        "tool_today": "53c339c5c",
        "core_yesterday": "82d795ec",
        "tool_yesterday": "53c339c5c",
    }
    assert actual == expected


def test_extract_scm_shas_valid(mock_scm_version_output):
    actual = extract_scm_shas(mock_scm_version_output, "today")
    assert actual == {"core_today": "82d795ec", "tool_today": "53c339c5c"}
