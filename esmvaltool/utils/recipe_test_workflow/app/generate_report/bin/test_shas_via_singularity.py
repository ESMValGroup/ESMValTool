import pytest

from esmvaltool.utils.recipe_test_workflow.app.generate_report.bin.shas_via_singularity import (
    get_shas_from_singularity,
    validate_all_shas,
)


def mock_scm_version_output():
    """
    A valid mock SCM version string.

    A function is used to allow safe in-test mutation and to allow the result
    to be passed as parametrized test values.

    Returns
    -------
    str
        A valid mock SCM version string.
    """
    return (
        "ESMValCore: 2.13.0.dev54+g82d795ec\n"
        "ESMValTool: 2.13.0.dev66+g53c339c5c\n"
    )


@pytest.mark.parametrize(
    "mock_day_version_today, mock_day_version_yesterday, expected",
    [
        (
            mock_scm_version_output(),
            None,
            {
                "ESMValCore": {"today": "82d795ec"},
                "ESMValTool": {"today": "53c339c5c"},
            },
        ),
        (
            mock_scm_version_output(),
            mock_scm_version_output(),
            {
                "ESMValCore": {"today": "82d795ec", "yesterday": "82d795ec"},
                "ESMValTool": {"today": "53c339c5c", "yesterday": "53c339c5c"},
            },
        ),
    ],
)
def test_get_shas_from_singularity_and_validate_for_valid_shas(
    mock_day_version_today, mock_day_version_yesterday, expected
):
    actual = get_shas_from_singularity(
        mock_day_version_today, mock_day_version_yesterday
    )
    # The unprocessed scm version strings are passed to the function purely
    # for error logging. Here 'None' is used.
    validate_all_shas(actual, None, None)
    assert actual == expected


@pytest.mark.parametrize(
    "shas, expected_message",
    [
        (
            {"ESMValCore": {}, "ESMValTool": {}},
            "Missing SHAs: dev_version_today=",
        ),
        (
            {"ESMValCore": {"today": "sha"}, "ESMValTool": {}},
            "Missing SHAs: dev_version_today=",
        ),
    ],
)
def test_validate_all_shas_for_invalid_shas(shas, expected_message):
    with pytest.raises(ValueError, match=expected_message):
        # The unprocessed scm version strings are passed to the function purely
        # for error logging. Here 'None' is used.
        validate_all_shas(shas, None, None)
