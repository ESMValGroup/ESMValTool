import pytest
from esmvalcore.config._config_validators import ValidationError

from esmvaltool.utils.recipe_test_workflow.app.configure.bin.configure import (
    validate_user_config_file,
)


def test_validate_user_config_file():
    mock_valid_config = {
        "output_dir": "~/esmvaltool_output",
        "auxiliary_data_dir": "~/auxiliary_data",
        "search_esgf": "never",
        "download_dir": "~/climate_data",
        "max_parallel_tasks": None,
        "log_level": "info",
        "exit_on_warning": True,
        "output_file_type": "png",
    }
    # No assert statement is needed. If the function call errors Pytest
    # considers the test failed.
    validate_user_config_file(mock_valid_config)


def test_validate_user_config_file_one_validation_error():
    mock_one_invalid_config = {
        "output_dir": "~/esmvaltool_output",
        "auxiliary_data_dir": "~/auxiliary_data",
        "search_esgf": "never",
        "download_dir": "~/climate_data",
        "max_parallel_tasks": None,
        "log_level": "info",
        "exit_on_warning": 100,
        "output_file_type": "png",
    }
    with pytest.raises(
        ValidationError,
        match='Validation error for EXIT_ON_WARNING with value "100"\n'
        "ERROR: Could not convert `100` to `bool`\n",
    ):
        validate_user_config_file(mock_one_invalid_config)


def test_validate_user_config_file_two_validation_errors():
    mock_two_invalids_config = {
        "output_dir": 111,
        "auxiliary_data_dir": "~/auxiliary_data",
        "search_esgf": "never",
        "download_dir": "~/climate_data",
        "max_parallel_tasks": None,
        "log_level": "info",
        "exit_on_warning": 100,
        "output_file_type": "png",
    }
    with pytest.raises(
        ValidationError,
        match='Validation error for OUTPUT_DIR with value "111"\nERROR: '
        "Expected a path, but got 111\n\nValidation error for "
        'EXIT_ON_WARNING with value "100"\nERROR: Could not convert `100` '
        "to `bool`\n",
    ):
        validate_user_config_file(mock_two_invalids_config)


def test_validate_user_config_file_key_error():
    mock_one_key_error = {
        "output_dir": "~/esmvaltool_output",
        "auxiliary_data_dir": "~/auxiliary_data",
        "search_esgf": "never",
        "download_dir": "~/climate_data",
        "one_rogue_field": 90210,
        "max_parallel_tasks": None,
        "log_level": "info",
        "exit_on_warning": True,
        "output_file_type": "png",
    }
    with pytest.raises(
        ValidationError,
        match="Key Error for ONE_ROGUE_FIELD. May not be a valid "
        "ESMValTool user configuration key\nERROR: 'one_rogue_field'\n",
    ):
        validate_user_config_file(mock_one_key_error)
