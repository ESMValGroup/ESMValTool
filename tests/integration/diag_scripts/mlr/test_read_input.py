"""Tests for reading data."""

import os
from unittest import mock

import pytest
import yaml

from esmvaltool.diag_scripts.mlr.models import MLRModel

EXCEPTIONS = {
    'ValueError': ValueError,
    'TypeError': TypeError,
}


def get_call_name(call):
    """Get name of a `mock.call` function."""
    call_str = str(call)
    call_str = call_str[call_str.find('call.') + len('call.'):]
    call_str = call_str[:call_str.find('(')]
    return call_str


def get_logger_msg(method_calls):
    """Get all important logger calls."""
    all_calls = [get_call_name(call) for call in method_calls]
    return [call for call in all_calls if call in ('warning', 'error')]


class SimplifiedMLRModel(MLRModel):
    """Test class to avoid calling the base class `__init__` method."""

    def __init__(self, cfg):
        """Very simplified constructor of the base class."""
        self._cfg = cfg
        self._data = {}
        self._datasets = {}
        self._classes = {}


with open(
        os.path.join(os.path.dirname(__file__), 'configs',
                     'test_load_input_datasets.yml')) as file_:
    CONFIG = yaml.safe_load(file_)


# TODO: Add tests for ancestors
@pytest.mark.parametrize('data', CONFIG)
@mock.patch('esmvaltool.diag_scripts.mlr.logger', autospec=True)
@mock.patch('esmvaltool.diag_scripts.mlr.models.logger', autospec=True)
def test_load_input_datasets(mock_models_logger, mock_mlr_logger, data):
    """Test loading of input datasets."""
    cfg = data['cfg']
    input_datasets = list(cfg['input_data'].values())
    output = data['output']
    mlr_model = SimplifiedMLRModel(cfg)

    # Load input datasets
    if 'EXCEPTION' in output:
        exc = output['EXCEPTION']
        with pytest.raises(EXCEPTIONS[exc['type']]) as exc_info:
            mlr_model._load_input_datasets(input_datasets)
            assert exc.get('value', '') in str(exc_info.value)
    else:
        mlr_model._load_input_datasets(input_datasets)
        assert mlr_model._datasets == output

    # Logger calls
    logger_calls = mock_models_logger.method_calls
    logger_calls.extend(mock_mlr_logger.method_calls)
    assert get_logger_msg(logger_calls) == data['logger']
