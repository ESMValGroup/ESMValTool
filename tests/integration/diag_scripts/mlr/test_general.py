"""General tests for the module :mod:`esmvaltool.diag_scripts.mlr.models`."""

import os
from unittest import mock

import pytest
import yaml

from esmvaltool.diag_scripts.mlr.models import MLRModel

# Load test configuration
with open(
        os.path.join(os.path.dirname(__file__), 'configs',
                     'test_general.yml')) as file_:
    CONFIG = yaml.safe_load(file_)


@mock.patch('esmvaltool.diag_scripts.mlr.models.logger', autospec=True)
class TestMLRModel():
    """Tests for the base class."""

    args = CONFIG['args']
    kwargs = CONFIG['kwargs']

    def test_direct_initialization(self, mock_logger):
        """Test direct initialization without factory function."""
        with pytest.raises(NotImplementedError):
            MLRModel(*self.args, **self.kwargs)
        assert mock_logger.mock_calls == []

    def test_register_mlr_model(self, mock_logger):
        """Test registering subclass."""
        MLRModel._MODELS = {}
        assert MLRModel._MODELS == {}

        @MLRModel.register_mlr_model('test_model')
        class MyMLRModel(MLRModel):
            """Subclass of `MLRModel`."""

        assert MLRModel._MODELS == {'test_model': MyMLRModel}
        assert MyMLRModel._MLR_MODEL_TYPE == 'test_model'
        mock_logger.debug.assert_called_once()
        MLRModel._MODELS = {}

    @mock.patch.object(MLRModel, '_load_mlr_models')
    @mock.patch.object(MLRModel, '__init__', autospec=True)
    def test_create(self, mock_mlr_model_init, mock_load_mlr_models,
                    mock_logger):
        """Test creating subclasses."""
        # No subclasses
        MLRModel._MODELS = {}
        assert MLRModel._MODELS == {}
        mock_mlr_model_init.return_value = None
        with pytest.raises(NotImplementedError):
            MLRModel.create('test_model', *self.args, **self.kwargs)
        mock_load_mlr_models.assert_called()
        mock_mlr_model_init.assert_not_called()
        mock_load_mlr_models.reset_mock()
        mock_mlr_model_init.reset_mock()
        mock_logger.reset_mock()

        # Wrong subclass
        @MLRModel.register_mlr_model('test_model')
        class MyMLRModel(MLRModel):
            """Subclass of `MLRModel`."""

        with pytest.raises(NotImplementedError):
            MLRModel.create('another_test_model', *self.args, **self.kwargs)
        mock_load_mlr_models.assert_called()
        mock_mlr_model_init.assert_not_called()
        mock_load_mlr_models.reset_mock()
        mock_mlr_model_init.reset_mock()
        mock_logger.reset_mock()
        MLRModel._MODELS = {}

        # Right subclass
        @MLRModel.register_mlr_model('test_model')
        class MyMLRModel1(MLRModel):
            """Subclass of `MLRModel`."""

        MLRModel.create('test_model', *self.args, **self.kwargs)
        mock_load_mlr_models.assert_called()
        mock_mlr_model_init.assert_called_with(mock.ANY, *self.args,
                                               **self.kwargs)
        mock_logger.info.assert_called()
        mock_load_mlr_models.reset_mock()
        mock_mlr_model_init.reset_mock()
        mock_logger.reset_mock()
        MLRModel._MODELS = {}
