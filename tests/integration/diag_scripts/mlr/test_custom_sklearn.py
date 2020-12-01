"""Integration tests for :mod:`esmvaltool.diag_scripts.mlr.custom_sklearn`."""

import numpy as np
import pytest
from sklearn.base import clone
from sklearn.compose import TransformedTargetRegressor
from sklearn.decomposition import PCA
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

from esmvaltool.diag_scripts.mlr.custom_sklearn import (
    AdvancedPipeline,
    AdvancedTransformedTargetRegressor,
    _get_fit_parameters,
)

X_TRAIN = np.array([[3.0], [6.0], [10.0]])
Y_TRAIN = np.array([10.0, 20.0, 30.0])
np.set_printoptions(precision=10)


STEPS_1 = [('a', 1)]
STEPS_2 = [('a', 1), ('b', 0)]
TEST_GET_FIT_PARAMETERS = [
    ({'a': 1}, STEPS_1, ValueError),
    ({'a': 1, 'a__b': 1}, STEPS_1, ValueError),
    ({'a__x': 1}, [], ValueError),
    ({'a__x': 1}, STEPS_1, {'a': {'x': 1}}),
    ({'a__x': 1, 'a__y': 2}, STEPS_1, {'a': {'x': 1, 'y': 2}}),
    ({'a__x': 1, 'a__y__z': 2}, STEPS_1, {'a': {'x': 1, 'y__z': 2}}),
    ({'a__x': 1, 'b__y': 2}, STEPS_1, ValueError),
    ({'a__x': 1, 'b__y': 2}, STEPS_2, {'a': {'x': 1}, 'b': {'y': 2}}),
]


@pytest.mark.parametrize('kwargs,steps,output', TEST_GET_FIT_PARAMETERS)
def test_get_fit_parameters(kwargs, steps, output):
    """Test retrieving of fit parameters."""
    if isinstance(output, type):
        with pytest.raises(output):
            _get_fit_parameters(kwargs, steps, 'x')
        return
    params = _get_fit_parameters(kwargs, steps, 'x')
    assert params == output


class StdLinearRegression(LinearRegression):
    """Expand :class:`sklearn.linear_model.CoolLinearRegression`."""

    def predict(self, x, return_std=False):
        """Expand :meth:`predict`."""
        pred = super().predict(x)
        if return_std:
            err = np.ones(x.shape[0], dtype=x.dtype)
            return (pred, err)
        return pred


class VarLinearRegression(LinearRegression):
    """Expand :class:`sklearn.linear_model.CoolLinearRegression`."""

    def predict(self, x, return_var=False, return_cov=False, err_2d=False):
        """Expand :meth:`predict`."""
        pred = super().predict(x)
        if return_var:
            err = np.ones(x.shape[0], dtype=x.dtype)
            if err_2d:
                err = err.reshape(-1, 1)
            return (pred, err)
        if return_cov:
            err = np.ones((x.shape[0], x.shape[0]), dtype=x.dtype)
            return (pred, err)
        return pred


class NonStandardScaler(StandardScaler):
    """Expand :class:`sklearn.preprocessing.StandardScaler`."""

    def fit(self, x, y=None, f=0.0):
        """Expand :meth:`fit`."""
        return_value = super().fit(x, y)
        if self.mean_ is not None:
            self.mean_ += f
        return return_value


class TestAdvancedPipeline():
    """Tests for :class:`esmvaltool.diag_scripts.mlr.AdvancedPipeline`."""

    AREG = AdvancedTransformedTargetRegressor(
        transformer=NonStandardScaler(),
        regressor=LinearRegression(),
    )
    REG = TransformedTargetRegressor(
        transformer=NonStandardScaler(),
        regressor=LinearRegression(),
    )
    STEPS = [
        [('t', NonStandardScaler())],
        [('t', NonStandardScaler()), ('r', LinearRegression())],
        [('t', NonStandardScaler()), ('r', REG)],
        [('t', NonStandardScaler()), ('r', AREG)],
        [('t', NonStandardScaler()), ('r', AREG)],
        [('t', NonStandardScaler()), ('r', AREG)],
        [('t', NonStandardScaler()), ('r', AREG)],
    ]
    PIPELINES = [AdvancedPipeline(step) for step in STEPS]
    KW_X0 = {'a': 1, 't__f': 2.0}
    KW_X1 = {'b__a': 1, 't__f': 2.0}
    KW_X2 = {'t__wrongparam': 1, 't__f': 2.0}
    KW_X3 = {'r__wrongparam': 1, 't__f': 2.0}
    KW_X4 = {'r__wrongstep__f': 1, 't__f': 2.0}
    KW_X5 = {'r__regressor__wrongparam': 1, 't__f': 2.0}
    KW_0 = {'t__f': 2.0}
    KW_1 = {'t__f': 2.0, 'r__sample_weight': np.arange(3.0)}
    KW_2 = {'t__f': 2.0, 'r__transformer__f': 3.0}

    TEST_CHECK_FINAL_STEP = zip(
        PIPELINES,
        [TypeError, TypeError, TypeError, True, True, True, True, True],
    )

    @pytest.mark.parametrize('pipeline,output', TEST_CHECK_FINAL_STEP)
    def test_check_final_step(self, pipeline, output):
        """Test checking if final step."""
        pipeline = clone(pipeline)
        if isinstance(output, type):
            with pytest.raises(output):
                pipeline._check_final_step()
            return
        assert pipeline._check_final_step() is None

    TEST_FIT_TARGET_TRANSFORMER_ONLY = zip(
        PIPELINES,
        [{}, {}, {}, KW_X3, KW_X4, KW_0, KW_2],
        [TypeError,
         TypeError,
         TypeError,
         ValueError,
         ValueError,
         (np.array([20.0]), np.array([200.0 / 3.0])),
         NotImplementedError],
    )

    @pytest.mark.parametrize('pipeline,kwargs,output',
                             TEST_FIT_TARGET_TRANSFORMER_ONLY)
    def test_fit_target_transformer_only(self, pipeline, kwargs, output):
        """Test fitting of target transformer only."""
        pipeline = clone(pipeline)
        if isinstance(output, type):
            with pytest.raises(output):
                pipeline.fit_target_transformer_only(Y_TRAIN, **kwargs)
            return
        pipeline.fit_target_transformer_only(Y_TRAIN, **kwargs)
        transformer = pipeline.steps[-1][1].transformer_
        np.testing.assert_allclose(transformer.mean_, output[0])
        np.testing.assert_allclose(transformer.var_, output[1])
        assert not hasattr(pipeline.steps[-1][1], 'regressor_')
        with pytest.raises(NotFittedError):
            pipeline.predict(X_TRAIN)
        with pytest.raises(NotFittedError):
            pipeline.steps[-1][1].predict(X_TRAIN)

    TEST_FIT_TRANSFORMERS_ONLY = zip(
        PIPELINES,
        [KW_0, KW_0, KW_1, {}, KW_X0, KW_X1, KW_2],
        [None,
         (np.array([8.333333]), np.array([8.222222])),
         (np.array([8.333333]), np.array([8.222222])),
         (np.array([6.333333]), np.array([8.222222])),
         ValueError,
         ValueError,
         (np.array([8.333333]), np.array([8.222222]))],
    )

    @pytest.mark.parametrize('pipeline,kwargs,output',
                             TEST_FIT_TRANSFORMERS_ONLY)
    def test_fit_transformers_only(self, pipeline, kwargs, output):
        """Test fitting transformers only."""
        pipeline = clone(pipeline)
        if isinstance(output, type):
            with pytest.raises(output):
                pipeline.fit_transformers_only(X_TRAIN, Y_TRAIN, **kwargs)
            return
        pipeline.fit_transformers_only(X_TRAIN, Y_TRAIN, **kwargs)
        transformer = pipeline.steps[0][1]
        if output is None:
            assert not hasattr(transformer, 'mean_')
            assert not hasattr(transformer, 'var_')
            return
        np.testing.assert_allclose(transformer.mean_, output[0])
        np.testing.assert_allclose(transformer.var_, output[1])
        with pytest.raises(NotFittedError):
            pipeline.predict(X_TRAIN)
        with pytest.raises(NotFittedError):
            pipeline.steps[-1][1].predict(X_TRAIN)

    TEST_TRANSFORM_ONLY = [
        (KW_X0, ValueError),
        (KW_X1, KeyError),
        ({}, np.array([[-1.1624763874], [-0.1162476387], [1.2787240262]])),
        (KW_0, np.array([[-3.1624763874], [-2.1162476387], [-0.7212759738]])),
    ]

    @pytest.mark.parametrize('kwargs,output', TEST_TRANSFORM_ONLY)
    def test_transform_only(self, kwargs, output):
        """Test transforming only."""
        pipeline = AdvancedPipeline([
            ('s', StandardScaler()),
            ('t', NonStandardScaler()),
            ('r', LinearRegression()),
        ])
        with pytest.raises(NotFittedError):
            pipeline.transform_only(X_TRAIN)
        if isinstance(output, type):
            with pytest.raises(output):
                pipeline.fit(X_TRAIN, Y_TRAIN, **kwargs)
            return
        pipeline.fit(X_TRAIN, Y_TRAIN, **kwargs)
        x_trans = pipeline.transform_only(X_TRAIN)
        np.testing.assert_allclose(x_trans, output)

    TEST_TRANSFORM_TARGET_ONLY = zip(
        PIPELINES,
        [{}, {}, {}, {}, KW_X2, KW_0, KW_X5],
        [TypeError,
         TypeError,
         TypeError,
         np.array([-1.22474487, 0.0, 1.22474487]),
         np.array([-1.22474487, 0.0, 1.22474487]),
         np.array([-1.22474487, 0.0, 1.22474487]),
         np.array([-1.22474487, 0.0, 1.22474487])],
    )

    @pytest.mark.parametrize('pipeline,kwargs,output',
                             TEST_TRANSFORM_TARGET_ONLY)
    def test_transform_target_only(self, pipeline, kwargs, output):
        """Test transforming of target only."""
        pipeline = clone(pipeline)
        if isinstance(output, type):
            with pytest.raises(output):
                pipeline.fit_target_transformer_only(Y_TRAIN, **kwargs)
            return
        with pytest.raises(NotFittedError):
            pipeline.transform_target_only(Y_TRAIN)
        pipeline.fit_target_transformer_only(Y_TRAIN, **kwargs)
        y_trans = pipeline.transform_target_only(Y_TRAIN)
        np.testing.assert_allclose(y_trans, output)
        assert not hasattr(pipeline.steps[-1][1], 'regressor_')
        with pytest.raises(NotFittedError):
            pipeline.predict(X_TRAIN)
        with pytest.raises(NotFittedError):
            pipeline.steps[-1][1].predict(X_TRAIN)


class TestAdvancedTransformedTargetRegressor():
    """Tests for class ``AdvancedTransformedTargetRegressor``."""

    AREG = AdvancedTransformedTargetRegressor(
        transformer=NonStandardScaler(),
        regressor=LinearRegression(),
    )
    FIT_KWARGS = [
        {'a': 1},
        {'b__a': 1, 't__f': 2.0},
        {'regressor__wrongparam': 1},
        {'transformer__fails': 1, 'regressor__a': 1, 'regressor__b': 1},
        {},
        {'regressor__sample_weight': np.arange(3.0)},
    ]

    TEST_FIT = zip(
        FIT_KWARGS,
        [ValueError,
         ValueError,
         TypeError,
         NotImplementedError,
         (np.array([20.0]), np.array([200.0 / 3.0]), np.array([0.34756273]),
          -2.2012306472308283,
          np.array([10.54054054, 19.05405405, 30.40540541])),
         (np.array([20.0]), np.array([200.0 / 3.0]), np.array([0.30618622]),
          -1.8371173070873827, np.array([12.5, 20.0, 30.0]))],
    )

    @pytest.mark.parametrize('kwargs,output', TEST_FIT)
    def test_fit(self, kwargs, output):
        """Test fitting with kwargs."""
        reg = clone(self.AREG)
        if isinstance(output, type):
            with pytest.raises(output):
                reg.fit(X_TRAIN, Y_TRAIN, **kwargs)
            return
        reg.fit(X_TRAIN, Y_TRAIN, **kwargs)
        transformer = reg.transformer_
        regressor = reg.regressor_
        np.testing.assert_allclose(transformer.mean_, output[0])
        np.testing.assert_allclose(transformer.var_, output[1])
        np.testing.assert_allclose(regressor.coef_, output[2])
        np.testing.assert_allclose(regressor.intercept_, output[3])
        np.testing.assert_allclose(reg.predict(X_TRAIN), output[4])

    Y_2D = np.array([[10.0], [20.0], [30.0]])
    TEST_FIT_TRANSFORMER_ONLY = zip(
        FIT_KWARGS,
        [ValueError,
         ValueError,
         (Y_2D, {'wrongparam': 1}, np.array([20.0]), np.array([200.0 / 3.0])),
         NotImplementedError,
         (Y_2D, {}, np.array([20.0]), np.array([200.0 / 3.0])),
         (Y_2D,
          {'sample_weight': np.arange(3.0)},
          np.array([20.0]), np.array([200.0 / 3.0]))],
    )

    @pytest.mark.parametrize('kwargs,output', TEST_FIT_TRANSFORMER_ONLY)
    def test_fit_transformer_only(self, kwargs, output):
        """Test fitting of transformer only."""
        reg = clone(self.AREG)
        if isinstance(output, type):
            with pytest.raises(output):
                reg.fit_transformer_only(Y_TRAIN, **kwargs)
            return
        (y_2d, reg_kwargs) = reg.fit_transformer_only(Y_TRAIN, **kwargs)
        np.testing.assert_allclose(y_2d, output[0])
        assert isinstance(reg_kwargs, dict)
        assert reg_kwargs.keys() == output[1].keys()
        for (key, val) in reg_kwargs.items():
            np.testing.assert_allclose(val, output[1][key])
        transformer = reg.transformer_
        np.testing.assert_allclose(transformer.mean_, output[2])
        np.testing.assert_allclose(transformer.var_, output[3])
        assert not hasattr(reg, 'regressor_')
        with pytest.raises(NotFittedError):
            reg.predict(X_TRAIN)

    VAR_AREG = AdvancedTransformedTargetRegressor(
        transformer=NonStandardScaler(),
        regressor=VarLinearRegression(),
    )
    STD_AREG = AdvancedTransformedTargetRegressor(
        transformer=NonStandardScaler(),
        regressor=StdLinearRegression(),
    )
    REGS = [
        AREG,
        VAR_AREG,
        STD_AREG,
        AREG,
        VAR_AREG,
        STD_AREG,
        AREG,
        VAR_AREG,
        STD_AREG,
        AREG,
        VAR_AREG,
        STD_AREG,
        AREG,
        VAR_AREG,
        STD_AREG,
        AREG,
        VAR_AREG,
        STD_AREG,
        AREG,
        VAR_AREG,
        STD_AREG,
        AREG,
        VAR_AREG,
        STD_AREG,
    ]
    PREDICT_KWARGS = [
        {},
        {},
        {},
        {'wrong_kwarg': 1},
        {'wrong_kwarg': 1},
        {'wrong_kwarg': 1},
        {'always_return_1d': False},
        {'always_return_1d': False},
        {'always_return_1d': False},
        {'always_return_1d': False, 'return_std': True},
        {'always_return_1d': False, 'return_std': True},
        {'always_return_1d': False, 'return_std': True},
        {'always_return_1d': False, 'return_var': True},
        {'always_return_1d': False, 'return_var': True},
        {'always_return_1d': False, 'return_var': True},
        {'always_return_1d': False, 'return_var': True, 'err_2d': True},
        {'always_return_1d': False, 'return_var': True, 'err_2d': True},
        {'always_return_1d': False, 'return_var': True, 'err_2d': True},
        {'always_return_1d': False, 'return_cov': True},
        {'always_return_1d': False, 'return_cov': True},
        {'always_return_1d': False, 'return_cov': True},
        {'return_var': True, 'err_2d': True},
        {'return_var': True, 'err_2d': True},
        {'return_var': True, 'err_2d': True},
        {'return_var': True, 'return_cov': True},
        {'return_var': True, 'return_cov': True},
        {'return_var': True, 'return_cov': True},
    ]
    PREDS_1D = [
        np.array([10.5405405405, 19.0540540541, 30.4054054054]),
        np.array([12.5, 20.0, 30.0]),
    ]
    ERR = np.full(3, 200.0 / 3.0)
    COV = np.full((3, 3), 200.0 / 3.0)
    PRED_OUTPUT_1D = [
        (PREDS_1D, None),
        (PREDS_1D, None),
        (PREDS_1D, None),
        TypeError,
        TypeError,
        TypeError,
        (PREDS_1D, None),
        (PREDS_1D, None),
        (PREDS_1D, None),
        TypeError,
        TypeError,
        NotImplementedError,
        TypeError,
        (PREDS_1D, ERR),
        TypeError,
        TypeError,
        (PREDS_1D, ERR),
        TypeError,
        TypeError,
        (PREDS_1D, COV),
        TypeError,
        TypeError,
        (PREDS_1D, ERR),
        TypeError,
        RuntimeError,
        RuntimeError,
        RuntimeError,
    ]

    TEST_PREDICT_1D = zip(REGS, PREDICT_KWARGS, PRED_OUTPUT_1D)

    @pytest.mark.parametrize('reg,kwargs,output', TEST_PREDICT_1D)
    def test_predict_1d(self, reg, kwargs, output):
        """Test prediction."""
        for (idx, fit_kwargs) in enumerate(
                ({}, {'regressor__sample_weight': [0.0, 1.0, 1.0]})):
            new_reg = clone(reg)
            with pytest.raises(NotFittedError):
                new_reg.predict(X_TRAIN)
            new_reg.fit(X_TRAIN, Y_TRAIN, **fit_kwargs)
            if isinstance(output, type):
                with pytest.raises(output):
                    new_reg.predict(X_TRAIN, **kwargs)
                return
            y_pred = new_reg.predict(X_TRAIN, **kwargs)
            if output[1] is None:
                assert y_pred.shape == output[0][idx].shape
                np.testing.assert_allclose(y_pred, output[0][idx])
            else:
                assert y_pred[0].shape == output[0][idx].shape
                assert y_pred[1].shape == output[1].shape
                np.testing.assert_allclose(y_pred[0], output[0][idx])
                np.testing.assert_allclose(y_pred[1], output[1])

    VAR_AREG_1 = AdvancedTransformedTargetRegressor(
        transformer=NonStandardScaler(with_std=False),
        regressor=VarLinearRegression(),
    )
    PCA_AREG = AdvancedTransformedTargetRegressor(
        transformer=PCA(),
        regressor=VarLinearRegression(),
    )
    REGS = [
        AREG,
        VAR_AREG_1,
        PCA_AREG,
        AREG,
        VAR_AREG_1,
        PCA_AREG,
        AREG,
        VAR_AREG_1,
        PCA_AREG,
        AREG,
        VAR_AREG_1,
        PCA_AREG,
        AREG,
        VAR_AREG_1,
        PCA_AREG,
        AREG,
        VAR_AREG_1,
        PCA_AREG,
        AREG,
        VAR_AREG_1,
        PCA_AREG,
        AREG,
        VAR_AREG_1,
        PCA_AREG,
    ]
    PREDS_2D = [
        np.array([[10.5405405405], [19.0540540541], [30.4054054054]]),
        np.array([[12.5], [20.0], [30.0]]),
    ]
    ERR_1D = np.ones(3)
    ERR_2D = np.ones((3, 1))
    COV_1 = np.ones((3, 3))
    PRED_OUTPUT_2D = [
        (PREDS_1D, None),
        (PREDS_1D, None),
        (PREDS_1D, None),
        TypeError,
        TypeError,
        TypeError,
        (PREDS_2D, None),
        (PREDS_2D, None),
        (PREDS_2D, None),
        TypeError,
        TypeError,
        TypeError,
        TypeError,
        (PREDS_2D, ERR_1D),
        NotImplementedError,
        TypeError,
        (PREDS_2D, ERR_2D),
        NotImplementedError,
        TypeError,
        (PREDS_2D, COV_1),
        NotImplementedError,
        TypeError,
        (PREDS_1D, ERR_1D),
        NotImplementedError,
        RuntimeError,
        RuntimeError,
        RuntimeError,
    ]

    TEST_PREDICT_2D = zip(REGS, PREDICT_KWARGS, PRED_OUTPUT_2D)

    @pytest.mark.parametrize('reg,kwargs,output', TEST_PREDICT_2D)
    def test_predict_2d(self, reg, kwargs, output):
        """Test prediction."""
        y_train = Y_TRAIN.reshape(-1, 1)
        for (idx, fit_kwargs) in enumerate(
                ({}, {'regressor__sample_weight': [0.0, 1.0, 1.0]})):
            new_reg = clone(reg)
            with pytest.raises(NotFittedError):
                new_reg.predict(X_TRAIN)
            new_reg.fit(X_TRAIN, y_train, **fit_kwargs)
            if isinstance(output, type):
                with pytest.raises(output):
                    new_reg.predict(X_TRAIN, **kwargs)
                return
            y_pred = new_reg.predict(X_TRAIN, **kwargs)
            if output[1] is None:
                assert y_pred.shape == output[0][idx].shape
                np.testing.assert_allclose(y_pred, output[0][idx])
            else:
                assert y_pred[0].shape == output[0][idx].shape
                assert y_pred[1].shape == output[1].shape
                np.testing.assert_allclose(y_pred[0], output[0][idx])
                np.testing.assert_allclose(y_pred[1], output[1])

    TEST_GET_FIT_PARAMS = zip(
        FIT_KWARGS[:-1] + [{'regressor__a': 1, 'regressor__b': 2}],
        [ValueError,
         ValueError,
         ({}, {'wrongparam': 1}),
         NotImplementedError,
         ({}, {}),
         ({}, {'a': 1, 'b': 2})],
    )

    @pytest.mark.parametrize('kwargs,output', TEST_GET_FIT_PARAMS)
    def test_get_fit_params(self, kwargs, output):
        """Test retrieving of fit kwargs."""
        if isinstance(output, type):
            with pytest.raises(output):
                self.AREG._get_fit_params(kwargs)
            return
        fit_params = self.AREG._get_fit_params(kwargs)
        assert fit_params == output

    TEST_TO_BE_SQUEEZED = [
        (np.array([0]), True, 1, False),
        (np.array([0]), True, 2, False),
        (np.array([0]), False, 1, False),
        (np.array([0]), False, 2, False),
        (np.array([[0]]), True, 1, True),
        (np.array([[0]]), True, 2, True),
        (np.array([[0]]), False, 1, True),
        (np.array([[0]]), False, 2, False),
        (np.array([[0, 0], [0, 0]]), True, 1, False),
        (np.array([[0, 0], [0, 0]]), True, 2, False),
        (np.array([[0, 0], [0, 0]]), False, 1, False),
        (np.array([[0, 0], [0, 0]]), False, 2, False),
    ]

    @pytest.mark.parametrize('array,always_1d,training_dim,output',
                             TEST_TO_BE_SQUEEZED)
    def test_to_be_squeezed(self, array, always_1d, training_dim, output):
        """Test check if array should be squeezed."""
        reg = clone(self.AREG)
        reg._training_dim = training_dim
        squeezed = reg._to_be_squeezed(array, always_return_1d=always_1d)
        assert squeezed == output
