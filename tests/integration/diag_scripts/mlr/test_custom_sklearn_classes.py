"""Integration tests for classes of custom :mod:`sklearn` functionalities.

Parts of this code have been copied from :mod:`sklearn`.

License: BSD 3-Clause License

Copyright (c) 2007-2020 The scikit-learn developers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

# pylint: disable=arguments-differ
# pylint: disable=attribute-defined-outside-init
# pylint: disable=invalid-name
# pylint: disable=no-self-use
# pylint: disable=protected-access
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-arguments

from copy import deepcopy

import numpy as np
import pytest
from sklearn.base import BaseEstimator, clone
from sklearn.compose import ColumnTransformer, TransformedTargetRegressor
from sklearn.decomposition import PCA
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LinearRegression
from sklearn.metrics import make_scorer, mean_absolute_error
from sklearn.preprocessing import FunctionTransformer, StandardScaler

from esmvaltool.diag_scripts.mlr.custom_sklearn import (
    _DEFAULT_TAGS,
    AdvancedPipeline,
    AdvancedRFE,
    AdvancedRFECV,
    AdvancedTransformedTargetRegressor,
    FeatureSelectionTransformer,
    _score_weighted,
)

# AdvancedPipeline


X_TRAIN = np.array([[3.0], [6.0], [10.0]])
Y_TRAIN = np.array([10.0, 20.0, 30.0])


class FeatureImportanceRegressor(BaseEstimator):
    """Estimator that has ``feature_importances_`` attribute."""

    def __init__(self):
        """Initialize instance."""
        super().__init__()
        self.feature_importances_ = 42

    def fit(self, *_):
        """Fit method."""
        return self


class StdLinearRegression(LinearRegression):
    """Expand :class:`sklearn.linear_model.LinearRegression`."""

    def predict(self, x, return_std=False):
        """Expand :meth:`predict`."""
        pred = super().predict(x)
        if return_std:
            err = np.ones(x.shape[0], dtype=x.dtype)
            return (pred, err)
        return pred


class VarLinearRegression(LinearRegression):
    """Expand :class:`sklearn.linear_model.LinearRegression`."""

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
    """Tests for ``AdvancedPipeline``."""

    def test_coef_(self):
        """Test ``coef_`` property."""
        pipeline = AdvancedPipeline(
            [('t', StandardScaler(with_std=False)), ('r', LinearRegression())],
        )
        pipeline.fit(np.arange(3).reshape(3, 1), np.arange(3))
        np.testing.assert_allclose(pipeline.coef_, [1.0])

    def test_feature_importances_(self):
        """Test ``feature_importances_`` property."""
        pipeline = AdvancedPipeline(
            [('t', StandardScaler()), ('r', FeatureImportanceRegressor())],
        )
        assert pipeline.feature_importances_ == 42

    def test_fit(self):
        """Test ``_fit``."""
        x_data = np.array([
            [0, 1000],
            [1, 0],
            [2, 3000],
            [0, -5000],
            [4, -3000],
            [4, -3000],
        ])
        y_data = np.array([1, 0, 3, -5, -3, -3])
        pipeline = AdvancedPipeline([
            ('t', StandardScaler()), ('r', LinearRegression()),
        ])
        sample_weights = np.array([0.0, 0.0, 0.0, 0.0, 1.0, 1.0])
        kwargs = {
            't': {'sample_weight': sample_weights},
            'r': {'sample_weight': sample_weights},
        }
        pipeline._fit(x_data, y_data, **kwargs)

        transformer_ = pipeline.steps[0][1]
        np.testing.assert_allclose(transformer_.scale_, [1.0, 1.0])
        np.testing.assert_allclose(transformer_.mean_, [4.0, -3000.0])

        regressor_ = pipeline.steps[1][1]
        with pytest.raises(NotFittedError):
            regressor_.predict([[0, 0]])

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


# AdvancedRFE


class NewLinearRegression(LinearRegression):
    """Expand ``LinearRegression``."""

    def predict(self, x_data, always_one=False):
        """Add dummy predict_kwargs to function."""
        if always_one:
            return 'one'
        return super().predict(x_data)


class TestAdvancedRFE():
    """Tests for ``AdvancedRFE``."""

    X_TRAIN = np.array(
        [[0.0, 0.0, 0.0],
         [2.0, 0.0, 1.0],
         [3.0, 0.0, -2.0]],
    )
    X_PRED = np.array(
        [[1000.0, 100.0, 10.0],
         [2000.0, 200.0, 20.0]],
    )

    Y_TRAIN = np.array([0.0, 1.0, -2.0])
    SAMPLE_WEIGHTS = np.array([1.0, 1.0, 0.0])

    def get_rfe(self, drop=False):
        """``AdvancedRFE`` object."""
        if drop:
            column_transformer_args = [[
                ('drop', 'drop', [2]),
                ('passthrough', 'passthrough', [0, 1]),
            ]]
        else:
            column_transformer_args = [[
                ('passthrough', 'passthrough', [0, 1, 2]),
            ]]
        pipeline_args = [[
            ('trans', ColumnTransformer(*column_transformer_args)),
            ('lin', NewLinearRegression()),
        ]]
        rfe_kwargs = {
            'estimator': AdvancedPipeline(*pipeline_args),
            'n_features_to_select': 1,
            'step': 1,
            'verbose': 1000,
        }
        return AdvancedRFE(**rfe_kwargs)

    @pytest.fixture
    def rfe(self):
        """``AdvancedRFE`` object."""
        return self.get_rfe(drop=False)

    @pytest.fixture
    def rfe_drop(self):
        """``AdvancedRFE`` object where features are dropped."""
        return self.get_rfe(drop=True)

    class NoCoefReg(BaseEstimator):
        """Estimator without ``coef_`` and ``feature_importances_``."""

        def fit(self, *_):
            """Fit method."""
            return self

    def test_advanced_rfe_fail(self, rfe_drop):
        """Test ``AdvancedRFE`` expected fail."""
        # Transformer that drops features
        with pytest.raises(NotImplementedError):
            rfe_drop.fit(self.X_TRAIN, self.Y_TRAIN)

        # Regressor without coef_ or feature_importances_
        msg = ("The classifier does not expose 'coef_' or "
               "'feature_importances_' attributes")
        fail_rfe = AdvancedRFE(self.NoCoefReg())
        with pytest.raises(RuntimeError, match=msg):
            fail_rfe.fit(np.arange(6).reshape(3, 2), np.arange(3))

        # Invalid step
        msg = "Step must be >0"
        fail_rfe = AdvancedRFE(LinearRegression(), step=-1)
        with pytest.raises(ValueError, match=msg):
            fail_rfe.fit(np.arange(6).reshape(3, 2), np.arange(3))

    class FIReg(BaseEstimator):
        """Estimator with working ``feature_importances_``."""

        def fit(self, x_data, *_):
            """Fit method."""
            self.feature_importances_ = np.full((4, x_data.shape[1]), 0.0)
            for idx in range(self.feature_importances_.shape[1]):
                self.feature_importances_[:, idx] = float(idx)
            print(self.feature_importances_)
            return self

    def test_feature_importances(self):
        """Test with ``feature_importances_``."""
        firfe = AdvancedRFE(self.FIReg())
        x_data = np.arange(3 * 6).reshape(3, 6)
        y_data = np.arange(3)
        firfe.fit(x_data, y_data)
        assert firfe.n_features_ == 3
        np.testing.assert_array_equal(firfe.ranking_, [4, 3, 2, 1, 1, 1])
        np.testing.assert_array_equal(firfe.support_,
                                      [False, False, False, True, True, True])

    def test_advanced_rfe_no_fit_kwargs(self, rfe):
        """Test ``AdvancedRFE`` without fit_kwargs."""
        rfe.fit(self.X_TRAIN, self.Y_TRAIN)
        assert rfe.n_features_ == 1
        np.testing.assert_array_equal(rfe.ranking_, [2, 3, 1])
        np.testing.assert_array_equal(rfe.support_, [False, False, True])
        est = rfe.estimator_
        assert isinstance(est, AdvancedPipeline)
        assert est.steps[0][1].transformers_ == [
            ('passthrough', 'passthrough', [0])]
        np.testing.assert_allclose(est.steps[1][1].coef_, [1.0])
        np.testing.assert_allclose(est.steps[1][1].intercept_, 0.0, atol=1e-10)
        pred = rfe.predict(self.X_PRED)
        np.testing.assert_allclose(pred, [10.0, 20.0])
        pred_one = rfe.predict(self.X_PRED, always_one=True)
        assert pred_one == 'one'

    def test_advanced_rfe_fit_kwargs(self, rfe):
        """Test ``AdvancedRFE`` with fit_kwargs."""
        rfe.fit(self.X_TRAIN, self.Y_TRAIN,
                lin__sample_weight=self.SAMPLE_WEIGHTS)
        assert rfe.n_features_ == 1
        np.testing.assert_array_equal(rfe.ranking_, [1, 3, 2])
        np.testing.assert_array_equal(rfe.support_, [True, False, False])
        est = rfe.estimator_
        assert isinstance(est, AdvancedPipeline)
        assert est.steps[0][1].transformers_ == [
            ('passthrough', 'passthrough', [0])]
        np.testing.assert_allclose(est.steps[1][1].coef_, [0.5])
        np.testing.assert_allclose(est.steps[1][1].intercept_, 0.0, atol=1e-10)
        pred = rfe.predict(self.X_PRED)
        np.testing.assert_allclose(pred, [500.0, 1000.0])
        pred_one = rfe.predict(self.X_PRED, always_one=True)
        assert pred_one == 'one'

    def step_score(self, estimator, features):
        """Score for a single step rfe."""
        x_test = np.arange(20).reshape(1, 20)
        y_test = np.arange(1)
        scorer = make_scorer(mean_absolute_error)
        return _score_weighted(estimator, x_test[:, features], y_test, scorer)

    def test_alternative_kwargs(self):
        """Test alternative kwargs."""
        rfe_kwargs = {
            'estimator': LinearRegression(),
            'n_features_to_select': None,
            'step': 0.1,
        }
        rfe = AdvancedRFE(**rfe_kwargs)
        zero_idx = np.array([1, 3, 5, 7, 9, 11, 13, 15, 17, 19])
        x_data = np.arange(3 * 20).reshape(3, 20)
        x_data[:, zero_idx] = 0.0
        y_data = np.arange(3)

        rfe._fit(x_data, y_data, step_score=self.step_score)

        assert rfe.n_features_ == 10
        assert len(rfe.ranking_) == 20
        assert len(rfe.support_) == 20
        expected_support = np.full(20, True)
        expected_support[zero_idx] = False
        np.testing.assert_array_equal(rfe.support_, expected_support)
        assert len(rfe.scores_) == 6


# AdvancedRFECV


class TestAdvancedRFECV():
    """Tests for ``AdvancedRFECV``."""

    @pytest.fixture
    def lin(self):
        """Return ``LinearRegression`` instance."""
        return LinearRegression()

    def test_init(self, lin):
        """Test ``__init__``."""
        rfecv = AdvancedRFECV(estimator=lin, step=2, min_features_to_select=3,
                              cv=5, scoring='neg_mean_absolute_error',
                              verbose=42, n_jobs=32)
        assert rfecv.estimator is lin
        assert rfecv.step == 2
        assert rfecv.min_features_to_select == 3
        assert rfecv.cv == 5
        assert rfecv.scoring == 'neg_mean_absolute_error'
        assert rfecv.verbose == 42
        assert rfecv.n_jobs == 32

    X_DATA = np.array([
        [0, 0, 0],
        [1, 1, 0],
        [2, 0, 0],
        [0, 3, 0],
        [0, 3, 0],
        [4, 4, 0],
        [4, 4, 0],
        [1000.0, 2000.0, 0.0],
    ])
    Y_DATA = np.array([1, 0, 3, -5, -5, -3, -3, -4])
    SAMPLE_WEIGHTS = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0])

    def test_fail(self, lin):
        """Test ``AdvancedRFECV`` expected fail."""
        msg = 'Step must be >0'
        rfecv = AdvancedRFECV(estimator=lin, step=-1)
        with pytest.raises(ValueError, match=msg):
            rfecv.fit(self.X_DATA, self.Y_DATA)

    def test_fit(self, lin):
        """Test ``fit``."""
        rfecv = AdvancedRFECV(estimator=lin, step=1, min_features_to_select=1,
                              cv=2, verbose=1000, n_jobs=2)
        rfecv.fit(self.X_DATA, self.Y_DATA, sample_weight=self.SAMPLE_WEIGHTS)
        assert rfecv.n_features_ == 2
        np.testing.assert_array_equal(rfecv.support_, [True, True, False])
        np.testing.assert_array_equal(rfecv.ranking_, [1, 1, 2])
        np.testing.assert_allclose(rfecv.grid_scores_,
                                   [-7.28912807, -0.69779194, -0.69779194])

        est = rfecv.estimator_
        assert isinstance(est, LinearRegression)
        np.testing.assert_allclose(est.coef_, [1.0, -2.0])
        np.testing.assert_allclose(est.intercept_, [1.0])

    def test_step_float(self, lin):
        """Test float for ``step``."""
        rfecv = AdvancedRFECV(estimator=lin, step=0.1, cv=2)
        rfecv.fit(self.X_DATA, self.Y_DATA)

        assert rfecv.n_features_ == 2
        np.testing.assert_array_equal(rfecv.support_, [True, True, False])
        np.testing.assert_array_equal(rfecv.ranking_, [1, 1, 2])
        np.testing.assert_allclose(
            rfecv.grid_scores_,
            [-3949286.19763361, -1630913.74908173, -1630913.74908173])

        est = rfecv.estimator_
        assert isinstance(est, LinearRegression)
        np.testing.assert_allclose(est.coef_, [0.99952835, -0.5006662])
        np.testing.assert_allclose(est.intercept_, -2.21009561525743)


# AdvancedTransformedTargetRegressor


class TestAdvancedTransformedTargetRegressor():
    """Tests for ``AdvancedTransformedTargetRegressor``."""

    def test_regressor_none(self):
        """Test ``regressor=None``."""
        areg = AdvancedTransformedTargetRegressor()
        areg.fit(np.arange(3).reshape(3, 1), np.arange(3))
        assert isinstance(areg.regressor_, LinearRegression)

    def test_coef_(self):
        """Test ``coef_`` property."""
        areg = AdvancedTransformedTargetRegressor()
        areg.fit(np.arange(3).reshape(3, 1), np.arange(3))
        np.testing.assert_allclose(areg.coef_, [1.0])

    def test_feature_importances_(self):
        """Test ``feature_importances_`` property."""
        areg = AdvancedTransformedTargetRegressor(
            regressor=FeatureImportanceRegressor())
        areg.fit(np.arange(3).reshape(3, 1), np.arange(3))
        assert areg.feature_importances_ == 42

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

    def identity(self, y_data):
        """Identity function."""
        return y_data

    def square(self, y_data):
        """Identity function."""
        return y_data**2

    def test_fit_transformer_fail(self):
        """Test ``_fit_transformer`` expected fail."""
        # Give transformer and func/inverse_func
        msg = ("'transformer' and functions 'func'/'inverse_func' cannot both "
               "be set.")
        areg = AdvancedTransformedTargetRegressor(
            transformer=StandardScaler(),
            func=self.identity,
            inverse_func=self.identity,
        )
        with pytest.raises(ValueError, match=msg):
            areg._fit_transformer(self.Y_2D)

        # Give func without inverse_func
        msg = "When 'func' is provided, 'inverse_func' must also be provided"
        areg = AdvancedTransformedTargetRegressor(
            func=self.identity,
            inverse_func=None,
        )
        with pytest.raises(ValueError, match=msg):
            areg._fit_transformer(self.Y_2D)

        # Warn if inverse_func is not true inverse of func
        msg = ("The provided functions or transformer are not strictly "
               "inverse of each other. If you are sure you want to proceed "
               "regardless, set 'check_inverse=False'")
        areg = AdvancedTransformedTargetRegressor(
            func=self.identity,
            inverse_func=self.square,
            check_inverse=True,
        )
        with pytest.warns(UserWarning, match=msg):
            areg._fit_transformer(self.Y_2D)

        # Do not warn if not specified
        areg = AdvancedTransformedTargetRegressor(
            func=self.identity,
            inverse_func=self.square,
            check_inverse=False,
        )
        with pytest.warns(None) as record:
            areg._fit_transformer(self.Y_2D)
        assert not record

    def test_fit_transformer_transformer(self):
        """Test ``_fit_transformer`` with transformer."""
        areg = AdvancedTransformedTargetRegressor(
            transformer=StandardScaler(),
        )
        areg._fit_transformer(self.Y_2D)
        assert isinstance(areg.transformer_, StandardScaler)
        np.testing.assert_allclose(areg.transformer_.scale_, [8.16496581])
        np.testing.assert_allclose(areg.transformer_.mean_, [20.0])

    def test_fit_transformer_func(self):
        """Test ``_fit_transformer`` with func."""
        areg = AdvancedTransformedTargetRegressor(
            func=self.identity,
            inverse_func=self.square,
            check_inverse=False,
        )
        areg._fit_transformer(self.Y_2D)
        assert isinstance(areg.transformer_, FunctionTransformer)
        np.testing.assert_allclose(
            areg.transformer_.transform([[42.0]]), [[42.0]])
        np.testing.assert_allclose(
            areg.transformer_.inverse_transform([[42.0]]), [[1764.0]])

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

    class Reg2DPrediction(BaseEstimator):
        """Estimator with 2D prediction output."""

        def fit(self, *_):
            """Fit method."""
            return self

        def predict(self, *_):
            """Predict method that returns 2D array."""
            return np.array([[42.0]])

    def test_predict_output_2d(self):
        """Test prediction."""
        areg = AdvancedTransformedTargetRegressor(
            transformer=StandardScaler(),
            regressor=self.Reg2DPrediction(),
        )
        areg.fit(np.arange(3).reshape(3, 1), np.arange(3))
        pred = areg.predict([[1]])
        np.testing.assert_allclose(areg.transformer_.scale_, [0.8164965809])
        np.testing.assert_allclose(areg.transformer_.mean_, [1.0])
        np.testing.assert_allclose(pred, [42.0 * 0.8164965809 + 1.0])

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


# FeatureSelectionTransformer


class TestFeatureSelectionTransformer():
    """Tests for ``FeatureSelectionTransformer``."""

    @pytest.fixture
    def fst(self):
        """Return ``FeatureSelectionTransformer`` instance."""
        return FeatureSelectionTransformer(grid_scores=1, n_features=2,
                                           ranking=3, support=4)

    def test_init(self, fst):
        """Test ``__init__``."""
        assert fst.grid_scores == 1
        assert fst.n_features == 2
        assert fst.ranking == 3
        assert fst.support == 4

    def test_fit(self, fst):
        """Test ``fit``."""
        output = fst.fit()
        assert output is fst
        output = fst.fit(1, 'a', valid_kwarg=2)
        assert output is fst

    def test_get_support_mask(self, fst):
        """Test ``_get_support_mask``."""
        mask = fst._get_support_mask()
        assert mask == 4

    def test_more_tags(self, fst):
        """Test ``_more_tags``."""
        tags = fst._more_tags()
        assert tags['allow_nan'] is True
        assert tags is not _DEFAULT_TAGS
        new_tags = deepcopy(_DEFAULT_TAGS)
        new_tags['allow_nan'] = True
        assert tags == new_tags
