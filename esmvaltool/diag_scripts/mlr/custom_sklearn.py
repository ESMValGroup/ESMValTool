"""Custom expansions of :mod:`sklearn` functionalities.

Note
----
This module provides custom expansions of some :mod:`sklearn` classes and
functions which are necessary to fit the purposes for the desired
functionalities of the :ref:`MLR module <api.esmvaltool.diag_scripts.mlr>`. As
long-term goal we would like to include these functionalities to the
:mod:`sklearn` package since we believe these additions might be helpful for
everyone. This module serves as interim solution. To ensure that all features
are properly working this module is also covered by extensive tests.

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
# pylint: disable=protected-access
# pylint: disable=super-init-not-called
# pylint: disable=too-many-arguments
# pylint: disable=too-many-instance-attributes
# pylint: disable=too-many-locals
# pylint: disable=too-many-return-statements

import itertools
import logging
import numbers
import os
import warnings
from contextlib import suppress
from copy import deepcopy
from inspect import getfullargspec
from traceback import format_exc

import numpy as np
import scipy.sparse as sp
from joblib import Parallel, delayed, effective_n_jobs
from sklearn.base import BaseEstimator, clone, is_classifier
from sklearn.compose import ColumnTransformer, TransformedTargetRegressor
from sklearn.exceptions import FitFailedWarning, NotFittedError
from sklearn.feature_selection import RFE, SelectorMixin
from sklearn.linear_model import LinearRegression
from sklearn.metrics import check_scoring
from sklearn.model_selection import check_cv
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer
from sklearn.utils import check_array, check_X_y, indexable, safe_sqr
from sklearn.utils.fixes import np_version, parse_version
from sklearn.utils.metaestimators import if_delegate_has_method
from sklearn.utils.validation import check_is_fitted

from esmvaltool.diag_scripts import mlr

logger = logging.getLogger(os.path.basename(__file__))


_DEFAULT_TAGS = {
    'non_deterministic': False,
    'requires_positive_X': False,
    'requires_positive_y': False,
    'X_types': ['2darray'],
    'poor_score': False,
    'no_validation': False,
    'multioutput': False,
    "allow_nan": False,
    'stateless': False,
    'multilabel': False,
    '_skip_test': False,
    '_xfail_checks': False,
    'multioutput_only': False,
    'binary_only': False,
    'requires_fit': True,
    'preserves_dtype': [np.float64],
    'requires_y': False,
    'pairwise': False,
}


def _determine_key_type(key, accept_slice=True):
    """Determine the data type of key."""
    err_msg = ("No valid specification of the columns. Only a scalar, list or "
               "slice of all integers or all strings, or boolean mask is "
               "allowed")

    dtype_to_str = {int: 'int', str: 'str', bool: 'bool', np.bool_: 'bool'}
    array_dtype_to_str = {'i': 'int', 'u': 'int', 'b': 'bool', 'O': 'str',
                          'U': 'str', 'S': 'str'}

    if key is None:
        return None
    if isinstance(key, tuple(dtype_to_str.keys())):
        try:
            return dtype_to_str[type(key)]
        except KeyError:
            raise ValueError(err_msg)
    if isinstance(key, slice):
        if not accept_slice:
            raise TypeError(
                'Only array-like or scalar are supported. '
                'A Python slice was given.'
            )
        if key.start is None and key.stop is None:
            return None
        key_start_type = _determine_key_type(key.start)
        key_stop_type = _determine_key_type(key.stop)
        if key_start_type is not None and key_stop_type is not None:
            if key_start_type != key_stop_type:
                raise ValueError(err_msg)
        if key_start_type is not None:
            return key_start_type
        return key_stop_type
    if isinstance(key, (list, tuple)):
        unique_key = set(key)
        key_type = {_determine_key_type(elt) for elt in unique_key}
        if not key_type:
            return None
        if len(key_type) != 1:
            raise ValueError(err_msg)
        return key_type.pop()
    if hasattr(key, 'dtype'):
        try:
            return array_dtype_to_str[key.dtype.kind]
        except KeyError:
            raise ValueError(err_msg)
    raise ValueError(err_msg)


def _array_indexing(array, key, key_dtype, axis):
    """Index an array or scipy.sparse consistently across numpy version."""
    if np_version < parse_version('1.12') or sp.issparse(array):
        if key_dtype == 'bool':
            key = np.asarray(key)
    if isinstance(key, tuple):
        key = list(key)
    return array[key] if axis == 0 else array[:, key]


def _list_indexing(x_data, key, key_dtype):
    """Index a python list."""
    if np.isscalar(key) or isinstance(key, slice):
        # key is a slice or a scalar
        return x_data[key]
    if key_dtype == 'bool':
        # key is a boolean array-like
        return list(itertools.compress(x_data, key))
    # key is a integer array-like of key
    return [x_data[idx] for idx in key]


def _pandas_indexing(x_data, key, key_dtype, axis):
    """Index a pandas dataframe or a series."""
    if hasattr(key, 'shape'):
        key = np.asarray(key)
        key = key if key.flags.writeable else key.copy()
    elif isinstance(key, tuple):
        key = list(key)
    # check whether we should index with loc or iloc
    indexer = x_data.iloc if key_dtype == 'int' else x_data.loc
    return indexer[:, key] if axis else indexer[key]


def _safe_indexing(x_data, indices, *_, axis=0):
    """Return rows, items or columns of x_data using indices."""
    if indices is None:
        return x_data

    if axis not in (0, 1):
        raise ValueError(
            "'axis' should be either 0 (to index rows) or 1 (to index "
            "column). Got {} instead.".format(axis)
        )

    indices_dtype = _determine_key_type(indices)

    if axis == 0 and indices_dtype == 'str':
        raise ValueError(
            "String indexing is not supported with 'axis=0'"
        )

    if axis == 1 and x_data.ndim != 2:
        raise ValueError(
            "'x_data' should be a 2D NumPy array, 2D sparse matrix or pandas "
            "dataframe when indexing the columns (i.e. 'axis=1'). "
            "Got {} instead with {} dimension(s).".format(type(x_data),
                                                          x_data.ndim)
        )

    if axis == 1 and indices_dtype == 'str' and not hasattr(x_data, 'loc'):
        raise ValueError(
            "Specifying the columns using strings is only supported for "
            "pandas DataFrames"
        )

    if hasattr(x_data, "iloc"):
        return _pandas_indexing(x_data, indices, indices_dtype, axis=axis)
    if hasattr(x_data, "shape"):
        return _array_indexing(x_data, indices, indices_dtype, axis=axis)
    return _list_indexing(x_data, indices, indices_dtype)


def _is_arraylike(input_array):
    """Check whether the input is array-like."""
    return (hasattr(input_array, '__len__') or
            hasattr(input_array, 'shape') or
            hasattr(input_array, '__array__'))


def _make_indexable(iterable):
    """Ensure iterable supports indexing or convert to an indexable variant."""
    if sp.issparse(iterable):
        return iterable.tocsr()
    if hasattr(iterable, "__getitem__") or hasattr(iterable, "iloc"):
        return iterable
    if iterable is None:
        return iterable
    return np.array(iterable)


def _num_samples(x_data):
    """Return number of samples in array-like x_data."""
    message = 'Expected sequence or array-like, got %s' % type(x_data)
    if hasattr(x_data, 'fit') and callable(x_data.fit):
        # Don't get num_samples from an ensembles length!
        raise TypeError(message)

    if not hasattr(x_data, '__len__') and not hasattr(x_data, 'shape'):
        if hasattr(x_data, '__array__'):
            x_data = np.asarray(x_data)
        else:
            raise TypeError(message)

    if hasattr(x_data, 'shape') and x_data.shape is not None:
        if len(x_data.shape) == 0:
            raise TypeError("Singleton array %r cannot be considered a valid "
                            "collection." % x_data)
        # Check that shape is returning an integer or default to len
        # Dask dataframes may not return numeric shape[0] value
        if isinstance(x_data.shape[0], numbers.Integral):
            return x_data.shape[0]

    try:
        return len(x_data)
    except TypeError as type_error:
        raise TypeError(message) from type_error


def _check_fit_params(x_data, fit_params, indices=None):
    """Check and validate the parameters passed during ``fit``."""
    fit_params_validated = {}
    for param_key, param_value in fit_params.items():
        if (not _is_arraylike(param_value) or
                _num_samples(param_value) != _num_samples(x_data)):
            # Non-indexable pass-through (for now for backward-compatibility).
            # https://github.com/scikit-learn/scikit-learn/issues/15805
            fit_params_validated[param_key] = param_value
        else:
            # Any other fit_params should support indexing
            # (e.g. for cross-validation).
            fit_params_validated[param_key] = _make_indexable(param_value)
            fit_params_validated[param_key] = _safe_indexing(
                fit_params_validated[param_key], indices
            )

    return fit_params_validated


def _safe_tags(estimator, key=None):
    """Safely get estimator tags."""
    if hasattr(estimator, "_get_tags"):
        tags_provider = "_get_tags()"
        tags = estimator._get_tags()
    elif hasattr(estimator, "_more_tags"):
        tags_provider = "_more_tags()"
        tags = {**_DEFAULT_TAGS, **estimator._more_tags()}
    else:
        tags_provider = "_DEFAULT_TAGS"
        tags = _DEFAULT_TAGS

    if key is not None:
        if key not in tags:
            raise ValueError(
                f"The key {key} is not defined in {tags_provider} for the "
                f"class {estimator.__class__.__name__}."
            )
        return tags[key]
    return tags


def _is_pairwise(estimator):
    """Return ``True`` if estimator is pairwise."""
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=FutureWarning)
        has_pairwise_attribute = hasattr(estimator, '_pairwise')
        pairwise_attribute = getattr(estimator, '_pairwise', False)
    pairwise_tag = _safe_tags(estimator, key="pairwise")

    if has_pairwise_attribute:
        if pairwise_attribute != pairwise_tag:
            warnings.warn(
                "_pairwise attribute is inconsistent with tags. Set the "
                "estimator tags of your estimator instead", FutureWarning,
            )
        return pairwise_attribute

    # Use pairwise tag when the attribute is not present
    return pairwise_tag


def _safe_split(estimator, x_data, y_data, indices, train_indices=None):
    """Create subset of dataset and properly handle kernels."""
    if _is_pairwise(estimator):
        if not hasattr(x_data, "shape"):
            raise ValueError("Precomputed kernels or affinity matrices have "
                             "to be passed as arrays or sparse matrices.")
        # x_data is a precomputed square kernel matrix
        if x_data.shape[0] != x_data.shape[1]:
            raise ValueError("x_data should be a square kernel matrix")
        if train_indices is None:
            x_subset = x_data[np.ix_(indices, indices)]
        else:
            x_subset = x_data[np.ix_(indices, train_indices)]
    else:
        x_subset = _safe_indexing(x_data, indices)

    if y_data is not None:
        y_subset = _safe_indexing(y_data, indices)
    else:
        y_subset = None

    return (x_subset, y_subset)


def _fit_and_score_weighted(estimator, x_data, y_data, scorer, train, test,
                            parameters, fit_params, error_score=np.nan,
                            sample_weights=None):
    """Expand :func:`sklearn.model_selection._validation._fit_and_score`."""
    # Adjust length of sample weights
    fit_params = fit_params if fit_params is not None else {}
    fit_params = _check_fit_params(x_data, fit_params, train)

    if parameters is not None:
        # clone after setting parameters in case any parameters
        # are estimators (like pipeline steps)
        # because pipeline doesn't clone steps in fit
        cloned_parameters = {}
        for (key, val) in parameters.items():
            cloned_parameters[key] = clone(val, safe=False)

        estimator = estimator.set_params(**cloned_parameters)

    (x_train, y_train) = _safe_split(estimator, x_data, y_data, train)
    (x_test, y_test) = _safe_split(estimator, x_data, y_data, test, train)
    if sample_weights is not None:
        sample_weights_test = sample_weights[test]
    else:
        sample_weights_test = None

    try:
        if y_train is None:
            estimator.fit(x_train, **fit_params)
        else:
            estimator.fit(x_train, y_train, **fit_params)
    except Exception:
        if error_score == 'raise':
            raise
        if isinstance(error_score, numbers.Number):
            test_score = error_score
            warnings.warn("Estimator fit failed. The score on this train-test "
                          "partition for these parameters will be set to %f. "
                          "Details: \n%s" % (error_score, format_exc()),
                          FitFailedWarning)
        else:
            raise ValueError("error_score must be the string 'raise' or a "
                             "numeric value. (Hint: if using 'raise', please "
                             "make sure that it has been spelled correctly.)")
    else:
        test_score = _score_weighted(estimator, x_test, y_test, scorer,
                                     sample_weights=sample_weights_test)

    return test_score


def _get_fit_parameters(fit_kwargs, steps, cls):
    """Retrieve fit parameters from ``fit_kwargs``."""
    params = {name: {} for (name, step) in steps if step is not None}
    step_names = list(params.keys())
    for (param_name, param_val) in fit_kwargs.items():
        param_split = param_name.split('__', 1)
        if len(param_split) != 2:
            raise ValueError(
                f"Fit parameters for {cls} have to be given in the form "
                f"'s__p', where 's' is the name of the step and 'p' the name "
                f"of the parameter, got '{param_name}'")
        try:
            params[param_split[0]][param_split[1]] = param_val
        except KeyError:
            raise ValueError(
                f"Expected one of {step_names} for step of fit parameter, got "
                f"'{param_split[0]}' for parameter '{param_name}'")
    return params


def _score_weighted(estimator, x_test, y_test, scorer, sample_weights=None):
    """Expand :func:`sklearn.model_selection._validation._score`."""
    if y_test is None:
        score = scorer(estimator, x_test, sample_weight=sample_weights)
    else:
        score = scorer(estimator, x_test, y_test, sample_weight=sample_weights)

    error_msg = ("Scoring must return a number, got %s (%s) instead. "
                 "(scorer=%s)")
    if hasattr(score, 'item'):
        with suppress(ValueError):
            # e.g. unwrap memmapped scalars
            score = score.item()
    if not isinstance(score, numbers.Number):
        raise ValueError(error_msg % (score, type(score), scorer))
    return score


def _split_fit_kwargs(fit_kwargs, train_idx, test_idx):
    """Get split fit kwargs for single CV step."""
    fit_kwargs_train = {}
    fit_kwargs_test = {}
    for (key, val) in fit_kwargs.items():
        if 'sample_weight' in key and 'sample_weight_eval_set' not in key:
            fit_kwargs_train[key] = deepcopy(val)[train_idx]
            fit_kwargs_test[key] = deepcopy(val)[test_idx]
        else:
            fit_kwargs_train[key] = deepcopy(val)
            fit_kwargs_test[key] = deepcopy(val)
    return (fit_kwargs_train, fit_kwargs_test)


def _rfe_single_fit(rfe, estimator, x_data, y_data, train, test, scorer,
                    **fit_kwargs):
    """Return the score for a fit across one fold."""
    (x_train, y_train) = _safe_split(estimator, x_data, y_data, train)
    (x_test, y_test) = _safe_split(estimator, x_data, y_data, test, train)
    (fit_kwargs_train, fit_kwargs_test) = _split_fit_kwargs(fit_kwargs, train,
                                                            test)
    if 'sample_weight' in fit_kwargs_test:
        fit_kwargs_test['sample_weights'] = fit_kwargs_test.pop(
            'sample_weight')

    def step_score(estimator, features):
        """Score for a single step in the recursive feature elimination."""
        return _score_weighted(estimator, x_test[:, features], y_test, scorer,
                               **fit_kwargs_test)

    return rfe._fit(x_train, y_train, step_score=step_score,
                    **fit_kwargs_train).scores_


def _map_features(features, support):
    """Map old features indices to new ones using boolean mask."""
    feature_mapping = {}
    new_idx = 0
    for (old_idx, supported) in enumerate(support):
        if supported:
            val = new_idx
            new_idx += 1
        else:
            val = None
        feature_mapping[old_idx] = val
    new_features = []
    for feature in features:
        new_feature = feature_mapping[feature]
        if new_feature is not None:
            new_features.append(new_feature)
    return new_features


def _update_transformers_param(estimator, support):
    """Update ``transformers`` argument of ``ColumnTransformer`` steps."""
    all_params = estimator.get_params()
    params = []
    for key in all_params:
        if key.endswith('transformers'):
            params.append(key)
            if isinstance(estimator, (Pipeline, AdvancedPipeline)):
                step = estimator.named_steps[key.split('__')[0]]
                if not isinstance(step, ColumnTransformer):
                    raise TypeError(
                        f"Found 'transformers' parameter ('{key}'), but the "
                        f"corresponding pipeline step is not a "
                        f"ColumnTransformer (got '{type(step)}')")
            else:
                raise TypeError(
                    f"Found 'transformers' parameter ('{key}'), but the "
                    f"corresponding estimator is not a Pipeline or "
                    f"AdvancedPipeline")
    new_params = {}
    for param in params:
        new_transformers = []
        for transformer in all_params[param]:
            new_columns = _map_features(transformer[2], support)
            new_transformers.append(
                (transformer[0], transformer[1], new_columns))
        new_params[param] = new_transformers
    estimator.set_params(**new_params)


def cross_val_score_weighted(estimator, x_data, y_data=None, groups=None,
                             scoring=None, cv=None, n_jobs=None, verbose=0,
                             fit_params=None, pre_dispatch='2*n_jobs',
                             error_score=np.nan, sample_weights=None):
    """Expand :func:`sklearn.model_selection.cross_val_score`."""
    scorer = check_scoring(estimator, scoring=scoring)
    (x_data, y_data, groups) = indexable(x_data, y_data, groups)

    cv = check_cv(cv, y_data, classifier=is_classifier(estimator))

    # We clone the estimator to make sure that all the folds are
    # independent, and that it is pickle-able.
    parallel = Parallel(n_jobs=n_jobs, verbose=verbose,
                        pre_dispatch=pre_dispatch)
    scores = parallel(
        delayed(_fit_and_score_weighted)(
            clone(estimator), x_data, y_data, scorer, train, test, None,
            fit_params, error_score=error_score, sample_weights=sample_weights)
        for train, test in cv.split(x_data, y_data, groups))
    return np.array(scores)


def get_rfecv_transformer(rfecv_estimator):
    """Get transformer step of RFECV estimator."""
    try:
        check_is_fitted(rfecv_estimator)
    except NotFittedError:
        raise NotFittedError(
            "RFECV instance used to initialize FeatureSelectionTransformer "
            "must be fitted")
    transformer = FeatureSelectionTransformer(
        grid_scores=rfecv_estimator.grid_scores_,
        n_features=rfecv_estimator.n_features_,
        ranking=rfecv_estimator.ranking_,
        support=rfecv_estimator.support_,
    )
    return transformer


def perform_efecv(estimator, x_data, y_data, **kwargs):
    """Perform exhaustive feature selection."""
    x_data, y_data = check_X_y(
        x_data, y_data, ensure_min_features=2, force_all_finite='allow-nan')
    n_all_features = x_data.shape[1]

    # Iterate over all possible feature combinations
    supports = list(itertools.product([False, True], repeat=n_all_features))
    supports.remove(tuple([False] * n_all_features))
    logger.info(
        "Testing all %i possible feature combinations for exhaustive feature "
        "selection", len(supports))
    grid_scores = []
    for support in supports:
        support = np.array(support)
        features = np.arange(n_all_features)[support]

        # Evaluate estimator on new subset of features
        new_estimator = clone(estimator)
        _update_transformers_param(new_estimator, support)
        scores = cross_val_score_weighted(new_estimator, x_data[:, features],
                                          y_data, **kwargs)
        grid_scores.append(np.mean(scores))
        logger.debug("Fitted estimator with %i features, CV score was %.5f",
                     support.sum(), np.mean(scores))

    # Final parameters
    grid_scores = np.array(grid_scores)
    best_idx = np.argmax(grid_scores)
    support = np.array(supports[best_idx])
    features = np.arange(n_all_features)[support]
    n_features = support.sum()
    ranking = np.where(support, 1, 2)
    transformer = FeatureSelectionTransformer(
        grid_scores=grid_scores, n_features=n_features, ranking=ranking,
        support=support)

    # Get final estimator
    best_estimator = clone(estimator)
    _update_transformers_param(best_estimator, support)
    best_estimator.fit(x_data[:, features], y_data,
                       **kwargs.get('fit_params', {}))

    logger.info("Found optimal score %.5f for %i features",
                grid_scores[best_idx], n_features)
    return (best_estimator, transformer)


class AdvancedPipeline(Pipeline):
    """Expand :class:`sklearn.pipeline.Pipeline`."""

    @property
    def coef_(self):
        """numpy.ndarray: Model coefficients."""
        return self.steps[-1][1].coef_

    @property
    def feature_importances_(self):
        """numpy.ndarray: Feature importances."""
        return self.steps[-1][1].feature_importances_

    def _check_final_step(self):
        """Check type of final step of pipeline."""
        final_step = self.steps[-1][1]
        if not isinstance(final_step, AdvancedTransformedTargetRegressor):
            raise TypeError(
                f"Expected estimator of type "
                f"{AdvancedTransformedTargetRegressor} for final step of "
                f"pipeline, got {final_step.__class__}")

    def fit_target_transformer_only(self, y_data, **fit_kwargs):
        """Fit only ``transform`` step of of target regressor."""
        self._check_final_step()
        reg = self.steps[-1][1]
        fit_params = _get_fit_parameters(fit_kwargs, self.steps,
                                         self.__class__)
        reg_fit_params = fit_params[self.steps[-1][0]]
        reg.fit_transformer_only(y_data, **reg_fit_params)

    def fit_transformers_only(self, x_data, y_data, **fit_kwargs):
        """Fit only ``transform`` steps of Pipeline."""
        fit_params = _get_fit_parameters(fit_kwargs, self.steps,
                                         self.__class__)
        return self._fit(x_data, y_data, **fit_params)

    def transform_only(self, x_data):
        """Only perform ``transform`` steps of Pipeline."""
        for (_, transformer) in self.steps[:-1]:
            x_data = transformer.transform(x_data)
        return x_data

    def transform_target_only(self, y_data):
        """Only perform ``transform`` steps of target regressor."""
        self._check_final_step()
        reg = self.steps[-1][1]
        if not hasattr(reg, 'transformer_'):
            raise NotFittedError(
                "Transforming target not possible, final regressor is not "
                "fitted yet, call fit() or fit_target_transformer_only() "
                "first")
        if y_data.ndim == 1:
            y_data = y_data.reshape(-1, 1)
        y_trans = reg.transformer_.transform(y_data)
        if y_trans.ndim == 2 and y_trans.shape[1] == 1:
            y_trans = y_trans.squeeze(axis=1)
        return y_trans


class AdvancedRFE(RFE):
    """Expand :class:`sklearn.feature_selection.RFE`."""

    def fit(self, x_data, y_data, **fit_kwargs):
        """Expand :meth:`fit` to accept kwargs."""
        return self._fit(x_data, y_data, **fit_kwargs)

    def _fit(self, x_data, y_data, step_score=None, **fit_kwargs):
        """Expand :meth:`_fit` to accept kwargs."""
        # Parameter step_score controls the calculation of self.scores_
        # step_score is not exposed to users
        # and is used when implementing AdvancedRFECV
        # self.scores_ will not be calculated when calling _fit through fit
        x_data, y_data = check_X_y(x_data, y_data, "csc",
                                   ensure_min_features=2,
                                   force_all_finite=False)

        # Initialization
        n_features = x_data.shape[1]
        if self.n_features_to_select is None:
            n_features_to_select = n_features // 2
        else:
            n_features_to_select = self.n_features_to_select

        if 0.0 < self.step < 1.0:
            step = int(max(1, self.step * n_features))
        else:
            step = int(self.step)
        if step <= 0:
            raise ValueError("Step must be >0")

        support_ = np.ones(n_features, dtype=bool)
        ranking_ = np.ones(n_features, dtype=int)

        if step_score:
            self.scores_ = []

        # Elimination
        while np.sum(support_) > n_features_to_select:
            # Remaining features
            features = np.arange(n_features)[support_]

            # Rank the remaining features
            estimator = clone(self.estimator)
            if self.verbose > 0:
                print("Fitting estimator with %d features." % np.sum(support_))

            _update_transformers_param(estimator, support_)
            estimator.fit(x_data[:, features], y_data, **fit_kwargs)

            # Get coefs (hasattr(estimator, 'coef_') raises a KeyError for
            # XGBRegressor models
            try:
                coefs = estimator.coef_
            except (AttributeError, KeyError):
                coefs = getattr(estimator, 'feature_importances_', None)
            if coefs is None:
                raise RuntimeError("The classifier does not expose "
                                   "'coef_' or 'feature_importances_' "
                                   "attributes")

            # Get ranks
            if coefs.ndim > 1:
                ranks = np.argsort(safe_sqr(coefs).sum(axis=0))
            else:
                ranks = np.argsort(safe_sqr(coefs))

            # Transformer steps that reduce number of features is not supported
            if len(ranks) != len(features):
                raise NotImplementedError(
                    f"Estimators that contain transforming steps that reduce "
                    f"the number of features are not supported in "
                    f"{self.__class__}, got {len(features):d} features for ",
                    f"fit(), but only {len(ranks):d} elements for 'coefs_' / "
                    f"'feature_importances_' are provided. Estimator:\n"
                    f"{estimator}")

            # for sparse case ranks is matrix
            ranks = np.ravel(ranks)

            # Eliminate the worse features
            threshold = min(step, np.sum(support_) - n_features_to_select)

            # Compute step score on the previous selection iteration
            # because 'estimator' must use features
            # that have not been eliminated yet
            if step_score:
                self.scores_.append(step_score(estimator, features))
            support_[features[ranks][:threshold]] = False
            ranking_[np.logical_not(support_)] += 1

        # Set final attributes
        features = np.arange(n_features)[support_]
        self.estimator_ = clone(self.estimator)
        _update_transformers_param(self.estimator_, support_)
        self.estimator_.fit(x_data[:, features], y_data, **fit_kwargs)

        # Compute step score when only n_features_to_select features left
        if step_score:
            self.scores_.append(step_score(self.estimator_, features))
        self.n_features_ = support_.sum()
        self.support_ = support_
        self.ranking_ = ranking_

        return self

    @if_delegate_has_method(delegate='estimator')
    def predict(self, x_data, **predict_kwargs):
        """Expand :meth:`predict()` to accept kwargs."""
        check_is_fitted(self)
        return self.estimator_.predict(self.transform(x_data),
                                       **predict_kwargs)


class AdvancedRFECV(AdvancedRFE):
    """Expand :class:`sklearn.feature_selection.RFECV`."""

    def __init__(self, estimator, step=1, min_features_to_select=1, cv=None,
                 scoring=None, verbose=0, n_jobs=None):
        """Original constructor of :class:`sklearn.feature_selection.RFECV`."""
        self.estimator = estimator
        self.step = step
        self.min_features_to_select = min_features_to_select
        self.cv = cv
        self.scoring = scoring
        self.verbose = verbose
        self.n_jobs = n_jobs

    def fit(self, x_data, y_data, groups=None, **fit_kwargs):
        """Expand :meth:`fit` to accept kwargs."""
        x_data, y_data = check_X_y(
            x_data, y_data, "csr", ensure_min_features=2,
            force_all_finite=False)

        # Initialization
        cv = check_cv(self.cv, y_data,
                      classifier=is_classifier(self.estimator))
        scorer = check_scoring(self.estimator, scoring=self.scoring)
        n_features = x_data.shape[1]

        if 0.0 < self.step < 1.0:
            step = int(max(1, self.step * n_features))
        else:
            step = int(self.step)
        if step <= 0:
            raise ValueError("Step must be >0")

        # Build an AdvancedRFE object, which will evaluate and score each
        # possible feature count, down to self.min_features_to_select
        rfe = AdvancedRFE(estimator=self.estimator,
                          n_features_to_select=self.min_features_to_select,
                          step=self.step, verbose=self.verbose)

        # Determine the number of subsets of features by fitting across
        # the train folds and choosing the "features_to_select" parameter
        # that gives the least averaged error across all folds.

        # Note that joblib raises a non-picklable error for bound methods
        # even if n_jobs is set to 1 with the default multiprocessing
        # backend.
        # This branching is done so that to
        # make sure that user code that sets n_jobs to 1
        # and provides bound methods as scorers is not broken with the
        # addition of n_jobs parameter.

        if effective_n_jobs(self.n_jobs) == 1:
            (parallel, func) = (list, _rfe_single_fit)
        else:
            parallel = Parallel(n_jobs=self.n_jobs)
            func = delayed(_rfe_single_fit)

        scores = parallel(
            func(rfe, self.estimator, x_data, y_data, train, test, scorer,
                 **fit_kwargs)
            for train, test in cv.split(x_data, y_data, groups))

        scores = np.sum(scores, axis=0)
        scores_rev = scores[::-1]
        argmax_idx = len(scores) - np.argmax(scores_rev) - 1
        n_features_to_select = max(
            n_features - (argmax_idx * step),
            self.min_features_to_select)

        # Re-execute an elimination with best_k over the whole set
        rfe = AdvancedRFE(estimator=self.estimator,
                          n_features_to_select=n_features_to_select,
                          step=self.step, verbose=self.verbose)

        rfe.fit(x_data, y_data, **fit_kwargs)

        # Set final attributes
        self.support_ = rfe.support_
        self.n_features_ = rfe.n_features_
        self.ranking_ = rfe.ranking_
        self.estimator_ = clone(self.estimator)
        _update_transformers_param(self.estimator_, self.support_)
        self.estimator_.fit(self.transform(x_data), y_data, **fit_kwargs)

        # Fixing a normalization error, n is equal to
        # get_n_splits(x_data, y_data) - 1 here, the scores are normalized by
        # get_n_splits(x_data, y_data)
        self.grid_scores_ = scores[::-1] / cv.get_n_splits(x_data, y_data,
                                                           groups)
        return self


class AdvancedTransformedTargetRegressor(TransformedTargetRegressor):
    """Expand :class:`sklearn.compose.TransformedTargetRegressor`."""

    @property
    def coef_(self):
        """numpy.ndarray: Model coefficients."""
        return self.regressor_.coef_

    @property
    def feature_importances_(self):
        """numpy.ndarray: Feature importances."""
        return self.regressor_.feature_importances_

    def fit(self, x_data, y_data, **fit_kwargs):
        """Expand :meth:`fit` to accept kwargs."""
        (y_2d,
         regressor_kwargs) = self.fit_transformer_only(y_data, **fit_kwargs)

        # Transform y and convert back to 1d array if necessary
        y_trans = self.transformer_.transform(y_2d)
        if y_trans.ndim == 2 and y_trans.shape[1] == 1:
            y_trans = y_trans.squeeze(axis=1)

        # Perform linear regression if regressor is not given
        if self.regressor is None:
            self.regressor_ = LinearRegression()
        else:
            self.regressor_ = clone(self.regressor)

        # Fit regressor with kwargs
        self.regressor_.fit(x_data, y_trans, **regressor_kwargs)
        return self

    def fit_transformer_only(self, y_data, **fit_kwargs):
        """Fit only ``transformer`` step."""
        y_data = check_array(y_data,
                             accept_sparse=False,
                             force_all_finite=True,
                             ensure_2d=False,
                             dtype='numeric')
        self._training_dim = y_data.ndim

        # Process kwargs
        (_, regressor_kwargs) = self._get_fit_params(fit_kwargs)

        # Transformers are designed to modify X which is 2D, modify y_data
        # FIXME: Transformer does NOT use transformer_kwargs
        if y_data.ndim == 1:
            y_2d = y_data.reshape(-1, 1)
        else:
            y_2d = y_data
        self._fit_transformer(y_2d)
        return (y_2d, regressor_kwargs)

    def predict(self, x_data, always_return_1d=True, **predict_kwargs):
        """Expand :meth:`predict()` to accept kwargs."""
        check_is_fitted(self)
        if not hasattr(self, 'regressor_'):
            raise NotFittedError(
                f"Regressor of {self.__class__} is not fitted yet, call fit() "
                f"first")

        # Kwargs for returning variance or covariance
        if ('return_std' in predict_kwargs and 'return_std' in getfullargspec(
                self.regressor_.predict).args):
            raise NotImplementedError(
                f"Using keyword argument 'return_std' for final regressor "
                f"{self.regressor_.__class__} is not supported yet, only "
                f"'return_var' is allowed. Expand the regressor to accept "
                f"'return_var' instead (see 'esmvaltool/diag_scripts/mlr"
                f"/models/gpr_sklearn.py' for an example)")
        mlr.check_predict_kwargs(predict_kwargs)
        return_var = predict_kwargs.get('return_var', False)
        return_cov = predict_kwargs.get('return_cov', False)

        # Prediction
        prediction = self.regressor_.predict(x_data, **predict_kwargs)
        if return_var or return_cov:
            pred = prediction[0]
        else:
            pred = prediction
        if pred.ndim == 1:
            pred_trans = self.transformer_.inverse_transform(
                pred.reshape(-1, 1))
        else:
            pred_trans = self.transformer_.inverse_transform(pred)
        if self._to_be_squeezed(pred_trans, always_return_1d=always_return_1d):
            pred_trans = pred_trans.squeeze(axis=1)
        if not (return_var or return_cov):
            return pred_trans

        # Return scaled variance or covariance if desired
        err = prediction[1]
        if not hasattr(self.transformer_, 'scale_'):
            raise NotImplementedError(
                f"Transforming of additional prediction output (e.g. by "
                f"'return_var' or 'return_cov') is not supported for "
                f"transformer {self.transformer_.__class__} yet, the "
                f"necessary attribute 'scale_' is missing")
        scale = self.transformer_.scale_
        if scale is not None:
            err *= scale**2
        if self._to_be_squeezed(err, always_return_1d=always_return_1d):
            err = err.squeeze(axis=1)
        return (pred_trans, err)

    def _get_fit_params(self, fit_kwargs):
        """Separate ``transformer`` and ``regressor`` kwargs."""
        steps = [
            ('transformer', self.transformer),
            ('regressor', self.regressor),
        ]
        fit_params = _get_fit_parameters(fit_kwargs, steps, self.__class__)
        fit_params.setdefault('transformer', {})
        fit_params.setdefault('regressor', {})

        # FIXME
        if fit_params['transformer']:
            raise NotImplementedError(
                f"Fit parameters {fit_params['transformer']} for transformer "
                f"{self.transformer.__class__} of {self.__class__} are not "
                f"supported at the moment")

        return (fit_params['transformer'], fit_params['regressor'])

    def _fit_transformer(self, y_data):
        """Check transformer and fit transformer."""
        if (self.transformer is not None and
                (self.func is not None or self.inverse_func is not None)):
            raise ValueError("'transformer' and functions 'func'/"
                             "'inverse_func' cannot both be set.")
        if self.transformer is not None:
            self.transformer_ = clone(self.transformer)
        else:
            if self.func is not None and self.inverse_func is None:
                raise ValueError(
                    "When 'func' is provided, 'inverse_func' must also be "
                    "provided")
            self.transformer_ = FunctionTransformer(
                func=self.func, inverse_func=self.inverse_func, validate=True,
                check_inverse=self.check_inverse)
        self.transformer_.fit(y_data)
        if self.check_inverse:
            idx_selected = slice(None, None, max(1, y_data.shape[0] // 10))
            y_sel = _safe_indexing(y_data, idx_selected)
            y_sel_t = self.transformer_.transform(y_sel)
            if not np.allclose(y_sel,
                               self.transformer_.inverse_transform(y_sel_t)):
                warnings.warn("The provided functions or transformer are "
                              "not strictly inverse of each other. If "
                              "you are sure you want to proceed regardless, "
                              "set 'check_inverse=False'", UserWarning)

    def _to_be_squeezed(self, array, always_return_1d=True):
        """Check if ``array`` should be squeezed or not."""
        squeeze = array.ndim == 2 and array.shape[1] == 1
        if not always_return_1d:
            squeeze = squeeze and self._training_dim == 1
        return squeeze


class FeatureSelectionTransformer(BaseEstimator, SelectorMixin):
    """Transformer step of a feature selection estimator."""

    def __init__(self, grid_scores, n_features, ranking, support):
        """Initialize feature selection transformer."""
        self.grid_scores = grid_scores
        self.n_features = n_features
        self.ranking = ranking
        self.support = support

    def fit(self, *_, **__):
        """Empty method."""
        return self

    def _get_support_mask(self):
        """Get support mask."""
        return self.support

    def _more_tags(self):
        """Additional estimator tags."""
        more_tags = deepcopy(_DEFAULT_TAGS)
        more_tags['allow_nan'] = True
        return more_tags
