"""Custom expansions of :mod:`sklearn` functionalities.

Note
----
This module provides custom expansions of some :mod:`sklearn` classes and
functions which are necessary to fit the purposes for the desired
functionalities of the :ref:`MLR module <api.esmvaltool.diag_scripts.mlr>`. As
long-term goal we would like to include these functionalities to the
:mod:`sklearn` package since we believe these additions might be helpful for
everyone. This module serves as interim solution. To ensure that all features
are properly working this module is also covered by tests, which will also be
expanded in the future.

"""

# pylint: disable=arguments-differ
# pylint: disable=attribute-defined-outside-init
# pylint: disable=protected-access
# pylint: disable=super-init-not-called
# pylint: disable=too-many-arguments
# pylint: disable=too-many-instance-attributes
# pylint: disable=too-many-locals

import itertools
import logging
import numbers
import os
import time
import warnings
from contextlib import suppress
from copy import deepcopy
from inspect import getfullargspec
from traceback import format_exception_only

import numpy as np
from joblib import Parallel, delayed, effective_n_jobs
from sklearn.base import BaseEstimator, clone, is_classifier
from sklearn.compose import ColumnTransformer, TransformedTargetRegressor
from sklearn.exceptions import FitFailedWarning, NotFittedError
from sklearn.feature_selection import RFE
from sklearn.feature_selection._base import SelectorMixin
from sklearn.linear_model import LinearRegression
from sklearn.metrics import check_scoring
from sklearn.metrics._scorer import (
    _check_multimetric_scoring,
    _MultimetricScorer,
)
from sklearn.model_selection import check_cv
from sklearn.model_selection._validation import _aggregate_score_dicts, _score
from sklearn.pipeline import Pipeline
from sklearn.utils import (
    _message_with_time,
    check_array,
    check_X_y,
    indexable,
    safe_sqr,
)
from sklearn.utils.metaestimators import _safe_split, if_delegate_has_method
from sklearn.utils.validation import _check_fit_params, check_is_fitted

from esmvaltool.diag_scripts import mlr

logger = logging.getLogger(os.path.basename(__file__))


def _fit_and_score_weighted(estimator, x_data, y_data, scorer, train, test,
                            verbose, parameters, fit_params,
                            error_score=np.nan, sample_weights=None):
    """Expand :func:`sklearn.model_selection._validation._fit_and_score`."""
    if verbose > 1:
        if parameters is None:
            msg = ''
        else:
            msg = '%s' % (', '.join('%s=%s' % (k, v)
                                    for k, v in parameters.items()))
        print("[CV] %s %s" % (msg, (64 - len(msg)) * '.'))

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

    start_time = time.time()

    x_train, y_train = _safe_split(estimator, x_data, y_data, train)
    x_test, y_test = _safe_split(estimator, x_data, y_data, test, train)
    if sample_weights is not None:
        sample_weights_test = sample_weights[test]
    else:
        sample_weights_test = None

    try:
        if y_train is None:
            estimator.fit(x_train, **fit_params)
        else:
            estimator.fit(x_train, y_train, **fit_params)

    except Exception as exc:
        # Note fit time as time until error
        fit_time = time.time() - start_time
        score_time = 0.0
        if error_score == 'raise':
            raise
        if isinstance(error_score, numbers.Number):
            if isinstance(scorer, dict):
                test_scores = {name: error_score for name in scorer}
            else:
                test_scores = error_score
            warnings.warn("Estimator fit failed. The score on this train-test"
                          " partition for these parameters will be set to %f. "
                          "Details: \n%s" %
                          (error_score, format_exception_only(type(exc),
                                                              exc)[0]),
                          FitFailedWarning)
        else:
            raise ValueError("error_score must be the string 'raise' or a"
                             " numeric value. (Hint: if using 'raise', please"
                             " make sure that it has been spelled correctly.)")

    else:
        fit_time = time.time() - start_time
        test_scores = _score_weighted(estimator, x_test, y_test, scorer,
                                      sample_weights=sample_weights_test)
        score_time = time.time() - start_time - fit_time
    if verbose > 2:
        if isinstance(test_scores, dict):
            for scorer_name in sorted(test_scores):
                msg += ", %s=" % scorer_name
                msg += "%.3f" % test_scores[scorer_name]
        else:
            msg += ", score="
            msg += "%.3f" % test_scores

    if verbose > 1:
        total_time = score_time + fit_time
        print(_message_with_time('CV', msg, total_time))

    return [test_scores]


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


def _rfe_single_fit(rfe, estimator, x_data, y_data, train, test, scorer,
                    **fit_kwargs):
    """Return the score for a fit across one fold."""
    (x_train, y_train) = _safe_split(estimator, x_data, y_data, train)
    (x_test, y_test) = _safe_split(estimator, x_data, y_data, test, train)
    (fit_kwargs_train, _) = _split_fit_kwargs(fit_kwargs, train, test)

    def step_score(estimator, features):
        """Score for a single step in the recursive feature elimination."""
        return _score(estimator, x_test[:, features], y_test, scorer)

    return rfe._fit(x_train, y_train, step_score=step_score,
                    **fit_kwargs_train).scores_


def _score_weighted(estimator, x_test, y_test, scorer, sample_weights=None):
    """Expand :func:`sklearn.model_selection._validation._score`."""
    if isinstance(scorer, dict):
        # will cache method calls if needed. scorer() returns a dict
        scorer = _MultimetricScorer(**scorer)
    if y_test is None:
        scores = scorer(estimator, x_test, sample_weight=sample_weights)
    else:
        scores = scorer(estimator, x_test, y_test,
                        sample_weight=sample_weights)

    error_msg = ("scoring must return a number, got %s (%s) "
                 "instead. (scorer=%s)")
    if isinstance(scores, dict):
        for name, score in scores.items():
            if hasattr(score, 'item'):
                with suppress(ValueError):
                    # e.g. unwrap memmapped scalars
                    score = score.item()
            if not isinstance(score, numbers.Number):
                raise ValueError(error_msg % (score, type(score), name))
            scores[name] = score
    else:  # scalar
        if hasattr(scores, 'item'):
            with suppress(ValueError):
                # e.g. unwrap memmapped scalars
                scores = scores.item()
        if not isinstance(scores, numbers.Number):
            raise ValueError(error_msg % (scores, type(scores), scorer))
    return scores


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
    scorer_name = 'score'
    scoring = {scorer_name: scorer}
    x_data, y_data, groups = indexable(x_data, y_data, groups)

    cv = check_cv(cv, y_data, classifier=is_classifier(estimator))
    scorers, _ = _check_multimetric_scoring(estimator, scoring=scoring)

    # We clone the estimator to make sure that all the folds are
    # independent, and that it is pickle-able.
    parallel = Parallel(n_jobs=n_jobs, verbose=verbose,
                        pre_dispatch=pre_dispatch)
    scores = parallel(
        delayed(_fit_and_score_weighted)(
            clone(estimator), x_data, y_data, scorers, train, test, verbose,
            None, fit_params, error_score=error_score,
            sample_weights=sample_weights)
        for train, test in cv.split(x_data, y_data, groups))

    test_scores = list(zip(*scores))[0]
    test_scores = _aggregate_score_dicts(test_scores)

    return np.array(test_scores[scorer_name])


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
        tags = self._get_tags()
        x_data, y_data = check_X_y(
            x_data, y_data, "csc", ensure_min_features=2,
            force_all_finite=not tags.get('allow_nan', True))

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

        support_ = np.ones(n_features, dtype=np.bool)
        ranking_ = np.ones(n_features, dtype=np.int)

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
                raise RuntimeError('The classifier does not expose '
                                   '"coef_" or "feature_importances_" '
                                   'attributes')

            # Get ranks
            if coefs.ndim > 1:
                ranks = np.argsort(safe_sqr(coefs).sum(axis=0))
            else:
                ranks = np.argsort(safe_sqr(coefs))

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
        self.cv = cv
        self.scoring = scoring
        self.verbose = verbose
        self.n_jobs = n_jobs
        self.min_features_to_select = min_features_to_select

    def fit(self, x_data, y_data, groups=None, **fit_kwargs):
        """Expand :meth:`fit` to accept kwargs."""
        x_data, y_data = check_X_y(
            x_data, y_data, "csr", ensure_min_features=2,
            force_all_finite=False)

        # Initialization
        cv = check_cv(self.cv, y_data, is_classifier(self.estimator))
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
        # addition of n_jobs parameter in version 0.18.

        if effective_n_jobs(self.n_jobs) == 1:
            parallel, func = list, _rfe_single_fit
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

        # FIXME
        if fit_params['transformer']:
            raise NotImplementedError(
                f"Fit parameters {fit_params['transformer']} for transformer "
                f"{self.transformer.__class__} of {self.__class__} are not "
                f"supported at the moment")

        return (fit_params['transformer'], fit_params['regressor'])

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
        more_tags = super()._more_tags()
        more_tags['allow_nan'] = True
        return more_tags
