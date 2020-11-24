"""A collection of functions to calibrate the sigma values."""
import logging
import os
from typing import Union

import numpy as np
import xarray as xr
from scipy.optimize import minimize
from utilities import (
    area_weighted_mean,
    calculate_model_distances,
    combine_ensemble_members,
    compute_overall_mean,
    read_metadata,
    read_model_data,
    weighted_quantile,
)

logger = logging.getLogger(os.path.basename(__file__))


def calibrate_independence_sigma(overall_independence: 'xr.DataArray',
                                 performance_contributions: list,
                                 cfg: dict) -> float:
    """
    Parameters
    ----------
    overall_independence : array_like, shape (N, N)
    # TODO
    force : bool, optional
        If force is False (default) this function will raise an error if the
        independence_contributions are selected such that that no clear
        separateion between dependent and independent models is possible. Set
        to True (not recommendet) to scip this check.

    Returns
    -------
    independence_sigma : float
    """
    # do stuff
    return 0.5


def mask_dependent_models(overall_performance,
                          overall_independence,
                          subset_size=20):
    """
    Remove (set to nan?) the models closest related to the perfect model for each perfect model
    so that there are the subset_size most independent models remaining
    Alternatively: remove random models to test sensitivity
    """
    # TODO: do stuff
    return overall_performance


def calculate_weights_data(performance, independence, performance_sigma,
                           independence_sigma):

    # for each perfect model we need to remove dependent models from indepndence
    idx_nan = np.isnan(performance)
    numerator = np.exp(-((performance / performance_sigma)**2))

    if independence is not None:
        exp = np.exp(-(
            (independence[~idx_nan, ~idx_nan] / independence_sigma)**2))
        # Note diagonal = exp(0) = 1, thus this is equal to 1 + sum(i!=j)
        # TODO: at the moment I can not be sure that the reference dim is actually the second one
        # since the matrix is symmetric it does not matter but if that ever changes it will be wrong!
        denominator = exp.sum(axis=1)
    else:
        denominator = 1

    weights = numerator / denominator
    weights /= weights.sum(where=~idx_nan)

    return weights


def calculate_weights_perfect_model(performance_matrix, independence,
                                    performance_sigma, independence_sigma):
    # TODO: this could/should probably be done in 'mask_dependent_models'
    performance_matrix.values[np.diag_indices(
        performance_matrix['model_ensemble'].size)] = np.nan
    weights_matrix = xr.apply_ufunc(
        calculate_weights_data,
        performance_matrix,
        independence,
        performance_sigma,
        independence_sigma,
        input_core_dims=[['model_ensemble'], [], [], []],
        output_core_dims=[['model_ensemble']],
        vectorize=True,
    )

    return weights_matrix


def calculate_percentiles(data, weights, percentiles=(10, 90)):
    idx_keep = ~np.isnan(weights)
    data = data[idx_keep]
    weights = weights[idx_keep]
    quantiles = np.array(percentiles) / 100

    return (weighted_quantile(data, quantiles),
            weighted_quantile(data, quantiles, weights))


# Option 1 to calibrate sigma
def get_confidence_score(target, percentiles, weights_matrix,
                         performance_sigma):
    percentiles_data, percentiles_data_weighted = xr.apply_ufunc(
        calculate_percentiles,
        target,
        weights_matrix,
        input_core_dims=[['model_ensemble'], ['model_ensemble']],
        output_core_dims=[['percentile'], ['percentile']],
        vectorize=True,
        kwargs={'percentiles': percentiles})

    target = target.rename({'model_ensemble': 'model_ensemble_reference'})

    # inside_count_unweighted = xr.ufuncs.logical_and(
    #     target >= percentiles_data.isel(percentile=0),
    #     target <= percentiles_data.isel(percentile=1))

    inside_count_weighted = xr.ufuncs.logical_and(
        target >= percentiles_data_weighted.isel(percentile=0),
        target <= percentiles_data_weighted.isel(percentile=1)).values

    # inside_ratio_unweighted = inside_count_unweighted.sum() / len(inside_count_unweighted)
    # NOTE: this is only strictly speaking true if the perfect model (nor any other model)
    # is not removed from the ensemble. A better reference might be the 'inside_ration_unweighted'
    inside_ratio_reference = (percentiles[1] -
                              percentiles[0]) / 100 - 50  # DEBUG: remove - 50
    inside_ratio_weighted = inside_count_weighted.sum() / len(
        inside_count_weighted)

    # NOTE: in the simplest case we just want the smalest sigma which is not overconfident
    # but there might be reasons to use something slightly different

    confident = inside_ratio_weighted >= inside_ratio_reference
    # either the smallest confident sigma or the largerst sigma in total wins
    if confident:
        return performance_sigma
    return 9999 - performance_sigma


def evaluate_future(performance_sigma, overall_performance, target,
                    overall_independence, independence_sigma, settings):
    """For each perfect_model and each sigma compare the weighted mean to the
    'truth'.

    # NOTE: for each perfect model at least one model (the perfect model
    itself) but potentially more models (see remove_related_models)
    should have weight nan, this needs to be ignored (probably happens
    automatically but just to keep in mind)
    """
    performance_sigma = performance_sigma[0]

    if False:  # settings.get('use_independence_weighting', True):
        weights_matrix = calculate_weights_perfect_model(
            overall_performance, overall_independence, performance_sigma,
            independence_sigma)
    else:
        weights_matrix = calculate_weights_perfect_model(
            overall_performance, None, performance_sigma, None)
    """
    NOTE: I realise that calculating the weights and percentiles at once and also
    immediately checking the inside_ration would probably be more concise (and
    faster?). But I'm trying to separate the perfect model evaluation and the
    metric we are evaluating here to enable easy adding of more metrics and
    potential combining of multiple metrics.
    """

    if settings.get('optimise_for', 'confidence') == 'confidence':
        percentiles = settings.get('percentiles', (10, 90))
        return get_confidence_score(target, percentiles, weights_matrix,
                                    performance_sigma)
    elif settings['omptimise_for'] == 'TEST_METRIC_2':
        raise NotImplementedError
    # ...


def calibrate_performance_sigma(performance_contributions: list,
                                overall_independence: 'xr.DataArray',
                                independence_sigma: Union[float, None],
                                cfg: dict) -> float:

    settings = cfg['calibrate_performance_sigma']
    if (settings['performance_sigma_range'][0] <= 0 or
        settings['performance_sigma_range'][1] >= 10):
        errmsg = 'performance_sigma_range has to be >0 and < 10'
        raise IOError(errmsg)
    if (settings['performance_sigma_range'][0] >=
        settings['performance_sigma_range'][1]):
        errmsg = ' '.join([
            'performance_sigma_range needs to have the form (min, max) with',
            'min < max'])
        raise IOError(errmsg)

    models, _ = read_metadata(cfg)

    performances = {}
    for variable_group in performance_contributions:

        logger.info('Reading model data for %s', variable_group)
        datasets_model = models[variable_group]
        model_data, _ = read_model_data(datasets_model)

        logger.info('Calculating performance for %s', variable_group)
        performance = calculate_model_distances(model_data)
        logger.debug(performance.values)
        performances[variable_group] = performance

    performance = xr.Dataset(performances)
    overall_performance = compute_overall_mean(
        performance, performance_contributions)  # (N, N) matrix

    target = models[settings['variable_group_target']]
    target_data, target_data_files = read_model_data(target)
    target_data = area_weighted_mean(target_data)
    # TODO: log_provenance ?

    if cfg['combine_ensemble_members']:
        # TODO: this combines both dimensions -> I think about if we want this!
        # other option: only average the model dimension (equivalent to the actual performance weighting)
        # then do something else with the perfect model dimension (but this will probably make things more
        # complicated.
        # - We could run the 80/80 test with all members and then average the in/out count
        # for the models -> but what about the other metrics?
        # - We could only use the first member
        # - We could do some random sampling/ bootstrapping

        overall_performance, groups_performance = combine_ensemble_members(
            overall_performance)
        target_data, _ = combine_ensemble_members(target_data)

    # NOTE: what do we do here if 'combine_ensemble_members' is False?
    overall_performance = mask_dependent_models(
        overall_performance, overall_independence,
        settings['most_independent_subset'])

    # TODO: a often used output is the inside_ration versus sigma plot this can
    # not be provided from this (I think). Use a flag to chose between this and speed?
    res = minimize(
        evaluate_future,
        x0=overall_performance.mean(),
        bounds=[settings['performance_sigma_range']],
        args=(overall_performance, target_data, overall_independence,
              independence_sigma, settings),
    )

    if res.fun > 999:
        logmsg = ' '.join([
            'No confident sigma could be found! Using largest possible value.',
            'Bad choice of predictors or target or too small sigma range?'
        ])
        logger.warning(logmsg)

    return res.x


def calibrate_sigmas(independence_sigma: Union[float, str, None],
                     performance_sigma: Union[float, str, None],
                     overall_independence: 'xr.DataArray',
                     performance_contributions: list,
                     cfg: dict) -> (Union[float, None], Union[float, None]):
    """Calibibrate independence and/or performance sigma if not set."""

    if independence_sigma == 'calibrate':
        independence_sigma = calibrate_independence_sigma(
            overall_independence, cfg)
    if performance_sigma == 'calibrate':
        # TODO: in some cases we need the independence weighting here, if we do not
        # calculate it this needs to rais an error
        performance_sigma = calibrate_performance_sigma(
            performance_contributions, overall_independence,
            independence_sigma, cfg)

    # TODO: potential cross checks, sanity tests, perfect model skill calculation

    return independence_sigma, performance_sigma
