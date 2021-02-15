"""A collection of functions to calibrate the shape parameters (sigmas)."""
import logging
import os
from typing import Union

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from core_functions import (
    area_weighted_mean,
    calculate_model_distances,
    calculate_weights,
    combine_ensemble_members,
    compute_overall_mean,
    weighted_quantile,
)
from io_functions import read_metadata, read_model_data
from scipy.optimize import brute

from esmvaltool.diag_scripts.shared import (
    get_diagnostic_filename,
    get_plot_filename,
)

logger = logging.getLogger(os.path.basename(__file__))

SIGMA_RANGE = (.1, 2)  # allow this to be set by the recipe later
PERCENTILES = [.1, .9]  # allow this to be set by the recipe later

confidence_test_values = {'baseline': {}}


def calculate_percentiles(target: 'xr.DataArray',
                          weights_matrix: 'xr.DataArray',
                          percentiles: list,
                          weighted: bool = True) -> ('xr.DataArray', float):
    """Calculate (equally) weighted percentiles based on each perfect model.

    Parameters
    ----------
    target : array_like, shape (N,)
        Array of model values that will be evaluated in order to estimate the
        optimal sigma value. For each perfect model, the perfect model (and
        potentially additional models) are excluded from the target, the
        rest is weighted. The perfect model is then used to evaluate the
        weighted distribution (see also weights_matrix).
    weights_matrix : array_like, shape (N, N)
        For each perfect model in the perfect_model_ensemble dimension
        the weights_matrix contains the respective model weights in the
        model_ensemble dimension based on this perfect model.

        Special feature: nan values in the model_ensemble dimension will lead
        to the model being excluded from the weights calculation. This
        is always the case for the perfect model itself (diagonal of the
        matrix) but might also be the case for other models. This is
        particularly important for a correct calculation of the independence
        weighting.
    percentiles : array_like, shape (2,)
        Lower and upper percentile to use in the confidence test. Has to
        satisfy 0 <= percentiles <=1 and percentiles[0] < percentiles[1]
    weighted : bool, optional
        If weighted is set to False (default: True) all values in the weights
        matrix will be set to 1, except for nan values which will be preserved.
        This can be used to calculate the unweighted baseline for each perfect
        model, consistent with the weighted case.

    Returns
    -------
    percentile_spread : array_like, shape (N,)
        Full range spanned by the two percentiles for each perfect model.
    inside_ratio: float
        Ratio of perfect models inside their respective percentile_spread.
    """
    percentiles = xr.DataArray(percentiles, dims='percentile')

    # need to rename for the inside_count comparison
    target_perfect = target.rename(
        {'model_ensemble': 'perfect_model_ensemble'})

    if not weighted:  # replace with equal weight but keep nans
        weights_matrix = 0 * weights_matrix + 1

    percentiles_data = xr.apply_ufunc(
        weighted_quantile,
        target,
        percentiles,
        weights_matrix,
        input_core_dims=[['model_ensemble'], ['percentile'],
                         ['model_ensemble']],
        output_core_dims=[['percentile']],
        vectorize=True,
    )

    inside_count = xr.ufuncs.logical_and(
        target_perfect >= percentiles_data.isel(percentile=0),
        target_perfect <= percentiles_data.isel(percentile=1)).values
    inside_ratio = inside_count.sum() / len(inside_count)

    percentiles_spread = (percentiles_data.isel(percentile=1) -
                          percentiles_data.isel(percentile=0))

    return percentiles_spread, inside_ratio


def compute_cost_function(target: 'xr.DataArray',
                          weights_matrix: 'xr.DataArray',
                          performance_sigma: float) -> float:
    """Optimize the performance sigma for confidence.

    Parameters
    ----------
    target : array_like, shape (N,)
        See calculate_percentiles for more information.
    weights_matrix : array_like, shape (N, N)
        See calculate_percentiles for more information.
    performance_sigma : float
        The performance sigma value used to calculate the weights.

    Returns
    -------
    cost_function_value : float
        Value of the cost function based on the given performance sigma value.
        The cost function is a discontinuous function distinguishing two cases:
              99 + abs(difference)  if overconfident
        f = {
              sigma                 else

    WARNING
    -------
    It is highly recommended to visually inspect the graphical output of this
    process to ensure the optimisation worked as intended.


    Additional information
    ----------------------
    After evaluating all possible sigma values the sigma which leads to the
    smallest cost function will be selected. Two different cases need to be
    disdinguished:
    * the cost function is > 99 for all cases:
        All sigma values lead to overconfident weighting. The script will
        plot a visualization of the test and will then raise an error. The user
        can decide to select a sigma value manually (e.g., if the test only
        failed by a narrow margin) and set it in the recipe.
    * the cost function is < 99 for at least one sigma:
        For this case the cost function increases linearly with sigma,
        therefore the smallest possible sigma where this is true will be
        selected on the fly and used in the computation of weights.
    """
    percentiles = PERCENTILES
    inside_ratio_reference = percentiles[1] - percentiles[0]

    # calculate the equally weighted case once as baseline
    if len(confidence_test_values['baseline']) == 0:
        percentiles_spread, inside_ratio = calculate_percentiles(
            target, weights_matrix, percentiles, weighted=False)
        confidence_test_values['baseline']['percentile_spread'] = (
            percentiles_spread)
        confidence_test_values['baseline']['inside_ratio'] = inside_ratio

    percentiles_spread, inside_ratio = calculate_percentiles(
        target, weights_matrix, percentiles)
    confidence_test_values[performance_sigma] = {
        'percentile_spread': percentiles_spread,
        'inside_ratio': inside_ratio
    }

    difference = inside_ratio - inside_ratio_reference

    if difference < 0:  # overconfident
        return 99 - difference
    return performance_sigma


def evaluate_target(performance_sigma: list,
                    overall_performance: 'xr.DataArray',
                    target: 'xr.DataArray',
                    overall_independence: 'xr.DataArray',
                    independence_sigma: float) -> float:
    """Evaluate the weighting in the target period.

    Parameters
    ----------
    performance_sigma : list of one float
        Performance weighting shape parameter, determines how strong the
        weighting for performance is (smaller values correspond to stronger
        weighting)
    overall_performance : array_like, shape (N, N)
        Contains the generalised distance for each model in the model_ensemble
        dimension for each perfect model in the perfect_model_ensemble
        dimension.
    target : array_like, shape (N,)
        See calculate_percentiles for more information.
    overall_independence : array_like, shape (N, N)
        Matrix containing model-model distances for independence.
    independence_sigma : float
        Independence weighting shape parameter.

    Returns
    -------
    cost_function_value : float
        See compute_cost_function for more information.
    """
    performance_sigma = performance_sigma[0]

    # exclude perfect model in each row by setting it to nan
    idx_diag = np.diag_indices(overall_performance['model_ensemble'].size)
    overall_performance.values[idx_diag] = np.nan

    weights_matrix = calculate_weights(overall_performance,
                                       overall_independence, performance_sigma,
                                       independence_sigma)

    cost_function_value = compute_cost_function(target, weights_matrix,
                                                performance_sigma)
    return cost_function_value


def visualize_save_calibration(performance_sigma, cfg, success):
    """Visualize a summary of the calibration."""
    percentiles = PERCENTILES
    inside_ratio_reference = percentiles[1] - percentiles[0]

    baseline = confidence_test_values.pop('baseline')
    sigmas = sorted(confidence_test_values)
    inside_ratios = [
        confidence_test_values[sigma]['inside_ratio'] for sigma in sigmas
    ]
    confidence = xr.Dataset(
        data_vars={
            'inside_ratio_reference': ((), inside_ratio_reference, {
                'units': '1'
            }),
            'inside_ratio': ('sigma', inside_ratios, {
                'units': '1'
            }),
            'sigma': ('sigma', sigmas, {
                'units': '1'
            }),
        })

    figure, axes = plt.subplots(figsize=(12, 8))
    axes.plot(sigmas,
              inside_ratios,
              color='k',
              lw=2,
              label='Perfect models inside the {:.0%} range'.format(
                  inside_ratio_reference),
              zorder=99)
    axes.axhline(inside_ratio_reference,
                 color='k',
                 ls='--',
                 label='Reference: {:.0%}'.format(inside_ratio_reference))
    axes.axhline(baseline['inside_ratio'],
                 color='k',
                 ls=':',
                 label='Unweighted baseline: {:.0%}'.format(
                     baseline['inside_ratio']))

    if success:
        axes.axvline(performance_sigma,
                     color='k',
                     ls='-.',
                     label='Selected performance sigma: {:.2f}'.format(
                         performance_sigma))
        axes.set_title('Performance sigma calibration')
    else:
        axes.axvline(
            performance_sigma,
            color='gray',
            ls='-.',
            label='Best performance sigma: {:.2f} (set manually to use)'.
            format(performance_sigma))
        axes.set_title('Performance sigma calibration (FAILED)')

    # optional: sharpness
    sharpness = xr.concat([
        confidence_test_values[sigma]['percentile_spread'].expand_dims(
            {'sigma': [sigma]}) for sigma in sigmas], dim='sigma')
    sharpness /= baseline['percentile_spread']
    axes.plot(sigmas,
              sharpness.mean('perfect_model_ensemble'),
              color='lightgray',
              label='Spread relative to unweighted (mean & 80% range)')
    axes.fill_between(sigmas,
                      *sharpness.quantile((.1, .9), 'perfect_model_ensemble'),
                      facecolor='lightgray',
                      edgecolor='none',
                      alpha=.3)
    axes.axhline(1, color='lightgray', ls='--')

    sharpness.attrs['units'] = '1'
    confidence['sharpness'] = sharpness

    axes.set_xlim(SIGMA_RANGE)
    axes.set_ylim(0, 1.3)
    axes.set_xlabel('sigma (1)')
    axes.set_ylabel('Ratio (1)')

    axes.legend()

    filename_plot = get_plot_filename('performance_sigma_calibration', cfg)
    figure.savefig(filename_plot, dpi=300, bbox_inches='tight')
    plt.close(figure)

    filename_data = get_diagnostic_filename('performance_sigma_calibration',
                                            cfg,
                                            extension='nc')
    confidence.to_netcdf(filename_data)


def calibrate_performance_sigma(
        performance_contributions: list,
        overall_independence: Union['xr.DataArray', None],
        independence_sigma: Union[float, None],
        cfg: dict) -> float:
    """Calibrate the performance sigma using a perfect model approach."""
    settings = cfg['calibrate_performance_sigma']
    models, _ = read_metadata(cfg)

    performances_matrix = {}
    for variable_group in performance_contributions:

        logger.info('Reading model data for %s', variable_group)
        datasets_model = models[variable_group]
        model_data, _ = read_model_data(datasets_model)

        logger.info('Calculating performance for %s', variable_group)
        performance_matrix = calculate_model_distances(
            model_data, 'perfect_model_ensemble')
        logger.debug(performance_matrix.values)
        performances_matrix[variable_group] = performance_matrix

    performance_matrix = xr.Dataset(performances_matrix)
    overall_performance = compute_overall_mean(performance_matrix,
                                               performance_contributions)

    target = models[settings['target']]
    target_data, _ = read_model_data(target)
    target_data = area_weighted_mean(target_data)

    if cfg['combine_ensemble_members']:
        overall_independence, _ = combine_ensemble_members(
            overall_independence,
            ['model_ensemble', 'model_ensemble_reference'])
        overall_performance, _ = combine_ensemble_members(
            overall_performance, ['model_ensemble', 'perfect_model_ensemble'])
        target_data, _ = combine_ensemble_members(target_data)

    performance_sigma, fval, _, _ = brute(
        evaluate_target,
        ranges=(SIGMA_RANGE, ),
        Ns=100,
        finish=None,
        args=(overall_performance, target_data, overall_independence,
              independence_sigma),
        full_output=True,
    )

    success = fval < 99
    visualize_save_calibration(performance_sigma, cfg, success=success)

    if success:
        logmsg = f'Found optimal performance sigma value: {performance_sigma}'
        logger.info(logmsg)
        return performance_sigma

    errmsg = ('No confident sigma could be found! Bad choice of predictors or '
              'target or too small sigma range?')
    raise ValueError(errmsg)
