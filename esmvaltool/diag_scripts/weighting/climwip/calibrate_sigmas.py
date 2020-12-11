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

SIGMA_RANGE = (.1, 2)
confidence_test_values = {'baseline': {}}


def optimize_confidence(target, weights_matrix, performance_sigma):
    """Optimize the performance_sigma for confidence."""
    percentiles = xr.DataArray([.1, .9], dims='percentile')
    inside_ratio_reference = (percentiles[1] - percentiles[0]).values

    # need to rename for the inside_count comparison
    target_perfect = target.rename(
        {'model_ensemble': 'perfect_model_ensemble'})

    # calculate the equally weighted case once as baseline
    if len(confidence_test_values['baseline']) == 0:
        equal_weights_matrix = weights_matrix * 0 + 1
        percentiles_data = xr.apply_ufunc(
            weighted_quantile,
            target,
            percentiles,
            equal_weights_matrix,
            input_core_dims=[['model_ensemble'], ['percentile'],
                             ['model_ensemble']],
            output_core_dims=[['percentile']],
            vectorize=True,
        )
        inside_count = xr.ufuncs.logical_and(
            target_perfect >= percentiles_data.isel(percentile=0),
            target_perfect <= percentiles_data.isel(percentile=1)).values
        inside_ratio_baseline = inside_count.sum() / len(inside_count)

        confidence_test_values['baseline'] = {
            'inside_ratio':
            inside_ratio_baseline,
            'inside_ratio_reference':
            inside_ratio_reference,
            'percentiles_spread': (percentiles_data.isel(percentile=1) -
                                   percentiles_data.isel(percentile=0)),
        }

    percentiles_data_weighted = xr.apply_ufunc(
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
        target_perfect >= percentiles_data_weighted.isel(percentile=0),
        target_perfect <= percentiles_data_weighted.isel(percentile=1)).values
    inside_ratio = inside_count.sum() / len(inside_count)
    ratio_difference = inside_ratio - inside_ratio_reference

    # save intermediate results
    percentiles_spread_weighted = (
        percentiles_data_weighted.isel(percentile=1) -
        percentiles_data_weighted.isel(percentile=0))
    confidence_test_values[performance_sigma] = {
        'inside_ratio':
        inside_ratio,
        'sharpness':
        (percentiles_spread_weighted /
         confidence_test_values['baseline']['percentiles_spread']),
    }

    if ratio_difference < 0:
        return np.abs(ratio_difference) + SIGMA_RANGE[1] - performance_sigma
    return performance_sigma - SIGMA_RANGE[1]


def evaluate_target(performance_sigma, overall_performance, target,
                    overall_independence, independence_sigma):
    """Evaluate the weighting in the target period."""
    performance_sigma = performance_sigma[0]

    # exclude perfect model in each row by setting it to nan
    idx_diag = np.diag_indices(overall_performance['model_ensemble'].size)
    overall_performance.values[idx_diag] = np.nan

    weights_matrix = calculate_weights(overall_performance,
                                       overall_independence, performance_sigma,
                                       independence_sigma)

    cost_function = optimize_confidence(target, weights_matrix,
                                        performance_sigma)
    return cost_function


def visualize_save_calibration(performance_sigma, costf, cfg):
    """Visualize a summary of the calibration."""
    baseline = confidence_test_values.pop('baseline')
    sigmas = sorted(confidence_test_values)
    inside_ratios = [
        confidence_test_values[sigma]['inside_ratio'] for sigma in sigmas
    ]
    confidence = xr.Dataset(
        data_vars={
            'inside_ratio_reference': (
                (), baseline['inside_ratio_reference'], {'units': '1'}),
            'inside_ratio': ('sigma', inside_ratios, {'units': '1'}),
            'sigma': ('sigma', sigmas, {'units': '1'}),
        })

    figure, axes = plt.subplots(figsize=(12, 8))
    axes.plot(sigmas,
              inside_ratios,
              color='k',
              lw=2,
              label='Perfect models inside the {:.0%} range'.format(
                  baseline['inside_ratio_reference']),
              zorder=99)
    axes.axhline(baseline['inside_ratio_reference'],
                 color='k',
                 ls='--',
                 label='Reference: {:.0%}'.format(
                     baseline['inside_ratio_reference']))
    axes.axhline(baseline['inside_ratio'],
                 color='k',
                 ls=':',
                 label='Unweighted baseline: {:.0%}'.format(
                     baseline['inside_ratio']))
    axes.axvline(performance_sigma, color='k', ls='-.',
                 label='Selected performance sigma: {:.2f}'.format(
                     performance_sigma))

    # optional: sharpness
    sharpness = xr.concat(
        [confidence_test_values[sigma]['sharpness'].expand_dims(
            {'sigma': [sigma]}) for sigma in sigmas],
        dim='sigma')
    axes.plot(sigmas,
              sharpness.mean('perfect_model_ensemble'),
              color='lightgray',
              label='Spread relative to unweighted (mean & 80% range)')
    axes.fill_between(
        sigmas,
        *sharpness.quantile((.1, .9), 'perfect_model_ensemble'),
        facecolor='lightgray',
        edgecolor='none',
        alpha=.3)
    axes.axhline(1, color='lightgray', ls='--')

    sharpness.attrs['units'] = '1'
    confidence['sharpness'] = sharpness
    # ---

    # optional: cost function on second yaxis
    axes2 = axes.twinx()
    axes2.plot(sigmas, costf, color='red', lw=.5)
    axes2.set_ylabel('Cost function (1)')
    axes2.yaxis.label.set_color('red')
    axes2.tick_params(axis='y', colors='red')

    costf = xr.DataArray(costf,
                         dims=['sigma'],
                         coords={'sigma': sigmas},
                         attrs={'units': '1'})
    confidence['cost_function'] = costf

    axes.set_xlim(SIGMA_RANGE)
    axes.set_ylim(0, 1.3)
    axes.set_xlabel('sigma (1)')
    axes.set_ylabel('Ratio (1)')
    axes.set_title('Performance sigma calibration')

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

    performance_sigma, _, _, fgrid = brute(
        evaluate_target,
        ranges=(SIGMA_RANGE, ),
        Ns=100,
        finish=None,
        args=(overall_performance, target_data, overall_independence,
              independence_sigma),
        full_output=True,
    )

    visualize_save_calibration(performance_sigma, fgrid, cfg)

    if performance_sigma == SIGMA_RANGE[1]:
        logmsg = ' '.join([
            'No confident sigma could be found! Using largest possible value.',
            'Bad choice of predictors or target or too small sigma range?'
        ])
        logger.warning(logmsg)
    else:
        logmsg = f'Found optimal performance sigma value: {performance_sigma}'
        logger.info(logmsg)

    return performance_sigma
