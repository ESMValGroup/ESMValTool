"""Implementation of the climwip weighting scheme.

Lukas Brunner et al. section 2.4
https://iopscience.iop.org/article/10.1088/1748-9326/ab492f
"""
import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import xarray as xr
from calibrate_sigmas import calibrate_performance_sigma
from core_functions import (
    area_weighted_mean,
    calculate_model_distances,
    calculate_weights,
    combine_ensemble_members,
    compute_overall_mean,
)
from io_functions import (
    log_provenance,
    read_input_data,
    read_metadata,
    read_model_data,
)

from esmvaltool.diag_scripts.shared import (
    get_diagnostic_filename,
    get_plot_filename,
    run_diagnostic,
)

logger = logging.getLogger(os.path.basename(__file__))


def read_observation_data(datasets: list) -> tuple:
    """Load observation data from list of metadata."""
    return read_input_data(datasets,
                           dim='obs_ensemble',
                           identifier_fmt='{dataset}')


def aggregate_obs_data(data_array: 'xr.DataArray',
                       operator: str = 'median') -> 'xr.DataArray':
    """Reduce data array along ensemble dimension.

    Apply the operator to squeeze the ensemble dimension by applying the
    `operator` along the ensemble (`obs_ensemble`) dimension. Returns an
    xarray.Dataset squeezed to 1D.
    """
    if operator == 'median':
        output = data_array.median(dim='obs_ensemble')
    else:
        raise ValueError(f'No such operator `{operator}`')

    return output


def visualize_and_save_independence(independence: 'xr.DataArray', cfg: dict,
                                    ancestors: list):
    """Visualize independence."""
    variable = independence.variable_group
    labels = list(independence.model_ensemble.values)

    figure, axes = plt.subplots(figsize=(15, 15),
                                subplot_kw={'aspect': 'equal'})
    chart = sns.heatmap(
        independence,
        linewidths=1,
        cmap="YlGn",
        xticklabels=labels,
        yticklabels=labels,
        cbar_kws={'label': f'Euclidean distance ({independence.units})'},
        ax=axes,
    )
    chart.set_title(f'Distance matrix for {variable}')

    filename_plot = get_plot_filename(f'independence_{variable}', cfg)
    figure.savefig(filename_plot, dpi=300, bbox_inches='tight')
    plt.close(figure)

    filename_data = get_diagnostic_filename(f'independence_{variable}',
                                            cfg,
                                            extension='nc')
    independence.to_netcdf(filename_data)

    caption = f'Euclidean distance matrix for variable {variable}'
    log_provenance(caption, filename_plot, cfg, ancestors)
    log_provenance(caption, filename_data, cfg, ancestors)


def calculate_performance(model_data: 'xr.DataArray',
                          obs_data: 'xr.DataArray') -> 'xr.DataArray':
    """Calculate performance.

    Calculate the area weighted mean between the given ensemble of model
    data, and observation data. The observation data must have the
    ensemble dimension squeezed or reduced. Returns an xarray.DataArray
    containing the same number of values as members of `model_data`.
    """
    diff = model_data - obs_data

    performance = area_weighted_mean(diff**2)**0.5

    performance.name = f'd{model_data.name}'
    performance.attrs['variable_group'] = model_data.name
    performance.attrs["units"] = model_data.units

    return performance


def barplot(metric: 'xr.DataArray', label: str, filename: str):
    """Visualize metric as barplot."""
    name = metric.name
    variable_group = metric.variable_group
    units = metric.units

    metric_df = metric.to_dataframe().reset_index()

    ylabel = f'{label} {variable_group} ({units})'

    figure, axes = plt.subplots(figsize=(15, 10))
    chart = sns.barplot(x='model_ensemble', y=name, data=metric_df, ax=axes)
    chart.set_xticklabels(chart.get_xticklabels(),
                          rotation=45,
                          horizontalalignment='right')
    chart.set_title(f'{label} for {variable_group}')
    chart.set_ylabel(ylabel)
    chart.set_xlabel('')

    figure.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close(figure)


def visualize_and_save_performance(performance: 'xr.DataArray', cfg: dict,
                                   ancestors: list):
    """Visualize performance."""
    label = 'RMS error'

    variable_group = performance.variable_group
    filename_plot = get_plot_filename(f'performance_{variable_group}', cfg)

    barplot(performance, label, filename_plot)

    filename_data = get_diagnostic_filename(f'performance_{variable_group}',
                                            cfg,
                                            extension='nc')
    performance.to_netcdf(filename_data)

    caption = f'Performance metric {label} for variable group {variable_group}'
    log_provenance(caption, filename_plot, cfg, ancestors)
    log_provenance(caption, filename_data, cfg, ancestors)


def split_ensemble_members(dataset: 'xr.DataArray',
                           groups: dict) -> 'xr.DataArray':
    """Split combined ensemble members of the same model."""
    model_ensemble = []
    nr_members = []
    for model in dataset['model_ensemble'].values:
        model_ensemble += groups[model]
        nr_members.append(len(groups[model]))

    data_scaled = dataset.values / nr_members
    data_expanded = np.repeat(data_scaled, nr_members)

    return xr.DataArray(data_expanded,
                        coords={'model_ensemble': model_ensemble},
                        dims='model_ensemble',
                        name=dataset.name,
                        attrs=dataset.attrs)


def visualize_and_save_weights(weights: 'xr.DataArray', cfg: dict,
                               ancestors: list):
    """Visualize weights."""
    label = 'Weights'

    filename_plot = get_plot_filename('weights', cfg)

    barplot(weights, label, filename_plot)

    filename_data = get_diagnostic_filename('weights', cfg, extension='nc')
    weights.to_netcdf(filename_data)

    caption = 'Weights'
    log_provenance(caption, filename_plot, cfg, ancestors)
    log_provenance(caption, filename_data, cfg, ancestors)


def parse_contributions_sigma(metric: str, cfg: dict) -> dict:
    """Return contributions > 0 and sigma for a given metric."""
    if cfg.get(f'{metric}_contributions') is None:  # not set or set to None
        contributions = {}
    else:
        contributions = {
            key: value
            for key, value in cfg[f'{metric}_contributions'].items()
            if value > 0
        }
    sigma = cfg.get(f'{metric}_sigma')
    return contributions, sigma


def main(cfg):
    """Perform climwip weighting method."""
    models, observations = read_metadata(cfg)

    independence_contributions, independence_sigma = parse_contributions_sigma(
        'independence', cfg)
    performance_contributions, performance_sigma = parse_contributions_sigma(
        'performance', cfg)

    if not (independence_contributions or performance_contributions):
        errmsg = ' '.join([
            'Either the independence_contributions or the',
            'performance_contributions field need to be set and contain at',
            'least one variable group with weight > 0 otherwise no weights',
            'can be calculated!'
        ])
        raise IOError(errmsg)

    model_ancestors = []
    obs_ancestors = []

    performances = {}
    independences = {}

    for variable_group in independence_contributions:

        logger.info('Reading model data for %s', variable_group)
        datasets_model = models[variable_group]
        model_data, model_data_files = read_model_data(datasets_model)

        logger.info('Calculating independence for %s', variable_group)
        independence = calculate_model_distances(model_data)
        visualize_and_save_independence(independence, cfg, model_data_files)
        logger.debug(independence.values)
        independences[variable_group] = independence

        model_ancestors.extend(model_data_files)

    for variable_group in performance_contributions:

        logger.info('Reading model data for %s', variable_group)
        datasets_model = models[variable_group]
        model_data, model_data_files = read_model_data(datasets_model)

        logger.info('Reading observation data for %s', variable_group)
        datasets_obs = observations[variable_group]
        obs_data, obs_data_files = read_observation_data(datasets_obs)
        obs_data = aggregate_obs_data(obs_data, operator='median')

        logger.info('Calculating performance for %s', variable_group)
        performance = calculate_performance(model_data, obs_data)
        visualize_and_save_performance(performance, cfg,
                                       model_data_files + obs_data_files)
        logger.debug(performance.values)
        performances[variable_group] = performance
        obs_ancestors.extend(obs_data_files)

        model_ancestors.extend(model_data_files)

    model_ancestors = list(set(model_ancestors))  # only keep unique items

    if independence_contributions:
        logger.info('Computing overall mean independence')
        independence = xr.Dataset(independences)
        overall_independence = compute_overall_mean(
            independence, independence_contributions)
        visualize_and_save_independence(overall_independence, cfg,
                                        model_ancestors)
        if independence_sigma is None:
            raise NotImplementedError('`independence_sigma` must be set if '
                                      '`independence_contributions` is set')
    else:
        overall_independence = None

    if performance_contributions:
        logger.info('Computing overall mean performance')
        performance = xr.Dataset(performances)
        overall_performance = compute_overall_mean(performance,
                                                   performance_contributions)
        visualize_and_save_performance(overall_performance, cfg,
                                       model_ancestors + obs_ancestors)
        if performance_sigma is None:
            performance_sigma = calibrate_performance_sigma(
                performance_contributions, overall_independence,
                independence_sigma, cfg)
    else:
        overall_performance = None

    if cfg['combine_ensemble_members']:
        overall_independence, groups_independence = combine_ensemble_members(
            overall_independence,
            ['model_ensemble', 'model_ensemble_reference'])
        overall_performance, groups_performance = combine_ensemble_members(
            overall_performance)
        # one of them could be empty if metric is not calculated
        groups = {**groups_independence, **groups_performance}

    logger.info('Calculating weights')
    weights = calculate_weights(overall_performance, overall_independence,
                                performance_sigma, independence_sigma)

    if cfg['combine_ensemble_members']:
        weights = split_ensemble_members(weights, groups)

    visualize_and_save_weights(weights,
                               cfg,
                               ancestors=model_ancestors + obs_ancestors)
    logger.debug(weights.values)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
