"""Implementation of the climwip weighting scheme.

Lukas Brunner et al. section 2.4
https://iopscience.iop.org/article/10.1088/1748-9326/ab492f
"""
import logging
import os
from collections import defaultdict
from datetime import datetime

import matplotlib.pyplot as plt
import natsort
import numpy as np
import seaborn as sns
import xarray as xr
from scipy.spatial.distance import pdist, squareform

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(os.path.basename(__file__))


def get_provenance_record(caption: str, ancestors: list):
    """Create a provenance record describing the diagnostic data and plots."""
    record = {
        'caption':
        caption,
        'domains': ['reg'],
        'authors': [
            'kalverla_peter',
            'smeets_stef',
            'brunner_lukas',
            'camphuijsen_jaro',
        ],
        'references': [
            'brunner2019',
            'lorenz2018',
            'knutti2017',
        ],
        'ancestors':
        ancestors,
    }
    return record


def log_provenance(caption: str, filename: str, cfg: dict, ancestors: list):
    """Log provenance info."""
    provenance_record = get_provenance_record(caption, ancestors=ancestors)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance_record)

    logger.info('Output stored as %s', filename)


def read_metadata(cfg: dict) -> tuple:
    """Read the metadata from the config file.

    Returns a two dicts, one for the model data and one for the
    observational data. They are split based on the value of the
    'obs_data' variable in the recipe. The dictionaries are sorted by
    the variable.
    """
    obs_ids = cfg['obs_data']
    if isinstance(obs_ids, str):
        obs_ids = [obs_ids]

    input_data = cfg['input_data'].values()
    projects = group_metadata(input_data, attribute='project')

    observations = defaultdict(list)
    models = defaultdict(list)

    for project, metadata in projects.items():

        for item in metadata:
            variable_group = item['variable_group']

            if project in obs_ids:
                observations[variable_group].append(item)
            else:
                models[variable_group].append(item)

    return models, observations


def make_standard_calendar(xrds: 'xr.Dataset'):
    """Make sure time coordinate uses the default calendar.

    Workaround for incompatible calendars 'standard' and 'no-leap'.
    Assumes yearly data.
    """
    try:
        years = xrds.time.dt.year.values
        xrds['time'] = [datetime(year, 7, 1) for year in years]
    except TypeError:
        # Time dimension is 0-d array
        pass
    except AttributeError:
        # Time dimension does not exist
        pass


def read_input_data(metadata: list,
                    dim: str = 'data_ensemble',
                    identifier_fmt: str = '{dataset}') -> tuple:
    """Load data from metadata.

    Read the input data from the list of given data sets. `metadata` is
    a list of metadata containing the filenames to load. Only returns
    the given `variable`. The datasets are stacked along the `dim`
    dimension. Returns an xarray.DataArray.
    """
    data_arrays = []
    identifiers = []
    input_files = []
    for info in metadata:
        filename = info['filename']
        short_name = info['short_name']
        variable_group = info['variable_group']

        xrds = xr.open_dataset(filename)
        make_standard_calendar(xrds)
        xrda = xrds[short_name]
        xrda = xrda.rename(variable_group)
        data_arrays.append(xrda)

        identifier = identifier_fmt.format(**info)
        identifiers.append(identifier)
        input_files.append(filename)

    diagnostic = xr.concat(data_arrays, dim=dim)
    diagnostic[dim] = identifiers

    # Clean up unnecessary coordinate info
    redundant_dims = np.setdiff1d(diagnostic.coords, diagnostic.dims)
    diagnostic = diagnostic.drop(redundant_dims)

    # Use natural sorting order
    sorting = natsort.natsorted(identifiers, alg=natsort.IC)
    diagnostic = diagnostic.sel(indexers={dim: sorting})

    return diagnostic, input_files


def read_model_data(datasets: list) -> tuple:
    """Load model data from list of metadata."""
    return read_input_data(datasets,
                           dim='model_ensemble',
                           identifier_fmt='{dataset}_{ensemble}_{exp}')


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


def area_weighted_mean(data_array: 'xr.DataArray') -> 'xr.DataArray':
    """Calculate area mean weighted by the latitude.

    Returns a data array consisting of N values, where N == number of
    ensemble members.
    """
    weights_lat = np.cos(np.radians(data_array.lat))
    means = data_array.weighted(weights_lat).mean(dim=['lat', 'lon'])

    return means


def distance_matrix(values: 'np.ndarray',
                    weights: 'np.ndarray' = None) -> 'np.ndarray':
    """Calculate the pairwise distance between model members.

    Takes a dataset with ensemble member/lon/lat. Flattens lon/lat
    into a single dimension. Calculates the distance between every
    ensemble member.

    If weights are passed, they should have the same shape as values.

    Returns 2D NxN array, where N == number of ensemble members.
    """
    n_members = values.shape[0]

    values = values.reshape(n_members, -1)

    # pdist does not work with NaN
    not_nan = np.where(np.all(np.isfinite(values), axis=0))[0]
    values = values[:, not_nan]

    if weights is not None:
        # Reshape weights to match values array
        weights = weights.reshape(n_members, -1)
        weights = weights[:, not_nan]
        weights = weights[0]  # Weights are equal along first dim

    d_matrix = squareform(pdist(values, metric='euclidean', w=weights))

    return d_matrix


def calculate_independence(data_array: 'xr.DataArray') -> 'xr.DataArray':
    """Calculate independence.

    The independence is calculated as a distance matrix between the
    datasets defined in the `data_array`. Returned is a square matrix
    with where the number of elements along each edge equals the number
    of ensemble members.
    """
    weights = np.cos(np.radians(data_array.lat))
    weights, _ = xr.broadcast(weights, data_array)

    diff = xr.apply_ufunc(
        distance_matrix,
        data_array,
        weights,
        input_core_dims=[['model_ensemble', 'lat', 'lon'],
                         ['model_ensemble', 'lat', 'lon']],
        output_core_dims=[['perfect_model_ensemble', 'model_ensemble']],
    )

    diff.name = f'd{data_array.name}'
    diff.attrs['variable_group'] = data_array.name
    diff.attrs["units"] = data_array.units
    diff['perfect_model_ensemble'] = diff.model_ensemble.values

    return diff


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


def compute_overall_mean(dataset):
    """Normalize all variables in a dataset and return their mean."""
    if 'perfect_model_ensemble' in dataset.dims:
        median_dim = ['perfect_model_ensemble', 'model_ensemble']
    else:
        median_dim = 'model_ensemble'
    normalized = dataset / dataset.median(dim=median_dim)

    overall_mean = normalized.to_array(
        dim='variable_group').mean('variable_group')
    overall_mean.name = 'overall_mean'
    overall_mean.attrs['variable_group'] = 'overall_mean'
    overall_mean.attrs['units'] = '1'
    return overall_mean


def calculate_weights(performance: 'xr.DataArray',
                      independence: 'xr.DataArray', sigma_performance: float,
                      sigma_independence: float) -> 'xr.DataArray':
    """Calculate the (NOT normalised) weights for each model N.

    Parameters
    ----------
    performance : array_like, shape (N,)
        Array specifying the model performance.
    independence : array_like, shape (N, N)
        Array specifying the model independence.
    sigma_performance : float
        Sigma value defining the form of the weighting function
            for the performance.
    sigma_independence : float
        Sigma value defining the form of the weighting function
            for the independence.

    Returns
    -------
    weights : ndarray, shape (N,)
    """
    numerator = np.exp(-((performance / sigma_performance)**2))
    exp = np.exp(-((independence / sigma_independence)**2))

    # Note diagonal = exp(0) = 1, thus this is equal to 1 + sum(i!=j)
    denominator = exp.sum(axis=0)

    weights = numerator / denominator

    # Normalize weights
    weights /= weights.sum()

    weights.name = 'weight'
    weights.attrs['variable_group'] = 'weight'  # used in barplot
    weights.attrs['units'] = '1'

    return weights


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


def main(cfg):
    """Perform climwip weighting method."""
    models, observations = read_metadata(cfg)

    variable_groups = models.keys()

    model_ancestors = []
    obs_ancestors = []

    performances = {}
    independences = {}

    for variable_group in variable_groups:

        logger.info('Reading model data for %s', variable_group)
        datasets_model = models[variable_group]
        model_data, model_data_files = read_model_data(datasets_model)

        logger.info('Reading observation data for %s', variable_group)
        datasets_obs = observations[variable_group]
        obs_data, obs_data_files = read_observation_data(datasets_obs)
        obs_data = aggregate_obs_data(obs_data, operator='median')

        if datasets_model[0].get('independence', True):
            logger.info('Calculating independence for %s', variable_group)
            independence = calculate_independence(model_data)
            visualize_and_save_independence(independence, cfg,
                                            model_data_files)
            logger.debug(independence.values)
            independences[variable_group] = independence

        if datasets_model[0].get('performance', True):
            logger.info('Calculating performance for %s', variable_group)
            performance = calculate_performance(model_data, obs_data)
            visualize_and_save_performance(performance, cfg,
                                           model_data_files + obs_data_files)
            logger.debug(performance.values)
            performances[variable_group] = performance
            obs_ancestors.extend(obs_data_files)

        model_ancestors.extend(model_data_files)

    logger.info('Computing overall mean independence')
    independence = xr.Dataset(independences)
    overall_independence = compute_overall_mean(independence)
    visualize_and_save_independence(overall_independence, cfg, model_ancestors)

    logger.info('Computing overall mean performance')
    performance = xr.Dataset(performances)
    overall_performance = compute_overall_mean(performance)
    visualize_and_save_performance(overall_performance, cfg,
                                   model_ancestors + obs_ancestors)

    logger.info('Calculating weights')
    sigma_performance = cfg['sigma_performance']
    sigma_independence = cfg['sigma_independence']
    weights = calculate_weights(overall_performance, overall_independence,
                                sigma_performance, sigma_independence)
    visualize_and_save_weights(weights,
                               cfg,
                               ancestors=model_ancestors + obs_ancestors)
    logger.debug(weights.values)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
