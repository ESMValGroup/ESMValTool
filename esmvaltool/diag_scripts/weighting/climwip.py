"""
Implementation of the climwip weighting scheme

Lukas Brunner et al. section 2.4
https://iopscience.iop.org/article/10.1088/1748-9326/ab492f
"""
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import xarray as xr
from scipy.spatial.distance import pdist, squareform
import yaml
from datetime import datetime

from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            get_plot_filename, run_diagnostic,
                                            select_metadata)


def get_provenance_record(caption, ancestors):
    """Create a provenance record describing the diagnostic data and plots."""

    record = {
        'caption': caption,
        'domains': ['reg'],
        'authors': [
            'kalverla_peter',
            # 'smeets_stef',
            # 'brunner_lukas',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': ancestors,
    }
    return record


def read_metadata(cfg, projects: list) -> dict:
    """Read the metadata from the config file.

    Project is the list of project identifiers (str) to select."""

    d = defaultdict(list)

    for project in projects:
        metadata = select_metadata(cfg['input_data'].values(), project=project)

        for item in metadata:
            variable = item['short_name']

            d[variable].append(item)

    return d


def make_standard_calendar(xrds):
    """Make sure time coordinate uses the default calendar.

    Workaround for imcompatible calendars 'standard' and 'no-leap'.
    """
    try:
        years = xrds.time.dt.year.values
        months = xrds.time.dt.month.values
        days = xrds.time.dt.day.values
        xrds['time'] = [datetime(year, month, day) for year, month, day in zip(years, months, days)]
    except TypeError:
        # Time dimension is 0-d array
        pass
    except AttributeError:
        # Time dimension does not exist
        pass


def read_input_data(metadata: list,
                    dim: str = 'data_ensemble',
                    identifier_fmt: str = '{dataset}') -> 'xr.DataArray':
    """Loads data from metadata.

    Read the input data from the list of given data sets. `metadata` is a list
    of metadata containing the filenames to load. Only returns the given `variable`.
    The datasets are stacked along the `dim` dimension. Returns an xarray.DataArray."""

    data_arrays = []
    identifiers = []
    ancestors = []
    for info in metadata:
        filename = info['filename']
        variable = info['short_name']
        xrds = xr.open_dataset(filename)
        make_standard_calendar(xrds)
        data_arrays.append(xrds[variable])
        identifier = identifier_fmt.format(**info)
        identifiers.append(identifier)
        ancestors.append(filename)

    diagnostic = xr.concat(data_arrays, dim=dim)
    diagnostic[dim] = identifiers

    # Clean up unnecessary coordinate info
    redundant_dims = np.setdiff1d(diagnostic.coords, diagnostic.dims)
    diagnostic = diagnostic.drop(redundant_dims)

    return diagnostic, ancestors


def read_model_data(datasets: list) -> 'xr.DataArray':
    """Loads model data from list of metadata."""
    return read_input_data(datasets,
                           dim='model_ensemble',
                           identifier_fmt='{dataset}_{ensemble}_{exp}')


def read_observation_data(datasets: list) -> 'xr.DataArray':
    """Loads observation data from list of metadata."""
    return read_input_data(datasets,
                           dim='obs_ensemble',
                           identifier_fmt='{dataset}')


def aggregate_obs_data(data_array: 'xr.DataArray',
                       operator: str = 'median') -> 'xr.DataArray':
    """Reduce data array along ensemble dimension.

    Apply the operator to squeeze the ensemble dimension by applying
    the `operator` along the ensemble (`obs_ensemble`) dimension. Returns
    an xarray.Dataset squeezed to 1D."""

    if operator == 'median':
        return data_array.median(dim='obs_ensemble')
    else:
        raise ValueError(f'No such operator `{operator}`')


def area_weighted_mean(data_array: 'xr.DataArray') -> 'xr.DataArray':
    """Calculate area mean weighted by the latitude.

    Returns a data array consisting of N values,
    where N == number of ensemble members"""

    weights_lat = np.cos(np.radians(data_array.lat))
    means = data_array.weighted(weights_lat).mean(dim=['lat', 'lon'])

    return means


def distance_matrix(values: 'numpy.ndarray', weights: 'numpy.ndarray' = None) -> 'numpy.ndarray':
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
    """Calculate_independence.

    The independence is calculated as a distance matrix between
    the datasets defined in the `data_array`. Returned is a square matrix with
    where the number of elements along each edge equals
    the number of ensemble members."""
    weights = np.cos(np.radians(data_array.lat))
    weights, _ = xr.broadcast(weights, data_array)

    diff = xr.apply_ufunc(
        distance_matrix,
        data_array,
        weights,
        input_core_dims=[['model_ensemble', 'lat', 'lon'], ['model_ensemble', 'lat', 'lon']],
        output_core_dims=[['perfect_model_ensemble', 'model_ensemble']],
        )

    diff.name = f'd{data_array.name}'
    diff.attrs['short_name'] = data_array.name
    diff.attrs["units"] = data_array.units
    diff['perfect_model_ensemble'] = diff.model_ensemble.values

    return diff


def visualize_independence(independence, cfg, provenance_info):
    """Visualize_independence."""
    # import IPython; IPython.embed(); exit()
    variable = independence.short_name
    labels = [x.replace('_', '\n') for x in independence.model_ensemble.values]

    chart = sns.heatmap(
        independence,
        annot=True,
        fmt='.3g',
        linewidths=1,
        cmap="YlGn",
        xticklabels=labels,
        yticklabels=labels,
        cbar_kws={'label': f'Euclidean distance ({independence.units})'})
    chart.set_title(f'Distance matrix for {variable}')

    filename = get_plot_filename(f'independence_{variable}', cfg)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    caption = f'Euclidean distance matrix for variable {variable}'
    provenance_record = get_provenance_record(caption,
                                              ancestors=provenance_info)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance_record)

    print(f'Output stored as {filename}')


def calculate_performance(model_data: 'xr.DataArray',
                          obs_data: 'xr.DataArray') -> 'xr.DataArray':
    """Calculate_performance.

    Calculate the area weighted mean between the given ensemble of model data,
    and observation data. The observation data must have the
    ensemble dimension squeezed or reduced. Returns an xarray.DataArray
    containing the same number of values as members of `model_data`"""

    diff = model_data - obs_data

    performance = area_weighted_mean(diff**2)**0.5

    performance.name = f'd{model_data.name}'
    performance.attrs['short_name'] = model_data.name
    performance.attrs["units"] = model_data.units

    return performance


def visualize_performance(performance, cfg, provenance_info):
    """Visualize performance."""

    name = performance.name
    variable = performance.short_name
    units = performance.units

    df = performance.to_dataframe().reset_index()
    df.model_ensemble = df.model_ensemble.map(lambda x: x.replace('_', '\n'))

    ylabel = f'RMS error {variable} ({units})'

    chart = sns.barplot(x='model_ensemble', y=name, data=df)
    chart.set_title(f'Performance metric for {variable}')
    chart.set_ylabel(ylabel)
    chart.set_xlabel('')

    filename = get_plot_filename(f'performance_{variable}', cfg)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    caption = f'Performance data for variable {variable}'
    provenance_record = get_provenance_record(caption,
                                              ancestors=provenance_info)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance_record)

    print(f'Output stored as {filename}')


def calculate_weights(performance: 'xr.DataArray',
                      independence: 'xr.DataArray', sigma_performance: float,
                      sigma_independence: float) -> 'xr.DataArray':
    """Calculates the (NOT normalised) weights for each model N.
    Parameters
    ----------
    performance : array_like, shape (N,)
        Array specifying the model performance.
    independence : array_like, shape (N, N)
        Array specifying the model independence.
    sigma_performance : float
        Sigma value defining the form of the weighting function for the performance.
    sigma_independence : float
        Sigma value defining the form of the weighting function for the independence.
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

    weights.name = f'w{performance.short_name}'
    weights.attrs['short_name'] = performance.short_name
    weights.attrs["units"] = 'a.u.'

    return weights


def visualize_weights(weights, cfg, provenance_info):
    """Visualize weights."""

    name = weights.name
    variable = weights.short_name
    units = weights.units

    df = weights.to_dataframe().reset_index()
    df.model_ensemble = df.model_ensemble.map(lambda x: x.replace('_', '\n'))

    ylabel = f'Weights {variable} ({units})'

    chart = sns.barplot(x='model_ensemble', y=name, data=df)
    chart.set_title(f'Weights for {variable}')
    chart.set_ylabel(ylabel)
    chart.set_xlabel('')

    filename = get_plot_filename(f'weights_{variable}', cfg)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    caption = f'Weights for variable {variable}'
    provenance_record = get_provenance_record(caption,
                                              ancestors=provenance_info)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance_record)

    print(f'Output stored as {filename}')


def save_weights(s, cfg, provenance_info):
    """Save the weights to a `.yml` file."""

    filename = get_diagnostic_filename(f'weights', cfg, extension='yml')

    d = s.to_dict()

    with open(filename, 'w') as f:
        yaml.dump(d, stream=f)

    caption = f'Weights for all variables'
    provenance_record = get_provenance_record(caption,
                                              ancestors=provenance_info)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance_record)


def main(cfg):
    """Perform climwip weighting method."""

    print("\nBEGIN DIAGNOSTIC\n")

    observations = read_metadata(cfg, projects=['native6'])
    models = read_metadata(cfg, projects=['CMIP5'])

    variables = models.keys()

    weights_dict = {}
    ancestors = []

    for variable in variables:

        print(f'Reading model data for {variable}')
        datasets_model = models[variable]
        model_data, model_provenance = read_model_data(datasets_model)

        print(f'Reading observation data for {variable}')
        datasets_obs = observations[variable]
        obs_data, obs_provenance = read_observation_data(datasets_obs)
        obs_data = aggregate_obs_data(obs_data, operator='median')

        print(f'Calculating independence for {variable}')
        independence = calculate_independence(model_data)
        visualize_independence(independence, cfg, model_provenance)
        print(independence.values)
        print()

        print(f'Calculating performance for {variable}')
        performance = calculate_performance(model_data, obs_data)
        visualize_performance(performance, cfg, model_provenance + obs_provenance)
        print(performance.values)
        print()

        print(f'Calculating weights for {variable}')
        sigma_performance = cfg['shape_params'][variable]['sigma_d']
        sigma_independence = cfg['shape_params'][variable]['sigma_s']
        weights = calculate_weights(performance, independence,
                                    sigma_performance, sigma_independence)
        visualize_weights(weights, cfg, model_provenance + obs_provenance)
        print(weights.values)
        print()

        weights_dict[variable] = weights
        ancestors.extend(model_provenance)
        ancestors.extend(obs_provenance)

    df = xr.Dataset(weights_dict).to_dataframe()
    s = df.mean(axis=1)  # Average over variables
    save_weights(s, cfg, provenance_info=ancestors)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
