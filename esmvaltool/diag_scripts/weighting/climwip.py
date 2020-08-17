"""
Implementation of step iii and iv of the climwip weighting scheme

Lukas Brunner et al. section 2.4
https://iopscience.iop.org/article/10.1088/1748-9326/ab492f
"""
from scipy import stats
from esmvaltool.diag_scripts.shared import run_diagnostic, select_metadata
import numpy as np
import xarray as xr
from scipy.spatial.distance import pdist, squareform
from collections import defaultdict


def weighted_distance_matrix(data_array: 'xr.DataArray'):
    pass


def area_weighted_mean(data_array: 'xr.DataArray') -> 'xr.DataArray':
    """Calculate area mean weighted by the latitude.

    Returns a data array consisting of N values,
    where N == number of ensemble members"""

    weights_lat = np.cos(np.radians(data_array.lat))
    means = data_array.weighted(weights_lat).mean(dim=['lat', 'lon'])

    return means


def distance_matrix(data_array: 'xr.DataArray') -> 'xr.DataArray':
    """Calculate the pairwise distance between model members.

    Takes a dataset with ensemble member/lon/lat. Flattens lon/lat
    into a single dimension. Calculates the distance between every
    ensemble member.

    Returns 2D NxN array, where N == number of ensemble members.
    """
    n_members = data_array.shape[0]

    data_array = data_array.reshape(n_members, -1)

    # pdist does not work with NaN
    idx = np.where(np.all(np.isfinite(data_array), axis=0))[0]
    data_array = data_array[:, idx]

    d_matrix = squareform(pdist(data_array, metric='euclidean'))

    return d_matrix


def calculate_independence(data_array: 'xr.DataArray') -> 'xr.DataArray':
    """Calculate_independence.

    The independence is calculated as a distance matrix between
    the datasets defined in the `data_array`. Returned is a square matrix with
    where the number of elements along each edge equals
    the number of ensemble members."""
    # TODO: use weighted_distance_matrix

    diff = xr.apply_ufunc(
        distance_matrix,
        data_array,
        input_core_dims=[['model_ensemble', 'lat', 'lon']],
        output_core_dims=[['perfect_model_ensemble', 'model_ensemble']])

    diff.name = 'data'
    # diff = diff.expand_dims({'diagnostic': [idx]})

    return diff


def visualize_independence(independence):
    """visualize_independence."""
    # TODO: complete this function
    returned_stuff = None
    return returned_stuff


def calculate_performance(model_data: 'xr.DataArray',
                          obs_data: 'xr.DataArray') -> 'xr.DataArray':
    """Calculate_performance.

    Calculate the area weighted mean between the given ensemble of model data,
    and observation data. The observation data must have the
    ensemble dimension squeezed or reduced. Returns an xarray.DataArray
    containing the same number of values as members of `model_data`"""

    diff = model_data - obs_data

    performance = area_weighted_mean(diff**2)**0.5

    return performance


def visualize_performance(performance):
    """visualize_performance."""

    # TODO: complete this function
    returned_stuff = None
    return returned_stuff


def calculate_weights(performance: 'xr.DataArray',
                      independence: 'xr.DataArray',
                      sigma_performance: float,
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

    return weights


def visualize_weighting():
    """visualize_weighting."""

    # TODO: complete this function
    returned_stuff = None
    return returned_stuff


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


def read_input_data(datasets: list,
                    dim: str = 'data_ensemble') -> 'xr.DataArray':
    """Loads data from metadata.

    Read the input data from the list of given data sets. `datasets` is a list
    of metadata containing the filenames to load. Only returns the given `variable`.
    The datasets are stacked along the `dim` dimension. Returns an xarray.DataArray."""

    to_concat = []

    for dataset in datasets:
        fn = dataset['filename']
        variable = dataset['short_name']
        xarr = xr.open_dataset(fn)
        to_concat.append(xarr[variable])

    diagnostic = xr.concat(to_concat, dim=dim)

    return diagnostic


def read_model_data(datasets: list) -> 'xr.DataArray':
    """Loads model data from list of metadata."""
    return read_input_data(datasets, dim='model_ensemble')


def read_observation_data(datasets: list) -> 'xr.DataArray':
    """Loads observation data from list of metadata."""
    return read_input_data(datasets, dim='obs_ensemble')


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


def main(cfg):
    """Perform climwip weighting method."""

    print("\nBEGIN DIAGNOSTIC\n")

    observations = read_metadata(cfg, projects=['native6'])
    models = read_metadata(cfg, projects=['CMIP5'])

    variables = models.keys()

    model_data_dict = {}
    obs_data_dict = {}
    independence_dict = {}
    performance_dict = {}

    for variable, datasets in models.items():
        print(f'Reading model data for {variable}')

        model_data = read_model_data(datasets)
        model_data_dict[variable] = model_data

    for variable, datasets in observations.items():
        print(f'Reading observation data for {variable}')

        obs_data = read_observation_data(datasets)
        obs_data = aggregate_obs_data(obs_data, operator='median')
        obs_data_dict[variable] = obs_data

    for variable, model_data in model_data_dict.items():
        print(f'Calculating independence for {variable}')

        independence = calculate_independence(model_data)
        independence_dict[variable] = independence

        visualize_independence(independence)  # TODO

        print(independence.values)
        print()

    for variable in variables:
        print(f'Calculating performance for {variable}')

        model_data = model_data_dict[variable]
        obs_data = obs_data_dict[variable]

        performance = calculate_performance(model_data, obs_data)
        performance_dict[variable] = performance

        visualize_performance(performance)  # TODO

        print(performance.values)
        print()

    for variable in variables:
        print(f'Calculating weights for {variable}')

        performance = performance_dict[variable]
        independence = independence_dict[variable]

        sigma_performance = cfg['shape_params'][variable]['sigma_d']
        sigma_independence = cfg['shape_params'][variable]['sigma_s']

        weights = calculate_weights(performance, independence,
                                    sigma_performance, sigma_independence)

        print(weights.values)
        print()

        # visualize_weighting(weights)

    # independences = xr.concat(independences, dim='diagnostic')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
