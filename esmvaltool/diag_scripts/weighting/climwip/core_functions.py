"""A collection of core functions."""
import logging
import os
from collections import defaultdict
from typing import Union

import numpy as np
import xarray as xr
from scipy.spatial.distance import pdist, squareform

logger = logging.getLogger(os.path.basename(__file__))


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


def compute_overall_mean(dataset: 'xr.Dataset',
                         weights: dict) -> 'xr.DataArray':
    """Normalize all variables in a dataset and return their weighted mean.

    Relative weights for each variable group are passed via the recipe.
    """
    if 'perfect_model_ensemble' in dataset.dims:
        median_dim = ['perfect_model_ensemble', 'model_ensemble']
    else:
        median_dim = 'model_ensemble'

    normalized = dataset / dataset.median(dim=median_dim)

    weights_selected = xr.DataArray(
        [weights[variable_group] for variable_group in dataset],
        coords={'variable_group': list(dataset)},
        dims='variable_group')
    overall_mean = normalized.to_array(
        dim='variable_group').weighted(weights_selected).mean('variable_group')
    overall_mean.name = 'overall_mean'
    overall_mean.attrs['variable_group'] = 'overall_mean'
    overall_mean.attrs['units'] = '1'
    return overall_mean


def combine_ensemble_members(
        dataset: Union['xr.DataArray', None]) -> (
            Union['xr.DataArray', None], dict):
    """Combine ensemble members of the same model."""
    if dataset is None:
        return None, {}

    groups = defaultdict(list)
    models = []
    for name in dataset['model_ensemble'].values:
        model = name.split('_')[0]
        groups[model].append(name)
        models.append(model)

    for dimn in ['model_ensemble', 'perfect_model_ensemble']:
        if dimn in dataset.dims:
            model = xr.DataArray(models, dims=dimn)
            dataset = dataset.groupby(model).mean(keep_attrs=True).rename(
                {'group': dimn})

    if 'perfect_model_ensemble' in dataset.dims:
        # need to set the diagonal elements back to zero after averaging
        dataset.values[np.diag_indices(dataset['model_ensemble'].size)] = 0

    return dataset, groups


def calculate_weights(
        performance: Union['xr.DataArray', None],
        independence: Union['xr.DataArray', None],
        performance_sigma: Union[float, None],
        independence_sigma: Union[float, None]) -> 'xr.DataArray':
    """Calculate normalized weights for each model N.

    Parameters
    ----------
    performance : array_like, shape (N,) or None
        Array specifying the model performance. None is mutually exclusive
        with independence being None.
    independence : array_like, shape (N, N) or None
        Array specifying the model independence. None is mutually exclusive
        with performance being None.
    performance_sigma : float or None
        Sigma value defining the form of the weighting function
        for the performance. Can be one only if performance is also None.
    independence_sigma : float or None
        Sigma value defining the form of the weighting function
            for the independence. Can be one only if independence is also None.

    Returns
    -------
    weights : ndarray, shape (N,)
    """
    if performance is not None:
        numerator = np.exp(-((performance / performance_sigma)**2))
    else:
        numerator = 1
    if independence is not None:
        exp = np.exp(-((independence / independence_sigma)**2))
        # Note diagonal = exp(0) = 1, thus this is equal to 1 + sum(i!=j)
        denominator = exp.sum('perfect_model_ensemble')
    else:
        denominator = 1

    weights = numerator / denominator

    # Normalize weights
    weights /= weights.sum()

    weights.name = 'weight'
    weights.attrs['variable_group'] = 'weight'  # used in barplot
    weights.attrs['units'] = '1'

    return weights


def weighted_quantile(values: list,
                      quantiles: list,
                      weights: list = None) -> 'np.array':
    """Calculate weighted quantiles.

    Analogous to np.quantile, but supports weights.

    Based on: https://stackoverflow.com/a/29677616/6012085

    Parameters
    ----------
    values: array_like
        List of input values.
    quantiles: array_like
        List of quantiles between 0.0 and 1.0.
    weights: array_like
        List with same length as `values` containing the weights.

    Returns
    -------
    np.array
        Numpy array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if weights is None:
        weights = np.ones(len(values))
    weights = np.array(weights)

    if not np.all((quantiles >= 0) & (quantiles <= 1)):
        raise ValueError('Quantiles should be between 0.0 and 1.0')

    idx = np.argsort(values)
    values = values[idx]
    weights = weights[idx]

    weighted_quantiles = np.cumsum(weights) - 0.5 * weights

    # Cast weighted quantiles to 0-1 To be consistent with np.quantile
    min_val = weighted_quantiles.min()
    max_val = weighted_quantiles.max()
    weighted_quantiles = (weighted_quantiles - min_val) / max_val

    return np.interp(quantiles, weighted_quantiles, values)
