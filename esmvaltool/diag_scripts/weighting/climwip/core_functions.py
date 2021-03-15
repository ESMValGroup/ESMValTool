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


def calculate_model_distances(
        data_array: 'xr.DataArray',
        dimension: str = 'model_ensemble_reference') -> 'xr.DataArray':
    """Calculate pair-wise distances between all values in data_array.

    Distances are calculated as the area weighted euclidean distance
    between each pair of models in data_array. Returned is a square matrix
    with where the number of elements along each edge equals the number
    of ensemble members.

    Parameters
    ----------
    data_array : array_like, shape (N,...)
        Array of (2 dimensional) model fields.
    dimension : string
        Name of the newly created reference dimension (default:
        'model_ensemble_reference'. Must not be equal to the existing
        model dimension ('model_ensemble')!

    Returns
    -------
    distances : array_like, shape (N, N)
        Symmetric matrix of pairwise model distances.
    """
    assert dimension != 'model_ensemble', f'{dimension} != "model_ensemble"'
    weights = np.cos(np.radians(data_array.lat))
    weights, _ = xr.broadcast(weights, data_array)

    diff = xr.apply_ufunc(
        distance_matrix,
        data_array,
        weights,
        input_core_dims=[['model_ensemble', 'lat', 'lon'],
                         ['model_ensemble', 'lat', 'lon']],
        output_core_dims=[[dimension, 'model_ensemble']],
    )

    diff.name = f'd{data_array.name}'
    diff.attrs['variable_group'] = data_array.name
    diff.attrs["units"] = data_array.units
    diff[dimension] = diff.model_ensemble.values

    return diff


def compute_overall_mean(dataset: 'xr.Dataset',
                         weights: dict) -> 'xr.DataArray':
    """Normalize all variables in a dataset and return their weighted mean.

    Relative weights for each variable group are passed via the recipe.
    """
    normalized = dataset / dataset.median()

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
        dataset: Union['xr.DataArray', None],
        dimensions: Union[str, list] = 'model_ensemble',
) -> (Union['xr.DataArray', None], dict):
    """Combine ensemble members of the same model.

    Parameters
    ----------
    dataset : None or data_array, shape (N,) or (N, N)
        A vector containing model-observations distances or a matrix containing
        model-model distances.
    dimensions : string or list of up to two strings
        Spezifies the dimensions along which ensemble members are combined.

    Returns
    -------
    dataset : None or data_array, shape (M,), (M, L) with M, L <= N
        data_array where ensemble members along the given dimensions are
        combined by averaging.
    groups : dict of form {string: list}
        Dictionary mapping the combined model names (keys) to the original
        ensemble member names (values).
    """
    if isinstance(dimensions, str):
        dimensions = [dimensions]
    assert len(
        dimensions) <= 2, 'dimensions can contain a maximum of two strings'

    if dataset is None:
        return None, {}

    groups = defaultdict(list)
    models = []
    for name in dataset['model_ensemble'].values:
        model = name.split('_')[0]
        groups[model].append(name)
        models.append(model)

    for dimension in dimensions:
        if dimension in dataset.dims:
            model = xr.DataArray(models, dims=dimension)
            dataset = dataset.groupby(model).mean(keep_attrs=True).rename(
                {'group': dimension})

    if len(dimensions) == 2:
        # need to set the diagonal elements back to zero after averaging
        dataset.values[np.diag_indices(dataset['model_ensemble'].size)] = 0

    return dataset, groups


def calculate_weights_data(
        performance: Union['np.array', None],
        independence: Union['np.array', None],
        performance_sigma: Union[float, None],
        independence_sigma: Union[float, None]) -> 'np.array':
    """Calculate normalized weights for each model N.

    Parameters
    ----------
    performance : array_like, shape (N,) or None
        Array specifying the model performance. None is mutually exclusive
        with independence being None. Single values in performance can be
        nan, then they will be excluded from the independence calculation as
        well (used for the perfect model test).
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
    numerator = 1
    not_nan = True
    denominator = 1

    if performance is not None:
        numerator = np.exp(-((performance / performance_sigma)**2))
        # nans in the performance vector indicate models to be excluded
        not_nan = np.isfinite(performance)
    if independence is not None:
        # don't consider nan models for independence of other models!
        exp = np.exp(-((independence[:, not_nan] / independence_sigma)**2))
        # Note diagonal = exp(0) = 1, thus this is equal to 1 + sum(i!=j)
        denominator = exp.sum(axis=1)

    weights = numerator / denominator
    weights /= weights.sum(where=not_nan)

    return weights


def calculate_weights(
        performance: Union['xr.DataArray', None],
        independence: Union['xr.DataArray', None],
        performance_sigma: Union[float, None],
        independence_sigma: Union[float, None]) -> 'xr.DataArray':
    """Xarray wrapper for calculate_weights_data."""
    performance_core_dims = [] if performance is None else ['model_ensemble']
    independence_core_dims = [] if independence is None else [
        'model_ensemble', 'model_ensemble_reference'
    ]

    weights = xr.apply_ufunc(
        calculate_weights_data,
        performance,
        independence,
        performance_sigma,
        independence_sigma,
        input_core_dims=[
            performance_core_dims, independence_core_dims, [], []
        ],
        output_core_dims=[['model_ensemble']],
        vectorize=True,
    )

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

    # remove nans
    not_nan = np.where((np.isfinite(values) & np.isfinite(weights)))[0]
    values = values[not_nan]
    weights = weights[not_nan]

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
