"""A collection of utility functions for dealing with weights."""
from collections import defaultdict

import numpy as np
import xarray as xr


def read_weights(filename: str) -> dict:
    """Read a `.nc` file into a weights DataArray."""
    weights_ds = xr.open_dataset(filename)
    return weights_ds['weight']


def read_metadata(cfg: dict, groupby: str = 'variable_group') -> dict:
    """Read the metadata from the config file."""
    datasets = defaultdict(list)

    metadata = cfg['input_data'].values()

    for item in metadata:
        variable = item[groupby]

        datasets[variable].append(item)

    return datasets


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


def calculate_percentiles(data: 'xr.DataArray',
                          percentiles: list,
                          weights: dict = None) -> 'xr.DataArray':
    """Calculate (weighted) percentiles.

    Calculate the (weighted) percentiles for the given data.

    Percentiles is a list of values between 0 and 100.

    The `model_ensemble` dimension in weights has to contain at
    least the same elements as in data.
    If `weights` is not specified, the non-weighted percentiles are calculated.

    Returns a DataArray with 'percentiles' as the dimension.
    """
    if weights is not None:
        weights = weights.sel(model_ensemble=data.model_ensemble)

    output = xr.apply_ufunc(weighted_quantile,
                            data,
                            input_core_dims=[['model_ensemble']],
                            output_core_dims=[['percentiles']],
                            kwargs={
                                'weights': weights,
                                'quantiles': percentiles / 100
                            },
                            vectorize=True)

    output['percentiles'] = percentiles

    return output
