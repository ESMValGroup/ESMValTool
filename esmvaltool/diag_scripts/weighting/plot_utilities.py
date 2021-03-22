"""A collection of utility functions for dealing with weights."""
from collections import defaultdict

import xarray as xr

from climwip.core_functions import weighted_quantile


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
