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
import pandas as pd
from pathlib import Path
import yaml

from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            get_plot_filename, run_diagnostic,
                                            select_metadata)

from climwip import read_metadata, read_model_data

def read_weights(filename: str):
    """Read a `.yml` file into a Pandas Dataframe."""
    with open(filename, 'r') as f:
        weights = yaml.safe_load(f)
    return weights


def weighted_quantile(values, quantiles, sample_weight=None,
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.

    source: https://stackoverflow.com/a/29677616/6012085

    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)


def main(cfg):
    # TODO: Finish this script

    # breakpoint()

    input_files = cfg['input_files']
    filename = cfg['weights']
    weights_path = Path(input_files[0]) / filename
    weights = read_weights(weights_path)

    models = read_metadata(cfg, projects=['CMIP5'])['tas']
    model_data, model_provenance = read_model_data(models)

    w = [weights[x] for x in model_data.model_ensemble.values]

    percentiles = np.array([25, 75])

    weighted_percentiles = xr.apply_ufunc(
        weighted_quantile,
        model_data,
        input_core_dims = [['model_ensemble']],
        output_core_dims= [['percentiles']],
        kwargs={
            'sample_weight': w,
            'quantiles': percentiles/100
            },
        vectorize=True
        )

    weighted_percentiles['percentiles'] = percentiles


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)

