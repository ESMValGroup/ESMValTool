"""Implementation of the climwip weighting scheme.

Lukas Brunner et al. section 2.4
https://iopscience.iop.org/article/10.1088/1748-9326/ab492f
"""
from collections import defaultdict
import logging
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from climwip import log_provenance, read_model_data
from climwip import get_diagnostic_filename

from esmvaltool.diag_scripts.shared import get_plot_filename, run_diagnostic

logger = logging.getLogger(os.path.basename(__file__))


def read_weights(filename: str) -> dict:
    """Read a `.nc` file into a weights DataArray."""
    weights_ds = xr.open_dataset(filename)
    return weights_ds.to_dataframe().to_dict()['weight']


def read_metadata(cfg: dict) -> dict:
    """Read the metadata from the config file."""
    datasets = defaultdict(list)

    metadata = cfg['input_data'].values()

    for item in metadata:
        variable = item['short_name']

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


def visualize_and_save_temperatures(temperature: 'xr.DataArray',
                                    iqr: 'xr.DataArray',
                                    iqr_weighted: 'xr.DataArray',
                                    cfg: dict,
                                    ancestors: list):
    """Visualize weighted temperature."""
    figure, axes = plt.subplots(dpi=300)

    def plot_shaded(xrange, upper, lower, color, **kwargs):
        axes.fill_between(
            xrange,
            upper,
            lower,
            facecolor=color,
            edgecolor='none',
            alpha=0.5,
            zorder=100,
            **kwargs,
        )

    color_non_weighted = 'red'
    color_weighted = 'green'
    color_data = 'gray'

    plot_shaded(iqr.time,
                iqr.data[:, 0],
                iqr.data[:, 1],
                color=color_non_weighted,
                label='Non-weighted inter-quartile range')

    plot_shaded(
        iqr_weighted.time,
        iqr_weighted.data[:, 0],
        iqr_weighted.data[:, 1],
        color=color_weighted,
        label='Weighted inter-quartile range',
    )

    for temp in temperature.data:
        axes.plot(temperature.time,
                  temp,
                  color=color_data,
                  lw=0.5,
                  alpha=0.5,
                  zorder=1,
                  label='Ensemble members')

    # Fix duplicate labels
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))  # dict removes dupes
    axes.legend(by_label.values(), by_label.keys())

    plt.title('Temperature anomaly relative to 1981-2010')
    plt.xlabel('Year')
    plt.ylabel(r'Temperature anomaly $\degree$C')

    filename_plot = get_plot_filename('temperature_anomaly_graph', cfg)
    figure.savefig(filename_plot, dpi=300, bbox_inches='tight')
    plt.close(figure)

    filename_data = get_diagnostic_filename(
        'temperature_anomalies', cfg, extension='nc')
    temperature.to_netcdf(filename_data)

    caption = 'Temperature anomaly relative to 1981-2010'
    log_provenance(caption, filename_plot, cfg, ancestors)
    log_provenance(caption, filename_data, cfg, ancestors)


def calculate_percentiles(data: 'xr.DataArray',
                          percentiles: list,
                          weights: dict = None):
    """Calculate (weighted) percentiles.

    Calculate the (weighted) percentiles for the given data.

    Percentiles is a list of values between 0 and 100.

    Weights is a dictionary where the keys represent the model members,
    and the values the weights.
    If `weights` is not specified, the non-weighted percentiles are calculated.

    Returns a DataArray with 'percentiles' as the dimension.
    """
    if weights:
        weights = [weights[member] for member in data.model_ensemble.values]

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


def main(cfg):
    """Plot weighted temperature graph."""
    input_files = cfg['input_files']
    filename = cfg['weights']
    weights_path = Path(input_files[0]) / filename
    weights = read_weights(weights_path)

    models = read_metadata(cfg)['tas']

    model_data, model_data_files = read_model_data(models)

    percentiles = np.array([25, 75])

    iqr = calculate_percentiles(
        model_data,
        percentiles,
    )

    iqr_weighted = calculate_percentiles(
        model_data,
        percentiles,
        weights=weights,
    )

    visualize_and_save_temperatures(
        model_data,
        iqr,
        iqr_weighted,
        cfg,
        model_data_files,
    )


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
