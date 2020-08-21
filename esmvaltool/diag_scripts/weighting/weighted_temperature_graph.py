"""
Implementation of the climwip weighting scheme

Lukas Brunner et al. section 2.4
https://iopscience.iop.org/article/10.1088/1748-9326/ab492f
"""
import logging
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import yaml

from climwip import log_provenance, read_metadata, read_model_data
from esmvaltool.diag_scripts.shared import get_plot_filename, run_diagnostic

logger = logging.getLogger(os.path.basename(__file__))


def read_weights(filename: str):
    """Read a `.yml` file into a Pandas Dataframe."""
    with open(filename, 'r') as f:
        weights = yaml.safe_load(f)
    return weights


def weighted_quantile(values,
                      quantiles,
                      sample_weight=None,
                      values_sorted=False,
                      old_style=False):
    """Very close to numpy.percentile, but supports weights.

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


def visualize_temperature_graph(
    temperature,
    iqr,
    iqr_weighted,
    cfg,
    provenance_info,
):
    """Visualize weighted temperature."""

    fig, ax = plt.subplots(dpi=300)

    def plot_shaded(t, upper, lower, color, **kwargs):
        ax.fill_between(
            t,
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

    for t in temperature.data:
        ax.plot(temperature.time,
                t,
                color=color_data,
                lw=0.5,
                alpha=0.5,
                zorder=1,
                label='Ensemble members')

    # Fix duplicate labels
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))  # dict removes dupes
    ax.legend(by_label.values(), by_label.keys())

    plt.title('Temperature anomaly relative to 1981-2010')
    plt.xlabel('Year')
    plt.ylabel('Temperature anomaly $\degree$C')

    filename = get_plot_filename(f'temperature_anomaly_graph', cfg)
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    caption = f'Temperature anomaly relative to 1981-2010'

    log_provenance(caption, filename, cfg, provenance_info)

    logger.info(f'Temperature anomaly plot stored as {filename}')


def calculate_percentiles(data: 'xr.DataArray',
                          percentiles: list,
                          weights: dict = None):
    """Calculate (weighted) percentiles.

    Calculate the (weighted) percentiles for the given data.

    Percentiles is a list of values between 0 and 100.

    Weights is a dictionary where the keys represent the model members, and the values the weights.
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
                                'sample_weight': weights,
                                'quantiles': percentiles / 100
                            },
                            vectorize=True)

    output['percentiles'] = percentiles

    return output


def main(cfg):
    input_files = cfg['input_files']
    filename = cfg['weights']
    weights_path = Path(input_files[0]) / filename
    weights = read_weights(weights_path)

    models = read_metadata(cfg, projects=['CMIP5'])['tas']
    model_data, model_provenance = read_model_data(models)

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

    visualize_temperature_graph(
        model_data,
        iqr,
        iqr_weighted,
        cfg,
        model_provenance,
    )


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
