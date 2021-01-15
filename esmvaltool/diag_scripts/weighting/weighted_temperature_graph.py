"""Implementation of the climwip weighting scheme.

Lukas Brunner et al. section 2.4
https://iopscience.iop.org/article/10.1088/1748-9326/ab492f
"""
import logging
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from climwip.io_functions import log_provenance, read_model_data

from esmvaltool.diag_scripts.shared import (
    get_diagnostic_filename,
    get_plot_filename,
    run_diagnostic,
)
from esmvaltool.diag_scripts.weighting.plot_utilities import (
    calculate_percentiles,
    read_metadata,
    read_weights,
)

logger = logging.getLogger(os.path.basename(__file__))


def visualize_and_save_temperatures(temperature: 'xr.DataArray',
                                    iqr: 'xr.DataArray',
                                    iqr_weighted: 'xr.DataArray', cfg: dict,
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

    filename_data = get_diagnostic_filename('temperature_anomalies',
                                            cfg,
                                            extension='nc')
    temperature.to_netcdf(filename_data)

    caption = 'Temperature anomaly relative to 1981-2010'
    log_provenance(caption, filename_plot, cfg, ancestors)
    log_provenance(caption, filename_data, cfg, ancestors)


def main(cfg):
    """Plot weighted temperature graph."""
    input_files = cfg['input_files']
    filename = cfg['weights']
    weights_path = Path(input_files[0]) / filename
    weights = read_weights(weights_path)

    models = read_metadata(cfg, 'short_name')['tas']
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
