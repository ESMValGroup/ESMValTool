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
                                    central_estimate: 'xr.DataArray',
                                    central_estimate_weighted: 'xr.DataArray',
                                    uncertainty_range: 'xr.DataArray',
                                    uncertainty_range_weighted: 'xr.DataArray',
                                    cfg: dict, ancestors: list):
    """Visualize weighted temperature."""
    figure, axes = plt.subplots(dpi=300)

    def plot_shaded(xrange, upper, lower, color, **kwargs):
        axes.fill_between(
            xrange.data,
            upper.data,
            lower.data,
            facecolor=color,
            edgecolor='none',
            alpha=0.3,
            zorder=100,
            **kwargs,
        )

    def plot_line(xrange, central, **kwargs):
        axes.plot(
            xrange,
            central,
            zorder=1000,
            **kwargs,
        )

    color_non_weighted = 'red'
    color_weighted = 'green'
    color_data = 'gray'
    central_string = cfg['settings'].get('central_estimate', 50)
    range_string = '{}-{}perc'.format(cfg['settings'].get('lower_bound', 25),
                                      cfg['settings'].get('upper_bound', 75))
    if not isinstance(central_string, str):
        central_string = f'{central_string}perc'

    plot_line(central_estimate.time,
              central_estimate,
              color=color_non_weighted,
              label='Non-weighted {}'.format(central_string))
    plot_shaded(uncertainty_range.time,
                uncertainty_range[:, 0],
                uncertainty_range[:, 1],
                color=color_non_weighted,
                label=f'Non-weighted {range_string} range')

    plot_line(central_estimate_weighted.time,
              central_estimate_weighted,
              color=color_weighted,
              label='Weighted {}'.format(central_string))
    plot_shaded(
        uncertainty_range_weighted.time,
        uncertainty_range_weighted[:, 0],
        uncertainty_range_weighted[:, 1],
        color=color_weighted,
        label='Weighted {} range'.format(range_string),
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

    start_year = cfg['settings']['start_year']
    end_year = cfg['settings']['end_year']
    caption = f'Temperature anomaly relative to {start_year}-{end_year}'
    plt.title(caption)
    plt.xlabel('Year')
    plt.ylabel(r'Temperature anomaly $\degree$C')

    filename_plot = get_plot_filename('temperature_anomaly_graph', cfg)
    figure.savefig(filename_plot, dpi=300, bbox_inches='tight')
    plt.close(figure)

    filename_data = get_diagnostic_filename('temperature_anomalies',
                                            cfg,
                                            extension='nc')
    temperature.to_netcdf(filename_data)

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

    settings = cfg['settings']
    central_estimate_var = settings.get('central_estimate', 50)
    if isinstance(central_estimate_var, (int, float)):
        central_estimate = calculate_percentiles(
            model_data,
            np.array([central_estimate_var]),
        )
        central_estimate_weighted = calculate_percentiles(
            model_data,
            np.array([central_estimate_var]),
            weights=weights,
        )
    elif central_estimate_var == 'mean':
        central_estimate = model_data.mean('model_ensemble')
        central_estimate_weighted = model_data.weighted(weights).mean(
            'model_ensemble')

    percentiles = np.array(
        [settings.get('lower_bound', 25),
         settings.get('upper_bound', 75)])

    uncertainty_range = calculate_percentiles(
        model_data,
        percentiles,
    )

    uncertainty_range_weighted = calculate_percentiles(
        model_data,
        percentiles,
        weights=weights,
    )

    visualize_and_save_temperatures(
        model_data,
        central_estimate,
        central_estimate_weighted,
        uncertainty_range,
        uncertainty_range_weighted,
        cfg,
        model_data_files,
    )


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
