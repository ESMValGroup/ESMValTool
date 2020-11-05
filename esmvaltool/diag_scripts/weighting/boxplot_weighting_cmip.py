"""Plotting boxplots to show weighted ensemble results weights provided by
ClimWIP weighting scheme.

Lukas Brunner et al. section 2.4
https://iopscience.iop.org/article/10.1088/1748-9326/ab492f
"""
import logging
import os
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import xarray as xr
from plot_utilities import boxplot
from climwip import (
    area_weighted_mean,
    get_diagnostic_filename,
    log_provenance,
    read_model_data,
)

from esmvaltool.diag_scripts.shared import get_plot_filename, run_diagnostic

logger = logging.getLogger(os.path.basename(__file__))


def read_weights(filename: str) -> dict:
    """Read a `.nc` file into a weights DataArray."""
    weights_ds = xr.open_dataset(filename)
    return weights_ds.to_dataframe().to_dict()['weight']


def read_metadata(cfg: dict) -> dict:
    """Read the metadata from the configure file."""
    datasets = defaultdict(list)

    metadata = cfg['input_data'].values()

    for item in metadata:
        variable = item['variable_group']

        datasets[variable].append(item)

    return datasets


def _visualize_and_save_percentiles(data: 'xr.DataArray', weights: dict,
                                    models_fut: list, cfg: dict,
                                    ancestors: list):
    """Visualize data in boxplot and save percentiles."""
    # ensure weighhts are sorted the same way as data and convert to list
    weights = [weights[member] for member in data.model_ensemble.values]

    figure, axes = plt.subplots(1, 1, figsize=(4, 10), dpi=300)
    figure.subplots_adjust(left=.1, right=.99, bottom=.22, top=.91)

    box_h1 = boxplot(axes,
                     0,
                     median=data,
                     mean=data,
                     box=data,
                     whisk=data,
                     width=.8,
                     color='darkgrey',
                     alpha=.3)

    box_h2 = boxplot(axes,
                     0,
                     median=data,
                     mean=data,
                     box=data,
                     whisk=data,
                     weights=weights,
                     width=.6,
                     color='tab:blue',
                     alpha=1)

    axes.set_ylabel('Temperature change relative to 1995-2015 [K]')
    axes.grid(axis='y')

    axes.set_xticklabels([])

    plt.ylim(0, 6)
    plt.legend((box_h1, box_h2), ('unweighted', 'weighted'))

    start = models_fut[0]['start_year']
    end = models_fut[0]['end_year']
    plt.title('Global temperature change by %s-%s' % (start, end))

    filename_plot = get_plot_filename('boxplot_temperature_change', cfg)
    figure.savefig(filename_plot, dpi=300, bbox_inches='tight')
    plt.close(figure)

    filename_data = get_diagnostic_filename('temperature_change_percentiles',
                                            cfg,
                                            extension='nc')
    data.to_netcdf(filename_data)

    caption = 'Temperature change boxplot'
    log_provenance(caption, filename_plot, cfg, ancestors)
    log_provenance(caption, filename_data, cfg, ancestors)


def main(cfg):
    """Plot weighted boxplot."""
    input_files = cfg['input_files']
    filename = cfg['weights']
    weights_path = Path(input_files[0]) / filename
    weights = read_weights(weights_path)

    models_fut = read_metadata(cfg)['tas_future']
    models_pres = read_metadata(cfg)['tas_reference']

    model_data1, model_data_files1 = read_model_data(models_fut)
    model_data2, model_data_files2 = read_model_data(models_pres)

    future_change = model_data1 - model_data2
    area_mean_change = area_weighted_mean(future_change)

    _visualize_and_save_percentiles(
        area_mean_change,
        weights,
        models_fut,
        cfg,
        [model_data_files1, model_data_files2],
    )


if __name__ == '__main__':
    with run_diagnostic() as configure:
        main(configure)
