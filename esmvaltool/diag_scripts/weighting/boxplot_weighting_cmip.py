"""Plotting boxplots to show weighted ensemble results weights provided by
ClimWIP weighting scheme.

Lukas Brunner et al. section 2.4
https://iopscience.iop.org/article/10.1088/1748-9326/ab492f
"""
import logging
import os
from collections import defaultdict
from pathlib import Path
import numpy as np

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
from esmvaltool.diag_scripts.weighting.plot_utilities import (
    calculate_percentiles,
    read_metadata,
    read_weights,
)

logger = logging.getLogger(os.path.basename(__file__))

def _visualize_and_save_percentiles(data: 'xr.DataArray', weights: dict,
                                    models_fut: list, cfg: dict,
                                    ancestors: list):
    """Visualize data in boxplot and save percentiles."""
    # ensure weighhts are sorted the same way as data and convert to list
    #import ipdb; ipdb.set_trace()
    weights = [weights.sel(model_ensemble = member) for member in data.model_ensemble.values]
    #weights = np.ndarray.tolist(weights.values)

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

    models = read_metadata(cfg)['tas']
    model_data, model_data_files = read_model_data(models)

    # if a historical period is given calculate the change
    if 'tas_reference' in read_metadata(cfg).keys():
        models_hist = read_metadata(cfg)['tas_reference']
        model_data_hist, model_data_files_hist = read_model_data(models_hist)
        model_data_files += model_data_files_hist
        model_data = model_data - model_data_hist

    area_mean_change = area_weighted_mean(model_data)

    _visualize_and_save_percentiles(
        area_mean_change,
        weights,
        models,
        cfg,
        model_data_files,
    )


if __name__ == '__main__':
    with run_diagnostic() as configure:
        main(configure)
