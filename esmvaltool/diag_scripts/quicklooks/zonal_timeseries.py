#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot zonal time series.

Description
-----------
Plot zonal time series. In quicklook mode, this diagnostic plots the
concatenated files.

Author
------
Lisa Bock (DLR, Germany)

Project
-------
CMIP6-DICAD

Configuration options in recipe
-------------------------------
latitude_range : list of float, optional
    Range for the latitde axis in the plots.
levels : list of float, optional
    Values for contour levels.
time_range : list of float, optional
    Range for the time axis in the plots.

"""

import logging
import os

import matplotlib.pyplot as plt

from esmvaltool.diag_scripts.quicklooks import (
    cube_time_to_float, get_grouped_input_data, get_provenance_record,
    load_cube, save_plot, set_plot_title, write_provenance)
from esmvaltool.diag_scripts.shared import (get_diagnostic_filename, io,
                                            run_diagnostic)

logger = logging.getLogger(os.path.basename(__file__))


def plot_single_dataset(cfg, dataset):
    """Plot zonal mean time series for single dataset.

    Parameters
    ----------
    cfg: dict
        Global configuration dictionary.
    dataset: list
        Metadata dictionary for desired dataset.

    """
    cube = load_cube(dataset, ['time', 'latitude'])
    if cube is None:
        return

    # Provenance
    caption = (
        f"Time series plot of zonal mean {dataset['short_name']} for dataset "
        f"{dataset['dataset']}.")
    provenance_record = get_provenance_record(caption, [dataset['filename']])

    # Transpose cube
    data = cube.data.transpose()
    times = cube_time_to_float(cube)

    # Plot
    plt.contourf(times,
                 cube.coord('latitude').points,
                 data,
                 extend='both',
                 levels=cfg.get('levels'))
    title = ' '.join([dataset['dataset'], dataset['long_name']])
    set_plot_title(title)
    plt.xlabel('year')
    plt.ylabel('latitude')
    if 'time_range' in cfg:
        plt.xlim(cfg['time_range'][0], cfg['time_range'][1])
    if 'latitude_range' in cfg:
        plt.ylim(cfg['latitude_range'][0], cfg['latitude_range'][1])

    # Make a colorbar for the ContourSet returned by the contourf call
    colorbar = plt.colorbar(orientation='horizontal', aspect=30)
    colorbar.set_label(f"{dataset['short_name']} / {cube.units}")

    # Save plot
    basename = '_'.join([
        dataset['short_name'],
        dataset['dataset'],
        'zonal_timeseries',
    ])
    plot_provenance = save_plot(basename, 'zonal', cfg)
    provenance_record.update(plot_provenance)

    # Write netcdf file and provenance
    netcdf_path = get_diagnostic_filename(basename, cfg)
    io.iris_save(cube, netcdf_path)
    write_provenance(netcdf_path, provenance_record, cfg)


def main(cfg):
    """Load config file and metadata, then pass them the plot making tools.

    Parameters
    ----------
    cfg: dict
        Global configuration dictionary.

    """
    grouped_data = get_grouped_input_data(cfg)

    # Iterate over data and plot
    for (short_name, datasets) in grouped_data.items():
        logger.info("Processing variable '%s'", short_name)
        for dataset in datasets:
            logger.info("Processing '%s'", dataset['filename'])

            # Create zonal time series for every dataset
            plot_single_dataset(cfg, dataset)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
