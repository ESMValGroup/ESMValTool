#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot global time series.

Description
-----------
Plot global time series. In quicklook mode, this diagnostic plots the
concatenated files.

Author
------
Lisa Bock (DLR, Germany)

Project
-------
CMIP6-DICAD

Configuration options in recipe
-------------------------------
multi_dataset_plot : bool, optional (default: False)
    If given, plot all given datasets for every variable in one plot.
read_all_available_datasets : bool, optional (default: False)
    If `True`, read all available datasets given to the diagnostic script. If
    `False`, only process variables given in the respective variable section of
    the diagnostic block (this settings only affects the quicklook mode).
time_range : list of float, optional
    Range for the time axis in the plots.
y_range : list of float, optional
    Range for the variable axis.

"""

import logging
import os

import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.quicklooks import (
    cube_time_to_float, get_grouped_input_data, get_provenance_record,
    load_cube, save_plot, set_plot_title, write_provenance)
from esmvaltool.diag_scripts.shared import (get_diagnostic_filename, io,
                                            run_diagnostic)

logger = logging.getLogger(os.path.basename(__file__))


def global_time_series_plot(cube, **kwargs):
    """Create a time series plot from the cube.

    Note that this function simple does the plotting, it does not save the
    image or do any of the complex work. This function also takes and of the
    key word arguments accepted by the matplotlib.pyplot.plot function.
    These arguments are typically, color, linewidth, linestyle, etc...

    If there's only one datapoint in the cube, it is plotted as a
    horizontal line.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube.

    """
    cube_data = np.ma.array(cube.data)
    if len(cube_data.compressed()) == 1:
        plt.axhline(cube_data.compressed(), **kwargs)
        return
    times = cube_time_to_float(cube)
    plt.plot(times, cube_data, **kwargs)


def plot_single_dataset(cfg, dataset):
    """Plot global mean time series for single dataset.

    Parameters
    ----------
    cfg: dict
        Global configuration dictionary.
    dataset: list
        Metadata dictionary for desired dataset.

    """
    cube = load_cube(dataset, ['time'])
    if cube is None:
        return

    # Provenance
    caption = (
        f"Time series plot of global mean {dataset['short_name']} for dataset "
        f"{dataset['dataset']}.")
    provenance_record = get_provenance_record(caption, [dataset['filename']])

    # Plot
    global_time_series_plot(cube)
    title = ' '.join([dataset['dataset'], dataset['long_name']])
    set_plot_title(title)
    plt.xlabel('year')
    plt.ylabel(f"{dataset['short_name']} / {cube.units}")
    if 'time_range' in cfg:
        plt.xlim(cfg['time_range'][0], cfg['time_range'][1])
    if 'y_range' in cfg:
        plt.ylim(cfg['y_range'][0], cfg['y_range'][1])

    # Save plot
    basename = '_'.join([
        dataset['short_name'],
        dataset['dataset'],
        'global_timeseries',
    ])
    plot_provenance = save_plot(basename, 'times', cfg)
    provenance_record.update(plot_provenance)

    # Write netcdf file and provenance
    netcdf_path = get_diagnostic_filename(basename, cfg)
    io.iris_save(cube, netcdf_path)
    write_provenance(netcdf_path, provenance_record, cfg)


def plot_multiple_datasets(cfg, datasets, short_name):
    """Create time series plot for multiple datasets.

    Parameters
    ----------
    cfg: dict
        Global configuration dictionary.
    datasets: list of dict
        List of metadata dictionaries for all datasets.
    short_name: str
        Short name of the variable to plot.

    """
    if not datasets:
        return
    logger.info("Plotting time series for multiple datasets for variable '%s'",
                short_name)
    cubes = {}
    for dataset in datasets:
        cube = load_cube(dataset, ['time'])
        if cube is None:
            continue
        global_time_series_plot(cube, ls='-', lw=2.0, label=dataset['dataset'])
        cubes[dataset['dataset']] = cube
    if not cubes:
        return

    # Provenance
    caption = f"Time series plot of global mean {short_name}."
    provenance_record = get_provenance_record(
        caption, [d['filename'] for d in datasets])

    # Plot appearance
    set_plot_title(datasets[0]['long_name'])
    legend = plt.legend(loc='center left',
                        bbox_to_anchor=[1.05, 0.5],
                        borderaxespad=0.0)
    plt.xlabel('year')
    plt.ylabel(f"{short_name} / {cube.units}")
    if 'time_range' in cfg:
        plt.xlim(cfg['time_range'][0], cfg['time_range'][1])
    if 'y_range' in cfg:
        plt.ylim(cfg['y_range'][0], cfg['y_range'][1])

    # Save plot
    basename = '_'.join([short_name, 'global_timeseries'])
    plot_provenance = save_plot(basename,
                                'times',
                                cfg,
                                additional_artists=[legend])
    provenance_record.update(plot_provenance)

    # Write netcdf file and provenance
    netcdf_path = get_diagnostic_filename(basename, cfg)
    var_attrs = {
        'short_name': short_name,
        'long_name': datasets[0]['long_name'],
        'units': datasets[0]['units'],
    }
    if datasets[0].get('standard_name'):
        var_attrs['standard_name'] = datasets[0]['standard_name']
    io.save_1d_data(cubes, netcdf_path, 'time', var_attrs)
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

            # Create global time series for every dataset
            plot_single_dataset(cfg, dataset)

        if cfg.get('multi_dataset_plot'):
            plot_multiple_datasets(cfg, datasets, short_name)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
