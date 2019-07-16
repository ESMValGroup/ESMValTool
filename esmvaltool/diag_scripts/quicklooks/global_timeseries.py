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
time_range : list of float, optional:
    Range for the time axis in the plots.
y_range : list of float, optional:
    Range for the variable axis.

"""

import logging
import os

import iris
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            get_plot_filename, group_metadata,
                                            io, run_diagnostic)

logger = logging.getLogger(os.path.basename(__file__))


def get_provenance_record(caption, ancestors):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_type': 'times',
        'authors': ['bock_ls'],
        'references': ['acknow_project'],
        'ancestors': ancestors,
    }
    return record


def load_cube(dataset):
    """Load cube and compute global average if necessary."""
    filename = dataset['filename']
    logger.debug("Loading '%s'", filename)
    cube = iris.load_cube(filename)
    coords = [coord.name() for coord in cube.coords(dim_coords=True)]

    # Check if cubes has desired coordinates
    if 'time' not in coords:
        raise iris.exceptions.CoordinateNotFoundError(
            f"File '{filename}' does not contain necessary coordinate 'time'")
    coords.remove('time')

    # Calculate global mean
    if cube.ndim > 1:
        if 'latitude' in coords:
            grid_areas = iris.analysis.cartography.area_weights(cube)
        else:
            grid_areas = None
        cube = cube.collapsed(coords, iris.analysis.MEAN, weights=grid_areas)

    return cube


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
    cubedata = np.ma.array(cube.data)
    if len(cubedata.compressed()) == 1:
        plt.axhline(cubedata.compressed(), **kwargs)
        return
    times = diagtools.cube_time_to_float(cube)
    plt.plot(times, cubedata, **kwargs)


def plot_single_dataset(cfg, dataset):
    """Plot global mean time series for single dataset.

    Parameters
    ----------
    cfg: dict
        Global configuration dictionary.
    dataset: list
        Metadata dictionary for desired dataset.

    """
    cube = load_cube(dataset)

    # Provenance
    provenance_record = get_provenance_record(
        f"Time series plot of global mean {dataset['short_name']} for dataset "
        f"{dataset['dataset']}.", [dataset['filename']])

    # Plot
    global_time_series_plot(cube)
    title = ' '.join([dataset['dataset'], dataset['long_name']])
    plt.title(title)
    plt.xlabel('year')
    plt.ylabel(f"{dataset['short_name']} / {cube.units}")
    if 'time_range' in cfg:
        plt.xlim(cfg['time_range'][0], cfg['time_range'][1])
    if 'y_range' in cfg:
        plt.ylim(cfg['y_range'][0], cfg['y_range'][1])

    # Save plot if desired
    if cfg['write_plots']:
        plot_path = get_plot_filename(
            '_'.join([
                dataset['short_name'],
                dataset['dataset'],
                'global_timeseries',
            ]), cfg)
        plt.savefig(plot_path, bbox_inches='tight', orientation='landscape')
        logger.info("Wrote %s", plot_path)
        provenance_record['plot_file'] = plot_path
    plt.close()

    # Write netcdf file
    netcdf_path = get_diagnostic_filename(
        '_'.join([
            dataset['short_name'],
            dataset['dataset'],
            'global_timeseries',
        ]), cfg)
    io.iris_save(cube, netcdf_path)

    # Write provenance
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)


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
        cube = load_cube(dataset)
        global_time_series_plot(cube, ls='-', lw=2.0, label=dataset['dataset'])
        cubes[dataset['dataset']] = cube

    # Provenance
    provenance_record = get_provenance_record(
        f"Time series plot of global mean {short_name}.",
        [d['filename'] for d in datasets])

    # Plot appearance
    plt.title(datasets[0]['long_name'])
    legend = plt.legend(loc='center left',
                        bbox_to_anchor=[1.05, 0.5],
                        borderaxespad=0.0)
    plt.xlabel('year')
    plt.ylabel(f"{short_name} / {cube.units}")
    if 'time_range' in cfg:
        plt.xlim(cfg['time_range'][0], cfg['time_range'][1])
    if 'y_range' in cfg:
        plt.ylim(cfg['y_range'][0], cfg['y_range'][1])

    # Save plot if desired
    if cfg['write_plots']:
        plot_path = get_plot_filename(
            '_'.join([short_name, 'global_timeseries']), cfg)
        plt.savefig(plot_path,
                    bbox_inches='tight',
                    orientation='landscape',
                    additional_artists=[legend])
        logger.info("Wrote %s", plot_path)
        provenance_record['plot_file'] = plot_path
    plt.close()

    # Write netcdf file
    netcdf_path = get_diagnostic_filename(
        '_'.join([short_name, 'global_timeseries']), cfg)
    var_attrs = {
        'short_name': short_name,
        'long_name': datasets[0]['long_name'],
        'units': datasets[0]['units'],
    }
    if datasets[0].get('standard_name'):
        var_attrs['standard_name'] = datasets[0]['standard_name']
    io.save_1d_data(cubes, netcdf_path, 'time', var_attrs)

    # Write provenance
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)


def main(cfg):
    """Load config file and metadata, then pass them the plot making tools.

    Parameters
    ----------
    cfg: dict
        Global configuration dictionary.

    """
    if cfg['quicklook']['active']:
        quicklook_dir = cfg['quicklook']['output_dir']
        logger.info("Reading data from quicklook directory %s", quicklook_dir)
        input_data = io.netcdf_to_metadata(cfg, root=quicklook_dir)
    else:
        logger.info("Reading data regular ESMValTool directory")
        input_data = cfg['input_data'].values()

    # Group data in terms of variables
    grouped_data = group_metadata(input_data, 'short_name')

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
