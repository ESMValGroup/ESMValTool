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

import iris
import matplotlib.pyplot as plt

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
        'plot_type': 'zonal',
        'authors': ['bock_ls'],
        'references': ['acknow_project'],
        'ancestors': ancestors,
    }
    return record


def load_cube(dataset):
    """Load cube and compute zonal average if necessary."""
    filename = dataset['filename']
    logger.debug("Loading '%s'", filename)
    cube = iris.load_cube(filename)
    coords = [coord.name() for coord in cube.coords(dim_coords=True)]

    # Check if cubes has desired coordinates
    for coord_name in ('time', 'latitude'):
        if coord_name not in coords:
            logger.warning(
                "File '%s' does not contain necessary coordinate '%s', "
                "skipping", filename, coord_name)
            return None
        coords.remove(coord_name)

    # Calculate zonal mean
    if cube.ndim > 2:
        cube = cube.collapsed(coords, iris.analysis.MEAN)

    return cube


def plot_single_dataset(cfg, dataset):
    """Plot zonal mean time series for single dataset.

    Parameters
    ----------
    cfg: dict
        Global configuration dictionary.
    dataset: list
        Metadata dictionary for desired dataset.

    """
    cube = load_cube(dataset)
    if cube is None:
        return

    # Provenance
    provenance_record = get_provenance_record(
        f"Time series plot of zonal mean {dataset['short_name']} for dataset "
        f"{dataset['dataset']}.", [dataset['filename']])

    # Transpose cube
    data = cube.data.transpose()
    times = diagtools.cube_time_to_float(cube)

    # Plot
    plt.contourf(times,
                 cube.coord('latitude').points,
                 data,
                 extend='both',
                 levels=cfg.get('levels'))
    title = ' '.join([dataset['dataset'], dataset['long_name']])
    plt.title(title)
    plt.xlabel('year')
    plt.ylabel('latitude')
    if 'time_range' in cfg:
        plt.xlim(cfg['time_range'][0], cfg['time_range'][1])
    if 'latitude_range' in cfg:
        plt.ylim(cfg['latitude_range'][0], cfg['latitude_range'][1])

    # Make a colorbar for the ContourSet returned by the contourf call.
    colorbar = plt.colorbar(orientation='horizontal', aspect=30)
    colorbar.set_label(f"{dataset['short_name']} / {cube.units}")

    # Saving files:
    if cfg['write_plots']:
        plot_path = get_plot_filename(
            '_'.join([
                dataset['short_name'],
                dataset['dataset'],
                'zonal_timeseries',
            ]), cfg)
        plt.savefig(plot_path, bbox_inches='tight', orientation='landscape')
        logger.info("Wrote %s", plot_path)
        provenance_record['plot_file'] = plot_path
    plt.close()

    # Write netcdf file for every plot
    netcdf_path = get_diagnostic_filename(
        '_'.join([
            dataset['short_name'],
            dataset['dataset'],
            'zonal_timeseries',
        ]), cfg)
    io.iris_save(cube, netcdf_path)

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

            # Create zonal time series for every dataset
            plot_single_dataset(cfg, dataset)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
