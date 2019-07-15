"""
Zonal time series diagnostics
=======================

Description
-----------
Plot time series of zonal mean
In quicklook mode (cfg['quicklook']['active']: True) this diagnostic plots
the concatinated file.

Author
------
Lisa Bock (DLR, Germany)

Project
-------
CMIP6-DICAD

Configuration options in recipe
-------------------------------
val_levs: values for contour levels
lat_int: min and max latitude values
time_int: min and max for time axis

"""


import logging
import os

import iris
import matplotlib.pyplot as plt

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools

from esmvaltool.diag_scripts.shared._base import (
    ProvenanceLogger, get_diagnostic_filename,
    run_diagnostic)

logger = logging.getLogger(os.path.basename(__file__))


def get_provenance_record(caption):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_type': 'zonal',
        'authors': ['bock_ls'],
        'references': ['acknow_project'],
    }
    return record


def make_zon_time_series_plots(
        cfg,
        metadata,
        filename,
):
    """
    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.
    filename: str
        The preprocessed model file

    """

    # Load cube
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    data = cube.data.transpose()

    #y_axis = 'latitude'
    times = diagtools.cube_time_to_float(cube)

    # set the levels for contour
    if 'val_levs' in cfg:
        levels = cfg['val_levs']
    else:
        levels = None

    fig = plt.contourf(times,
                       cube.coord('latitude').points,
                       data, extend='both', levels=levels)

    # Add title, legend to plots
    title = ' '.join([metadata['dataset'], metadata['long_name']])
    plt.title(title)
    plt.ylabel('lat')

    # set the limits for time axis
    if 'time_int' in cfg:
        plt.gca().set_xlim(cfg['time_int'])
    # set the lat-limits
    if 'lat_int' in cfg:
        plt.gca().set_ylim(cfg['lat_int'])

    # Make a colorbar for the ContourSet returned by the contourf call.
    cbar = plt.colorbar(fig)
    cbar.ax.set_ylabel(str(cube.units))

    # Determine image filename:
    path = diagtools.get_image_path(
        cfg,
        metadata,
        prefix='Model',
        suffix='zonal_timeseries' + image_extention,
        metadata_id_list=[
            'dataset', 'field', 'short_name'
        ],
    )

    # Saving files:
    if cfg['write_plots']:

        logger.info('Saving plots to %s', path)
        plt.savefig(path)

    plt.close()

    # Write netcdf file for every plot
    diagname = '_'.join([metadata['dataset'], metadata['short_name'],
                         'zonal_timeseries'])
    diagnostic_file = get_diagnostic_filename(diagname, cfg)
    logger.info("Saving analysis results to %s", diagnostic_file)
    iris.save(cube, target=diagnostic_file)

    # Provenance
    provenance_record = get_provenance_record(
        "Timeseries of zonal mean of {} for dataset {}."
        .format(metadata['short_name'], metadata['dataset']))
    provenance_record.update({
        'plot_file': path,
    })

    return (diagnostic_file, provenance_record)


def main(cfg):
    """
    Load the config file and some metadata, then pass them the plot making
    tools.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.

    """

    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info('metadata filename:\t%s', metadata_filename)

        metadatas = diagtools.get_input_files(cfg, index=index)

        for filename in sorted(metadatas):

            logger.info('-----------------')
            logger.info(
                'model filenames:\t%s',
                filename,
            )

            metadata = metadatas[filename]

            if cfg['quicklook']['active']:
                # if quicklook mode - plotting of concatinating file
                # path to concatinated file
                quicklook_dir = cfg['quicklook']['output_dir']
                filename = quicklook_dir + '/'
                filename += '_'.join([metadata['dataset'],
                                      metadata['short_name'] + '.nc'])
                logger.info('concatinated filename:\t%s', filename)

            # Time series of individual model
            (path, provenance_record) = make_zon_time_series_plots(
                cfg, metadata, filename)

            # Provenance
            if path is not None:
                provenance_record['ancestors'] = filename
                with ProvenanceLogger(cfg) as provenance_logger:
                    provenance_logger.log(path, provenance_record)

    logger.info('Success')


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
