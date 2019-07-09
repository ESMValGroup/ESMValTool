"""
Time series diagnostics
=======================
"""

import logging
import os
from pprint import pformat

import iris
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

logger = logging.getLogger(os.path.basename(__file__))


def compute_diagnostic(filename):
    """Compute global average."""
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)

    grid_areas = iris.analysis.cartography.area_weights(cube)

    # from zonal mean to global mean
    return cube.collapsed('latitude', iris.analysis.MEAN, weights=grid_areas)


def timeplot(cube, **kwargs):
    """
    Create a time series plot from the cube.

    Note that this function simple does the plotting, it does not save the
    image or do any of the complex work. This function also takes and of the
    key word arguments accepted by the matplotlib.pyplot.plot function.
    These arguments are typically, color, linewidth, linestyle, etc...

    If there's only one datapoint in the cube, it is plotted as a
    horizontal line.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube

    """
    cubedata = np.ma.array(cube.data)
    if len(cubedata.compressed()) == 1:
        plt.axhline(cubedata.compressed(), **kwargs)
        return

    times = diagtools.cube_time_to_float(cube)
    plt.plot(times, cubedata, **kwargs)


def make_time_series_plots(
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

    # Load and compute cube
    cube = compute_diagnostic(filename)

    # Set up units
    cube = diagtools.bgc_units(cube, metadata['short_name'])

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    timeplot(cube, label=metadata['dataset'])

    # Add title, legend to plots
    title = ' '.join([metadata['dataset'], metadata['long_name']])
    plt.title(title)
    #plt.legend(loc='best')
    plt.ylabel(str(cube.units))
 
    # set the y-limits
    if 'y_min' in cfg:
        plt.ylim(bottom=cfg['y_min'])
    if 'y_max' in cfg:
        plt.ylim(top=cfg['y_max'])

    # Determine image filename:
    path = diagtools.get_image_path(
        cfg,
        metadata,
        prefix='Model',
        suffix='global_timeseries' + image_extention,
        metadata_id_list=[
            'dataset', 'field', 'short_name', 'start_year', 'end_year'
        ],
    )

    # Saving files:
    if cfg['write_plots']:

        logger.info('Saving plots to %s', path)
        plt.savefig(path)

    plt.close()


def multi_model_time_series(
        cfg,
        metadata,
):
    """
    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.

    """

    model_cube = {}
    for filename in sorted(metadata):
        # Load and compute cube
        cube = compute_diagnostic(filename)

        # Set up units
        cube = diagtools.bgc_units(cube, metadata[filename]['short_name'])

        model_cube[filename] = cube

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Make the plot
    title = ''
    plot_details = {}
    cmap = plt.cm.get_cmap('viridis')

    # Plot each file in the group
    for index, filename in enumerate(sorted(metadata)):
        if len(metadata) > 1:
            color = cmap(index / (len(metadata) - 1.))
        else:
            color = 'blue'

        cube = model_cube[filename]

        timeplot(
            cube,
            c=color,
            # label=metadata[filename]['dataset'])
            ls='-',
            lw=2.,
        )
        plot_details[filename] = {
            'c': color,
            'ls': '-',
            'lw': 2.,
            'label': metadata[filename]['dataset']
        }

        title = metadata[filename]['long_name']
    # Add title, legend to plots
    plt.title(title)
    plt.legend(loc='best')
    plt.ylabel(str(model_cube[filename].units))

    # set the y-limits
    if 'y_min' in cfg:
        plt.ylim(bottom=cfg['y_min'])
    if 'y_max' in cfg:
        plt.ylim(top=cfg['y_max'])

    # Saving files:
    if cfg['write_plots']:
        path = diagtools.get_image_path(
            cfg,
            metadata[filename],
            prefix='MultiModel',
            suffix='global_timeseries' + image_extention,
            metadata_id_list=[
                'field', 'short_name', 'start_year', 'end_year'
            ],
        )

    # Resize and add legend outside thew axes.
    plt.gcf().set_size_inches(9., 6.)
    diagtools.add_legend_outside_right(
        plot_details, plt.gca(), column_width=0.15)

    logger.info('Saving plots to %s', path)
    plt.savefig(path)
    plt.close()



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

        if len(metadatas) > 1:
            # Time series plot with all models
            multi_model_time_series(
                cfg,
                metadatas,
            )

        for filename in sorted(metadatas):

            logger.info('-----------------')
            logger.info(
                'model filenames:\t%s',
                filename,
            )

            # Time series of individual model
            make_time_series_plots(cfg, metadatas[filename], filename)
    logger.info('Success')


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
