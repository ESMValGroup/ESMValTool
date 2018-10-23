"""
Diagnostic transect:

Diagnostic to produce images of a transect.
These plost show either latitude or longitude against depth, and the cube value
is used as the colour scale.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has no time component, and one of the latitude or
longitude coordinates has been reduced to a single value.

An approproate preprocessor for a 3D+time field would be:
preprocessors:
  prep_transect:
    time_average:
    extract_slice: # Atlantic Meridional Transect
      latitude: [-50.,50.]
      longitude: 332.

This tool is part of the ocean diagnostic tools package in the ESMValTool.

Author: Lee de Mora (PML)
        ledm@pml.ac.uk
"""
import logging
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')  # noqa
import iris

import iris.quickplot as qplt
import matplotlib.pyplot as plt

import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def determine_transect_str(cube):
    """
    Determine the Transect String

    Takes a guess at a string to describe the transect.
    """
    options = ['latitude', 'longitude']
    for option in options:
        coord = cube.coord(option)
        if len(coord.points) > 1:
            continue
        value = coord.points.mean()
        value = round(value, 2)
        if option == 'latitude':
            return str(value) + ' N'
        if option == 'longitude':
            if value > 180.:
                return str(value - 360.) + ' W'
            return str(value) + ' E'
    return ''


def make_transects_plots(
        cfg,
        metadata,
        filename,
):
    """
    Make a simple plot of the transect for an indivudual model.

    The cfg is the opened global config,
    metadata is the metadata dictionairy
    filename is the preprocessing model file.
    """
    # Load cube and set up units
    cube = iris.load_cube(filename)
    cube = diagtools.bgc_units(cube, metadata['short_name'])

    # Is this data is a multi-model dataset?
    multi_model = metadata['dataset'].find('MultiModel') > -1

    # Make a dict of cubes for each layer.

    qplt.contourf(cube, 25, linewidth=0, rasterized=True)
    plt.axes().set_yscale('log')

    # Add title to plot
    title = ' '.join(
        [metadata['dataset'], metadata['long_name'],
         determine_transect_str(cube)])
    plt.title(title)

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    if multi_model:
        path = diagtools.folder(
            cfg['plot_dir']) + os.path.basename(filename).replace(
                '.nc', '_transect' + image_extention)
    else:
        path = diagtools.get_image_path(
            cfg,
            metadata,
            suffix='transect' + image_extention,
        )

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)

    plt.close()


def add_sea_floor(cube):
    """Add a simple sea floor line from the cube mask."""
    land_cube = cube.copy()
    mask = 1. * land_cube.data.mask
    # mask = np.ma.masked_where(mask==0, mask)
    land_cube.data = np.ma.masked_where(mask == 0, mask)
    land_cube.data.mask = np.zeros_like(mask)
    qplt.contour(land_cube, 2,
                 cmap='Greys_r',
                 rasterized=True)


def make_transect_contours(
        cfg,
        metadata,
        filename,
):
    """
    Make a simple contour plot of the transect for an indivudual model.

    The cfg is the opened global config,
    metadata is the metadata dictionairy
    filename is the preprocessing model file.
    """
    # Load cube and set up units
    cube = iris.load_cube(filename)
    cube = diagtools.bgc_units(cube, metadata['short_name'])

    # Is this data is a multi-model dataset?
    multi_model = metadata['dataset'].find('MultiModel') > -1

    # Load threshold/thresholds.
    plot_details = {}
    colours = []
    if 'threshold' in cfg.keys():
        thresholds = [float(cfg['threshold']), ]
    elif 'thresholds' in cfg.keys():
        thresholds = [float(thres) for thres in cfg['thresholds']]
    linewidths = [1 for thres in thresholds]
    linestyles = ['-' for thres in thresholds]
    for itr, thres in enumerate(thresholds):
        if len(thresholds) > 1:
            colour = plt.cm.jet(float(itr) / float(len(thresholds) - 1.))
        else:
            colour = plt.cm.jet(0)
        label = str(thres) + ' ' + str(cube.units)
        colours.append(colour)
        plot_details[thres] = {'c': colour, 'lw': 1, 'ls': '-', 'label': label}

    qplt.contour(cube,
                 thresholds,
                 colors=colours,
                 linewidths=linewidths,
                 linestyles=linestyles,
                 rasterized=True)
    plt.axes().set_yscale('log')

    add_sea_floor(cube)

    # Add legend
    diagtools.add_legend_outside_right(plot_details,
                                       plt.gca(),
                                       column_width=0.08,
                                       loc='below')

    # Add title to plot
    title = ' '.join(
        [metadata['dataset'], metadata['long_name'],
         determine_transect_str(cube)])
    plt.title(title)

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    if multi_model:
        path = diagtools.folder(
            cfg['plot_dir']) + os.path.basename(filename)
        path.replace('.nc', '_transect_contour' + image_extention)
    else:
        path = diagtools.get_image_path(
            cfg,
            metadata,
            suffix='transect_contour' + image_extention,
        )

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)

    plt.close()


def multi_model_contours(
        cfg,
        metadata,
):
    """
    Make a contour transect showing several models.

    The cfg is the opened global config,
    metadata is the metadata dictionairy.
    """
    ####
    # Load the data for each layer as a separate cube
    model_cubes = {}
    for filename in sorted(metadata):
        cube = iris.load_cube(filename)
        model_cubes[filename] = cube

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)
    # Load threshold/thresholds.
    if 'threshold' in cfg.keys():
        thresholds = [float(cfg['threshold']), ]
    elif 'thresholds' in cfg.keys():
        thresholds = [float(thres) for thres in cfg['thresholds']]

    # Make a plot for each layer and each threshold
    for threshold in thresholds:
        print('\n', thresholds, threshold)
        title = ''
        plot_details = {}
        cmap = plt.cm.get_cmap('jet')

        # Plot each file in the group
        for index, filename in enumerate(sorted(metadata)):
            if len(metadata) > 1:
                color = cmap(index / (len(metadata) - 1.))
            else:
                color = 'blue'
            linewidth = 1.
            linestyle = '-'
            # Determine line style for MultiModel statistics:
            if 'MultiModel' in metadata[filename]['dataset']:
                linewidth = 2.
                linestyle = ':'
            # Determine line style for Observations
            if metadata[filename]['project'] in diagtools.get_obs_projects():
                color = 'black'
                linewidth = 1.7
                linestyle = '-'

            qplt.contour(model_cubes[filename],
                         [threshold, ],
                         colors=[color, ],
                         linewidths=linewidth,
                         linestyles=linestyle,
                         rasterized=True)

            plot_details[filename] = {'c': color,
                                      'ls': linestyle,
                                      'lw': linewidth,
                                      'label': metadata[filename]['dataset']}

            plt.axes().set_yscale('log')

            title = metadata[filename]['long_name']
            units = str(model_cubes[filename].units)

            add_sea_floor(model_cubes[filename])

        # Add title, threshold, legend to plots
        title = ' '.join([title,
                          str(threshold),
                          units,
                          determine_transect_str(model_cubes[filename])])
        plt.title(title)
        plt.legend(loc='best')

        # Saving files:
        if cfg['write_plots']:
            path = diagtools.get_image_path(
                cfg,
                metadata[filename],
                prefix='MultipleModels',
                suffix='_'.join(['contour_tramsect',
                                 str(threshold) + image_extention]),
                metadata_id_list=[
                    'field', 'short_name', 'preprocessor', 'diagnostic',
                    'start_year', 'end_year'
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
    Load the config file, and send it to the plot maker.

    The cfg is the opened global config.
    """
    #####
    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info(
            'metadata filename:\t%s',
            metadata_filename,
        )

        metadata = diagtools.get_input_files(cfg, index=index)

        #######
        # Multi model contour plots
        multi_model_contours(
            cfg,
            metadata,
        )

        for filename in sorted(metadata.keys()):

            logger.info('-----------------')
            logger.info(
                'model filenames:\t%s',
                filename,
            )

            ######
            # Time series of individual model
            make_transects_plots(cfg, metadata[filename], filename)

            ######
            # Contour maps of individual model
            if 'threshold' in cfg.keys() or 'thresholds' in cfg.keys():
                make_transect_contours(cfg, metadata[filename], filename)

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
