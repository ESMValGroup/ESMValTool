"""
Diagnostic Maps:

Diagnostic to produce images of a map with coastlines from a cube.
These plost show latitude vs longitude and the cube value is used as the colour
scale.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has no time component, a small number of depth layers,
and a latitude and longitude coordinates.

An approproate preprocessor for a 3D+time field would be:
preprocessors:
  prep_map:
    extract_levels:
      levels:  [100., ]
      scheme: linear_extrap
    time_average:

This tool is part of the ocean diagnostic tools package in the ESMValTool.

Author: Lee de Mora (PML)
        ledm@pml.ac.uk
"""
import logging
import os
import sys
from itertools import product
import matplotlib
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt

import iris
import iris.quickplot as qplt
import cartopy

from esmvaltool.preprocessor._regrid import _stock_cube

import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def regrid_irregulars(cube,scheme='nearest'):
    """Regrid irregular grids."""
    lats = cube.coord('latitude')
    if lats.ndim == 1: return cube
    print(cube)
    horizontal_schemes = dict(
        linear=iris.analysis.Linear(extrapolation_mode='mask'),
        nearest=iris.analysis.Nearest(extrapolation_mode='mask'),
        area_weighted=iris.analysis.AreaWeighted(),
        unstructured_nearest=iris.analysis.UnstructuredNearest())
    target_grid = _stock_cube('1x1')
    print(target)
    return cube.regrid(target_grid, horizontal_schemes[scheme])


def make_map_plots(
        cfg,
        metadata,
        filename,
):
    """
    Make a simple map plot for an individual model.

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
    cubes = diagtools.make_cube_layer_dict(cube)

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Making plots for each layer
    for layer_index, (layer, cube_layer) in enumerate(cubes.items()):
        layer = str(layer)

        # cube_layer = regrid_irregulars(cube_layer)
        qplt.contourf(cube_layer, 25, linewidth=0, rasterized=True)

        try:
            plt.gca().coastlines()
        except AttributeError:
            logger.warning('Not able to add coastlines')
        try:
                plt.gca().add_feature(cartopy.feature.LAND,
                              zorder=10,
                              facecolor=[0.8, 0.8, 0.8])
        except AttributeError:
           logger.warning('Not able to add coastlines')

        # Add title to plot
        title = ' '.join([metadata['dataset'], metadata['long_name']])
        if layer:
            title = ' '.join(
                [title, '(', layer,
                 str(cube_layer.coords('depth')[0].units), ')'])
        plt.title(title)

        # Determine image filename:
        if multi_model:
            path = diagtools.folder(
                cfg['plot_dir']) + os.path.basename(filename).replace(
                    '.nc', '_map_' + str(layer_index) + image_extention)
        else:
            path = diagtools.get_image_path(
                cfg,
                metadata,
                suffix='map_' + str(layer_index) + image_extention,
            )

        # Saving files:
        if cfg['write_plots']:

            logger.info('Saving plots to %s', path)
            plt.savefig(path)

        plt.close()


def make_map_contour(
        cfg,
        metadata,
        filename,
):
    """
    Make a simple contour map plot for an individual model.

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
    cubes = diagtools.make_cube_layer_dict(cube)

    # Load image format extention and threshold.thresholds.
    image_extention = diagtools.get_image_format(cfg)

    # Load threshold/thresholds.
    plot_details = {}
    colours = []
    if 'threshold' in cfg.keys():
        thresholds = [float(cfg['threshold']), ]
    elif 'thresholds' in cfg.keys():
        thresholds = [float(thres) for thres in cfg['thresholds']]

    for itr, thres in enumerate(thresholds):
        if len(thresholds) > 1:
            colour = plt.cm.jet(float(itr) / float(len(thresholds) - 1.))
        else:
            colour = plt.cm.jet(0)
        label = str(thres) + ' ' + str(cube.units)
        colours.append(colour)
        plot_details[thres] = {'c': colour,
                               'lw': 1,
                               'ls': '-',
                               'label': label}

    linewidths = [1 for thres in thresholds]
    linestyles = ['-' for thres in thresholds]
    # Making plots for each layer
    for layer_index, (layer, cube_layer) in enumerate(cubes.items()):
        layer = str(layer)
        # cube_layer = regrid_irregulars(cube_layer)
        qplt.contour(cube_layer,
                     thresholds,
                     colors=colours,
                     linewidths=linewidths,
                     linestyles=linestyles,
                     rasterized=True)

        try:
            plt.gca().coastlines()
        except AttributeError:
            logger.warning('Not able to add coastlines')
        try:
            plt.gca().add_feature(cartopy.feature.LAND,
                              zorder=10,
                              facecolor=[0.8, 0.8, 0.8])
        except AttributeError:
            logger.warning('Not able to add coastlines')
        # Add legend
        diagtools.add_legend_outside_right(plot_details,
                                           plt.gca(),
                                           column_width=0.02,
                                           loc='below')

        # Add title to plot
        title = ' '.join([metadata['dataset'], metadata['long_name']])
        if layer:
            title = ' '.join(
                [title, '(', layer,
                 str(cube_layer.coords('depth')[0].units), ')'])
        plt.title(title)

        # Determine image filename:
        if multi_model:
            path = diagtools.folder(cfg['plot_dir'])
            path = path + os.path.basename(filename)
            path = path.replace('.nc', '_contour_map_' + str(layer_index))
            path = path + image_extention
        else:
            path = diagtools.get_image_path(
                cfg,
                metadata,
                suffix='_contour_map_' + str(layer_index) + image_extention,
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
    Make a contour map showing several models.

    The cfg is the opened global config,
    metadata is the metadata dictionairy.
    """
    ####
    # Load the data for each layer as a separate cube
    model_cubes = {}
    layers = {}
    for filename in sorted(metadata):
        cube = iris.load_cube(filename)

        cubes = diagtools.make_cube_layer_dict(cube)
        model_cubes[filename] = cubes
        for layer in cubes:
            layers[layer] = True

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Load threshold/thresholds.
    if 'threshold' in cfg.keys():
        thresholds = [float(cfg['threshold']), ]
    elif 'thresholds' in cfg.keys():
        thresholds = [float(thres) for thres in cfg['thresholds']]

    # Make a plot for each layer and each threshold
    for layer, threshold in product(layers, thresholds):

        title = ''
        z_units = ''
        plot_details = {}
        cmap = plt.cm.get_cmap('jet')
        land_drawn = False

        # Plot each file in the group
        for index, filename in enumerate(sorted(metadata)):
            if len(metadata) > 1:
                color = cmap(index / (len(metadata) - 1.))
            else:
                color = 'blue'

            cube =   model_cubes[filename][layer]
        #     cube = regrid_irregulars(cube)

            print(index, os.path.basename(filename), layer, threshold, model_cubes[filename][layer].shape)
            if 'MultiModel' in metadata[filename]['dataset']:
                qplt.contour(cube,
                             [threshold, ],
                             colors=[color, ],
                             linewidths=2.,
                             linestyles=':',
                             rasterized=True)
                plot_details[filename] = {
                    'c': color,
                    'ls': ':',
                    'lw': 2.,
                    'label': metadata[filename]['dataset']
                }
            else:
                qplt.contour(cube,
                             [threshold, ],
                             colors=[color, ],
                             linewidths=1.,
                             linestyles='-',
                             rasterized=True)
                plot_details[filename] = {
                    'c': color,
                    'ls': '-',
                    'lw': 2.,
                    'label': metadata[filename]['dataset']
                }

            if not land_drawn:
                try:
                    plt.gca().coastlines()
                except AttributeError:
                    logger.warning('Not able to add coastlines')
                plt.gca().add_feature(cartopy.feature.LAND,
                                      zorder=10,
                                      facecolor=[0.8, 0.8, 0.8])
                land_drawn = True

            title = metadata[filename]['long_name']
            if layer != '':
                z_units = model_cubes[filename][layer].coords('depth')[0].units
            units = str(model_cubes[filename][layer].units)

        # Add title, threshold, legend to plots
        title = ' '.join([title, str(threshold), units])
        if layer:
            title = ' '.join([title,  '(', str(layer), str(z_units), ')'])
        plt.title(title)
        plt.legend(loc='best')

        # Saving files:
        if cfg['write_plots']:
            path = diagtools.get_image_path(
                cfg,
                metadata[filename],
                prefix='MultipleModels_',
                suffix='_'.join(['_contour_map_',
                                 str(threshold),
                                 str(layer) + image_extention]),
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
    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info(
            'metadata filename:\t%s',
            metadata_filename,
        )

        metadatas = diagtools.get_input_files(cfg, index=index)

        #######
        # Multi model contour plots
        multi_model_contours(
            cfg,
            metadatas,
        )

        for filename in sorted(metadatas.keys()):

            logger.info('-----------------')
            logger.info(
                'model filenames:\t%s',
                filename,
            )

            ######
            # Contour maps of individual model
            if 'threshold' in cfg.keys() or 'thresholds' in cfg.keys():
                make_map_contour(cfg, metadatas[filename], filename)

            ######
            # Mps of individual model
            make_map_plots(cfg, metadatas[filename], filename)

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
