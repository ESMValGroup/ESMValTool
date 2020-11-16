"""
Maps diagnostics
================

Diagnostic to produce images of a map with coastlines from a cube.
These plost show latitude vs longitude and the cube value is used as the colour
scale.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has no time component, a small number of depth layers,
and a latitude and longitude coordinates.

An approproate preprocessor for a 3D+time field would be::

  preprocessors:
    prep_map:
      extract_levels:
        levels:  [100., ]
         scheme: linear_extrap
      climate_statistics:
        operator: mean


Note that this recipe may not function on machines with no access to the
internet, as cartopy may try to download the shapefiles. The solution to
this issue is the put the relevant cartopy shapefiles on a disk visible to your
machine, then link that path to ESMValTool via the `auxiliary_data_dir`
variable. The cartopy masking files can be downloaded from::

  https://www.naturalearthdata.com/downloads/

Here, cartopy uses the 1:10, physical coastlines and land files::

      110m_coastline.dbf  110m_coastline.shp  110m_coastline.shx
      110m_land.dbf  110m_land.shp  110m_land.shx

This tool is part of the ocean diagnostic tools package in the ESMValTool.

Author: Lee de Mora (PML)
        ledm@pml.ac.uk
"""
import itertools
import logging
import os
import sys

import cartopy
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def make_map_plots(
        cfg,
        metadata,
        filename,
):
    """
    Make a simple map plot for an individual model.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.
    metadata: dict
        the metadata dictionary
    filename: str
        the preprocessed model file.

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

        qplt.contourf(cube_layer, 25, linewidth=0, rasterized=True)

        try:
            plt.gca().coastlines()
        except AttributeError:
            logger.warning('Not able to add coastlines')

        # Add title to plot
        title = ' '.join([metadata['dataset'], metadata['long_name']])
        if layer:
            title = ' '.join([
                title, '(', layer,
                str(cube_layer.coords('depth')[0].units), ')'
            ])
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

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.
    metadata: dict
        the metadata dictionary
    filename: str
        the preprocessed model file.

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
    thresholds = diagtools.load_thresholds(cfg, metadata)

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
        depth_units = str(cube_layer.coords('depth')[0].units)
        if layer:
            title = '{} ({} {})'.format(title, layer, depth_units)
        plt.title(title)

        # Determine image filename:
        if multi_model:
            path = os.path.join(diagtools.folder(cfg['plot_dir']),
                                os.path.basename(filename))
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

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.
    metadata: dict
        the metadata dictionary.

    """
    ####
    # Load the data for each layer as a separate cube
    model_cubes = {}
    layers = {}
    for filename in sorted(metadata):
        cube = iris.load_cube(filename)
        cube = diagtools.bgc_units(cube, metadata[filename]['short_name'])

        cubes = diagtools.make_cube_layer_dict(cube)
        model_cubes[filename] = cubes
        for layer in cubes:
            layers[layer] = True

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Load threshold/thresholds.
    thresholds = diagtools.load_thresholds(cfg, metadata)

    # Make a plot for each layer and each threshold
    for layer, threshold in itertools.product(layers, thresholds):

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
            linewidth = 1.
            linestyle = '-'

            # Determine line style for Observations
            if metadata[filename]['project'] in diagtools.get_obs_projects():
                color = 'black'
                linewidth = 1.7
                linestyle = '-'

            # Determine line style for MultiModel statistics:
            if 'MultiModel' in metadata[filename]['dataset']:
                color = 'black'
                linestyle = ':'
                linewidth = 1.4

            cube = model_cubes[filename][layer]
            qplt.contour(cube,
                         [threshold, ],
                         colors=[color, ],
                         linewidths=linewidth,
                         linestyles=linestyle,
                         rasterized=True)
            plot_details[filename] = {
                'c': color,
                'ls': linestyle,
                'lw': linewidth,
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
            title = ' '.join([title, '(', str(layer), str(z_units), ')'])
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
    Load the config file, and send it to the plot makers.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.

    """
    cartopy.config['data_dir'] = cfg['auxiliary_data_dir']

    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info(
            'metadata filename:\t%s',
            metadata_filename,
        )

        metadatas = diagtools.get_input_files(cfg, index=index)
        thresholds = diagtools.load_thresholds(cfg, metadatas)

        if thresholds:
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
            if thresholds:
                make_map_contour(cfg, metadatas[filename], filename)

            ######
            # Maps of individual model
            make_map_plots(cfg, metadatas[filename], filename)

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
