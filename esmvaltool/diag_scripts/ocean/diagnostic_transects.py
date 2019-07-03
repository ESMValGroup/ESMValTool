"""
Transects diagnostics
=====================

Diagnostic to produce images of a transect. These plost show either latitude or
longitude against depth, and the cube value is used as the colour scale.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has no time component, and one of the latitude or
longitude coordinates has been reduced to a single value.

An approproate preprocessor for a 3D+time field would be::

  preprocessors:
    prep_transect:
      climate_statistics:
        operator: mean
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
from itertools import product

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def titlify(title):
    """
    Check whether a title is too long then add it to current figure.

    Parameters
    ----------
    title: str
        The title for the figure.
    """
    cutoff = 40
    if len(title) > cutoff:
        # Find good mid point
        titles = title.split(' ')
        length = 0
        for itr, word in enumerate(titles):
            length += len(word)
            if length > cutoff:
                titles[itr] += '\n'
                length = 0.
        title = ' '.join(titles)
    plt.title(title)


def determine_transect_str(cube, region=''):
    """
    Determine the Transect String.

    Takes a guess at a string to describe the transect.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube to use to determine the transect name.

    """
    if region:
        return region

    options = ['latitude', 'longitude']
    cube_dims = [c.standard_name for c in cube.coords()]
    for option in options:
        if option not in cube_dims:
            continue
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


def make_depth_safe(cube):
    """
    Make the depth coordinate safe.

    If the depth coordinate has a value of zero or above, we replace the
    zero with the average point of the first depth layer.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube to make the depth coordinate safe

    Returns
    ----------
    iris.cube.Cube:
        Output cube with a safe depth coordinate

    """
    depth = cube.coord('depth')

    # it's fine
    if depth.points.min() * depth.points.max() > 0.:
        return cube

    if depth.attributes['positive'] != 'down':
        raise Exception('The depth field is not set up correctly')

    depth_points = []
    bad_points = depth.points <= 0.
    for itr, point in enumerate(depth.points):
        if bad_points[itr]:
            depth_points.append(depth.bounds[itr, :].mean())
        else:
            depth_points.append(point)

    cube.coord('depth').points = depth_points
    return cube


def make_cube_region_dict(cube):
    """
    Take a cube and return a dictionairy region: cube.

    Each item in the dict is a layer with a separate cube for each layer.
    ie: cubes[region] = cube from specific region

    Cubes with no region component are returns as:
    cubes[''] = cube with no region component.

    This is based on the method diagnostics_tools.make_cube_layer_dict,
    however, it wouldn't make sense to look for depth layers here.

    Parameters
    ----------
    cube: iris.cube.Cube
        the opened dataset as a cube.

    Returns
    ---------
    dict
        A dictionairy of layer name : layer cube.
    """
    #####
    # Check layering:
    coords = cube.coords()
    layers = []
    for coord in coords:
        if coord.standard_name in ['region', ]:
            layers.append(coord)

    cubes = {}
    if layers == []:
        cubes[''] = cube
        return cubes

    # iris stores coords as a list with one entry:
    layer_dim = layers[0]
    if len(layer_dim.points) in [1, ]:
        cubes[''] = cube
        return cubes

    if layer_dim.standard_name == 'region':
        coord_dim = cube.coord_dims('region')[0]
        for layer_index, layer in enumerate(layer_dim.points):
            slices = [slice(None) for index in cube.shape]
            slices[coord_dim] = layer_index
            layer = layer.replace('_', ' ').title()
            cubes[layer] = cube[tuple(slices)]
    return cubes


def determine_set_y_logscale(cfg, metadata):
    """
    Determine whether to use a log scale y axis.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.

    Returns
    ----------
    bool:
        Boolean to flag whether to plot as a log scale.
    """
    set_y_logscale = True

    if 'set_y_logscale' in cfg:
        set_y_logscale = cfg['set_y_logscale']

    if 'set_y_logscale' in metadata:
        set_y_logscale = metadata['set_y_logscale']

    return set_y_logscale


def make_transects_plots(
        cfg,
        metadata,
        filename,
):
    """
    Make a simple plot of the transect for an indivudual model.

    This tool loads the cube from the file, checks that the units are
    sensible BGC units, checks for layers, adjusts the titles accordingly,
    determines the ultimate file name and format, then saves the image.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.
    filename: str
        The preprocessed model file.

    """
    # Load cube and set up units
    cube = iris.load_cube(filename)
    cube = diagtools.bgc_units(cube, metadata['short_name'])

    # Is this data is a multi-model dataset?
    multi_model = metadata['dataset'].find('MultiModel') > -1

    cube = make_depth_safe(cube)
    cubes = make_cube_region_dict(cube)

    # Determine y log scale.
    set_y_logscale = determine_set_y_logscale(cfg, metadata)

    for region, cube in cubes.items():
        # Make a dict of cubes for each layer.
        qplt.contourf(cube, 15, linewidth=0, rasterized=True)

        if set_y_logscale:
            plt.axes().set_yscale('log')

        if region:
            region_title = region
        else:
            region_title = determine_transect_str(cube, region)

        # Add title to plot
        title = ' '.join(
            [metadata['dataset'], metadata['long_name'], region_title])
        titlify(title)

        # Load image format extention
        image_extention = diagtools.get_image_format(cfg)

        # Determine image filename:
        if multi_model:
            path = diagtools.folder(
                cfg['plot_dir']) + os.path.basename(filename).replace(
                    '.nc', region + '_transect' + image_extention)
        else:
            path = diagtools.get_image_path(
                cfg,
                metadata,
                suffix=region + 'transect' + image_extention,
            )

        # Saving files:
        if cfg['write_plots']:
            logger.info('Saving plots to %s', path)
            plt.savefig(path)

        plt.close()


def add_sea_floor(cube):
    """
    Add a simple sea floor line from the cube mask.

    Parameters
    ----------
    cube: iris.cube.Cube
        Input cube to use to produce the sea floor.

    """
    land_cube = cube.copy()
    land_cube.data = np.ma.array(land_cube.data)
    mask = 1. * land_cube.data.mask
    if mask.shape == ():
        mask = np.zeros_like(land_cube.data)
    land_cube.data = np.ma.masked_where(mask == 0, mask)
    land_cube.data.mask = mask
    qplt.contour(land_cube, 2, cmap='Greys_r', rasterized=True)


def make_transect_contours(
        cfg,
        metadata,
        filename,
):
    """
    Make a contour plot of the transect for an indivudual model.

    This tool loads the cube from the file, checks that the units are
    sensible BGC units, checks for layers, adjusts the titles accordingly,
    determines the ultimate file name and format, then saves the image.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.
    filename: str
        The preprocessed model file.

    """
    # Load cube and set up units
    cube = iris.load_cube(filename)
    cube = diagtools.bgc_units(cube, metadata['short_name'])
    cube = make_depth_safe(cube)

    # Load threshold/thresholds.
    plot_details = {}
    colours = []
    thresholds = diagtools.load_thresholds(cfg, metadata)
    linewidths = [1 for thres in thresholds]
    linestyles = ['-' for thres in thresholds]

    cubes = make_cube_region_dict(cube)
    for region, cube in cubes.items():
        for itr, thres in enumerate(thresholds):
            colour = diagtools.get_colour_from_cmap(itr, len(thresholds))
            label = str(thres) + ' ' + str(cube.units)
            colours.append(colour)
            plot_details[thres] = {
                'c': colour,
                'lw': 1,
                'ls': '-',
                'label': label
            }

        qplt.contour(
            cube,
            thresholds,
            colors=colours,
            linewidths=linewidths,
            linestyles=linestyles,
            rasterized=True)

        # Determine y log scale.
        if determine_set_y_logscale(cfg, metadata):
            plt.axes().set_yscale('log')

        add_sea_floor(cube)

        # Add legend
        diagtools.add_legend_outside_right(
            plot_details, plt.gca(), column_width=0.08, loc='below')

        # Add title to plot
        title = ' '.join([
            metadata['dataset'], metadata['long_name'],
            determine_transect_str(cube, region)
        ])
        titlify(title)

        # Load image format extention
        image_extention = diagtools.get_image_format(cfg)

        # Determine image filename:
        if metadata['dataset'].find('MultiModel') > -1:
            path = diagtools.folder(
                cfg['plot_dir']) + os.path.basename(filename)
            path.replace('.nc', region + '_transect_contour' + image_extention)
        else:
            path = diagtools.get_image_path(
                cfg,
                metadata,
                suffix=region + 'transect_contour' + image_extention,
            )

        # Saving files:
        if cfg['write_plots']:
            logger.info('Saving plots to %s', path)
            plt.savefig(path)

        plt.close()


def multi_model_contours(
        cfg,
        metadatas,
):
    """
    Make a multi model comparison plot showing several transect contour plots.

    This tool loads several cubes from the files, checks that the units are
    sensible BGC units, checks for layers, adjusts the titles accordingly,
    determines the ultimate file name and format, then saves the image.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadatas: dict
        The metadatas dictionairy for a specific model.

    """
    ####
    # Load the data for each layer as a separate cube
    model_cubes = {}
    regions = {}
    thresholds = {}
    set_y_logscale = True

    for filename in sorted(metadatas):
        cube = iris.load_cube(filename)
        cube = diagtools.bgc_units(cube, metadatas[filename]['short_name'])
        cube = make_depth_safe(cube)
        cubes = make_cube_region_dict(cube)
        model_cubes[filename] = cubes
        for region in model_cubes[filename]:
            regions[region] = True

        # Determine y log scale.
        set_y_logscale = determine_set_y_logscale(cfg, metadatas[filename])

        # Load threshold/thresholds.
        tmp_thresholds = diagtools.load_thresholds(cfg, metadatas[filename])
        for threshold in tmp_thresholds:
            thresholds[threshold] = True

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Make a plot for each layer and each threshold
    for region, threshold in product(regions, thresholds):
        logger.info('plotting threshold: \t%s', threshold)
        title = ''
        plot_details = {}

        # Plot each file in the group
        for index, filename in enumerate(sorted(metadatas)):
            color = diagtools.get_colour_from_cmap(index, len(metadatas))
            linewidth = 1.
            linestyle = '-'
            # Determine line style for MultiModel statistics:
            if 'MultiModel' in metadatas[filename]['dataset']:
                linewidth = 2.
                linestyle = ':'
            # Determine line style for Observations
            if metadatas[filename]['project'] in diagtools.get_obs_projects():
                color = 'black'
                linewidth = 1.7
                linestyle = '-'

            qplt.contour(
                model_cubes[filename][region], [
                    threshold,
                ],
                colors=[
                    color,
                ],
                linewidths=linewidth,
                linestyles=linestyle,
                rasterized=True)

            plot_details[filename] = {
                'c': color,
                'ls': linestyle,
                'lw': linewidth,
                'label': metadatas[filename]['dataset']
            }

            if set_y_logscale:
                plt.axes().set_yscale('log')

            title = metadatas[filename]['long_name']
            units = str(model_cubes[filename][region].units)

            add_sea_floor(model_cubes[filename][region])

        # Add title, threshold, legend to plots
        title = ' '.join([
            title,
            str(threshold), units,
            determine_transect_str(model_cubes[filename][region], region)
        ])
        titlify(title)
        plt.legend(loc='best')

        # Saving files:
        if cfg['write_plots']:
            path = diagtools.get_image_path(
                cfg,
                metadatas[filename],
                prefix='MultipleModels',
                suffix='_'.join([
                    'contour_tramsect', region,
                    str(threshold) + image_extention
                ]),
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
    Load the config file and some metadata, then pass them the plot making
    tools.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.

    """
    #####
    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info(
            'metadata filename:\t%s',
            metadata_filename,
        )

        metadatas = diagtools.get_input_files(cfg, index=index)

        thresholds = diagtools.load_thresholds(cfg,
                                               next(iter(metadatas.values())))

        #######
        # Multi model contour plots
        if thresholds:
            multi_model_contours(
                cfg,
                metadatas,
            )

        for filename in sorted(metadatas):

            logger.info('-----------------')
            logger.info(
                'model filenames:\t%s',
                filename,
            )

            ######
            # Time series of individual model
            make_transects_plots(cfg, metadatas[filename], filename)

            ######
            # Contour maps of individual model
            if thresholds:
                make_transect_contours(cfg, metadatas[filename], filename)

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
