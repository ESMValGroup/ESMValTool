"""
Diagnostic tools:

This module contains several python tools used elsewhere by the ocean
diagnostics package.

This tool is part of the ocean diagnostic tools package in the ESMValTool.

Author: Lee de Mora (PML)
    ledm@pml.ac.uk
"""
import logging
import os
import sys
import yaml
import matplotlib
matplotlib.use('Agg')  # noqa

import matplotlib.pyplot as plt

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def folder(name):
    """
    Make a directory out of a string or list or strings.

    Take a string or a list of strings, convert it to a directory style,
    then make the folder and the string.
    Returns folder string and final character is always os.sep. ('/')
    """
    sep = os.sep
    if isinstance(name, list):
        name = os.sep.join(name)
    if name[-1] != sep:
        name = name + sep
    if os.path.exists(name) is False:
        os.makedirs(name)
        logger.info('Making new directory:\t%s', str(name))
    return name


def get_input_files(cfg, index=0):
    """
    Load input configuration file.

    Get a dictionary with input files from the metadata.yml files.
    """
    metadata_file = cfg['input_files'][index]
    with open(metadata_file) as input_file:
        metadata = yaml.safe_load(input_file)
    return metadata


def bgc_units(cube, name):
    """
    Convert the cubes into some friendlier units.

    This is because many CMIP standard units are not the standard units
    used by the BGC community (ie, Celsius is prefered over Kelvin, etc.)
    """
    new_units = ''

    if name in ['tos', 'thetao']:
        new_units = 'celsius'

    if name in [
            'no3',
    ]:
        new_units = 'mmol m-3'

    if name in [
            'chl',
    ]:
        new_units = 'mg m-3'

    if new_units != '':
        logger.info(' '.join(
            ["Changing units from",
             str(cube.units), 'to', new_units]))
        cube.convert_units(new_units)

    return cube


def timecoord_to_float(times):
    """
    Convert from time coordinate into decimal time.

    Takes an iris time coordinate and returns a list of floats.
    """
    dtimes = times.units.num2date(times.points)
    floattimes = []
    for dtime in dtimes:
        # TODO: it would be better to have a calendar dependent value
        # for daysperyear, as this is not accurate for 360 day calendars.
        daysperyear = 365.25
        floattime = dtime.year + dtime.dayofyr / daysperyear + dtime.hour / (
            24. * daysperyear)
        if dtime.minute:
            floattime += dtime.minute / (24. * 60. * daysperyear)
        floattimes.append(floattime)
    return floattimes


def add_legend_outside_right(plot_details, ax1, column_width=0.1):
    """
    Add a legend outside the plot, to the right.

    plot_details is a 2 level dict,
    where the first level is some key (which is hidden)
    and the 2nd level contains the keys:
        'c': color
        'lw': line width
        'label': label for the legend.
    ax1 is the axis where the plot was drawn.
    """
    #####
    # Create dummy axes:
    legend_size = len(plot_details.keys()) + 1
    ncols = int(legend_size / 25) + 1
    box = ax1.get_position()
    ax1.set_position(
        [box.x0, box.y0, box.width * (1. - column_width * ncols), box.height])

    # Add emply plots to dummy axis.
    for index in sorted(plot_details.keys()):

        plt.plot(
            [], [],
            c=plot_details[index]['c'],
            lw=plot_details[index]['lw'],
            ls=plot_details[index]['ls'],
            label=plot_details[index]['label'])

    legd = ax1.legend(
        loc='center left',
        ncol=ncols,
        prop={'size': 10},
        bbox_to_anchor=(1., 0.5))
    legd.draw_frame(False)
    legd.get_frame().set_alpha(0.)


def get_image_format(cfg, default='png'):
    """
    Load the image format from the global config file.

    Current tested options are svg, png.

    The cfg is the opened global config.
    The default format is used if no specific format is requested.
    The default is set in the user config.yml
    Individual diagnostics can set their own format which will
    supercede the main config.yml.
    """
    image_extention = default

    # Load format from config.yml and set it as default
    if 'output_file_type' in cfg.keys():
        image_extention = cfg['output_file_type']

    # Load format from config.yml and set it as default
    if 'image_format' in cfg.keys():
        image_extention = cfg['image_format']

    matplotlib_image_formats = plt.gcf().canvas.get_supported_filetypes()
    if image_extention not in matplotlib_image_formats:
        logger.warning(' '.join(['Image format ', image_extention,
                                 'not in matplot:',
                                 ', '.join(matplotlib_image_formats)]))

    image_extention = '.' + image_extention
    image_extention = image_extention.replace('..', '.')
    return image_extention


def get_image_path(cfg,
                   metadata,
                   prefix='diag',
                   suffix='image',
                   metadata_id_list='default',):
    """
    Produce a path to the final location of the image.

    The cfg is the opened global config,
    metadata is the metadata dictionairy (for the individual dataset file)
    """
    #####
    if metadata_id_list == 'default':
        metadata_id_list = ['project', 'dataset', 'mip', 'exp', 'ensemble',
                            'field', 'short_name', 'preprocessor',
                            'diagnostic', 'start_year', 'end_year', ]

    path = folder(cfg['plot_dir'])
    if prefix:
        path += prefix + '_'
    # Check that the keys are in the dict.
    intersection = [va for va in metadata_id_list if va in metadata.keys()]
    path += '_'.join([str(metadata[b]) for b in intersection])
    if suffix:
        path += '_' + suffix

    image_extention = get_image_format(cfg)

    if path.find(image_extention) == -1:
        path += image_extention

    logger.info("Image path will be: %s", path)
    return path


def make_cube_layer_dict(cube):
    """
    Take a cube and return a dictionairy layer:cube

    Each item in the dict is a layer with a separate cube for each layer.
    ie:
        cubes[depth] = cube from specific layer

    Cubes with no depth component are returns as:
        cubes[''] = cube with no depth component.
    """
    #####
    # Check layering:
    depth = cube.coords('depth')
    cubes = {}

    if depth == []:
        cubes[''] = cube
    else:
        # iris stores coords as a list with one entry:
        depth = depth[0]
        if len(depth.points) in [
                1,
        ]:
            cubes[''] = cube
        else:
            coord_dim = cube.coord_dims('depth')[0]
            for layer_index, layer in enumerate(depth.points):
                slices = [slice(None) for index in cube.shape]
                slices[coord_dim] = layer_index
                cubes[layer] = cube[tuple(slices)]
    return cubes
