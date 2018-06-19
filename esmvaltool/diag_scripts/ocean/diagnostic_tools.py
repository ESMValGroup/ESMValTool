"""Python example diagnostic."""
import inspect
import logging
import os
import sys

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import yaml

from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def folder(name):
    """
        This snippet takes a string, makes the folder and the string.
        It also accepts lists of strings.
        """
    if isinstance(name, list):
        name = '/'.join(name)
    if name[-1] != '/':
        name = name + '/'
    if os.path.exists(name) is False:
        os.makedirs(name)
        logger.info('Making new directory:', name)
    return name


def get_input_files(cfg, index=0):
    """Get a dictionary with input files from metadata.yml files."""
    metadata_file = cfg['input_files'][index]
    with open(metadata_file) as file:
        metadata = yaml.safe_load(file)
    return metadata


def sensibleUnits(cube, name):
    """
    Convert the cubes into some friendlier units for the range of
    values typically seen in BGC.
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
            ["Changing units from ",
             str(cube.units), 'to', new_units]))
        cube.convert_units(new_units)

    return cube


def timecoord_to_float(times):
    """
       converts an iris time coordinate into a list of floats.
    """
    dtimes = times.units.num2date(times.points)
    floattimes = []
    daysperyear = 365.25
    for dt in dtimes:
        floattime = dt.year + dt.dayofyr / daysperyear + dt.hour / (
            24. * daysperyear)
        if dt.minute:
            floattime += dt.minute / (24. * 60. * daysperyear)
        floattimes.append(floattime)
    return floattimes


def add_legend_outside_right(plotDetails, ax1):
    """
           Add a legend outside the plot, to the right.
           PlotDetails is a 2 level dict,
           where the first level is some key (which is hidden)
           and the 2nd level contains the keys:
               'c': color
               'lw': line width
               'label': label for the legend.
           ax1 is the axis where the plot was drawn.
        """
    #####
    # Create dummy axes:
    legendSize = len(plotDetails.keys()) + 1
    ncols = int(legendSize / 25) + 1
    box = ax1.get_position()
    ax1.set_position(
        [box.x0, box.y0, box.width * (1. - 0.1 * ncols), box.height])

    # Add emply plots to dummy axis.
    for i in sorted(plotDetails.keys()):

        plt.plot(
            [], [],
            c=plotDetails[i]['c'],
            lw=plotDetails[i]['lw'],
            ls=plotDetails[i]['ls'],
            label=plotDetails[i]['label'])

    legd = ax1.legend(
        loc='center left',
        ncol=ncols,
        prop={'size': 10},
        bbox_to_anchor=(1., 0.5))
    legd.draw_frame(False)
    legd.get_frame().set_alpha(0.)


def get_image_path(cfg,
                   md,
                   prefix='',
                   suffix='',
                   image_extention='png',
                   basenamelist=[
                       'project', 'model', 'mip', 'exp', 'ensemble', 'field',
                       'short_name', 'preprocessor', 'diagnostic',
                       'start_year', 'end_year'
                   ]):
    """
        This produces a path to the final location of the image.
        The cfg is the opened global config,
        md is the metadata dictionairy (for the individual model file)
        """

    path = folder(cfg['plot_dir'])
    if prefix:
        path += prefix + '_'
    path += '_'.join([str(md[b]) for b in basenamelist])
    if suffix:
        path += '_' + suffix
    path += '.' + image_extention
    logger.info("Image path will be: %s", path)
    return path


def make_cube_layer_dict(cube):
    """
        This method takes a cube and return a dictionairy
        with a cube for each layer as it's own item. ie:
          cubes[depth] = cube from specific layer
        Also, cubes with no depth component are returns as:
          cubes[''] = cube with no depth component.
        """

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
            for l, layer in enumerate(depth.points):
                cubes[layer] = cube[:, l]
    return cubes
