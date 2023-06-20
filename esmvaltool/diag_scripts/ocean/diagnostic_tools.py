"""
Diagnostic tools
================

This module contains several python tools used elsewhere by the ocean
diagnostics package.

This tool is part of the ocean diagnostic tools package in the ESMValTool.

Author: Lee de Mora (PML)
    ledm@pml.ac.uk
"""
import logging
import os
import sys
from pathlib import Path
import iris

import numpy as np
import cftime
import matplotlib.pyplot as plt
import yaml

from esmvaltool.diag_scripts.shared._base import _get_input_data_files

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
logging.getLogger('matplotlib.font_manager').disabled = True


ipcc_colours={
    'historical':'blue',
    'ssp126': 'green',
    'ssp245': 'gold',
    'ssp370': 'orange',
    'ssp585': 'red',}
CMIP5_blue = '#2551cc'
CMIP6_red = '#cc2323' # (204, 35, 35)
histnat_green= '#004F00' # 0,79,0
historical_beige = '#c47900' #



def get_obs_projects():
    """
    Return a list of strings with the names of observations projects.

    Please keep this list up to date, or replace it with something more
    sensible.

    Returns
    ---------
    list
        Returns a list of strings of the various types of observational data.
    """
    obs_projects = [
        'obs4MIPs',
    ]
    return obs_projects


def folder(name):
    """
    Make a directory out of a string or list or strings.

    Take a string or a list of strings, convert it to a directory style,
    then make the folder and the string.
    Returns folder string and final character is always os.sep. ('/')

    Arguments
    ---------
    name: list or string
        A list of nested directories, or a path to a directory.

    Returns
    ---------
    str
        Returns a string of a full (potentially new) path of the directory.
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


def get_input_files(cfg, index=''):
    """
    Load input configuration file as a Dictionairy.

    Get a dictionary with input files from the metadata.yml files.
    This is a wrappper for the _get_input_data_files function from
    diag_scripts.shared._base.

    Arguments
    ---------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.
    index: int
        the index of the file in the cfg file.

    Returns
    ---------
    dict
        A dictionary of the input files and their linked details.
    """
    if isinstance(index, int):
        metadata_file = cfg['input_files'][index]
        with open(metadata_file) as input_file:
            metadata = yaml.safe_load(input_file)
        return metadata
    return _get_input_data_files(cfg)


def bgc_units(cube, name):
    """
    Convert the cubes into some friendlier units.

    This is because many CMIP standard units are not the standard units
    used by the BGC community (ie, Celsius is prefered over Kelvin, etc.)

    Parameters
    ----------
    cube: iris.cube.Cube
        the opened dataset as a cube.
    name: str
        The string describing the data field.

    Returns
    -------
    iris.cube.Cube
        the cube with the new units.
    """
    new_units = ''
    if name in ['tos', 'thetao']:
        new_units = 'celsius'

    if name in ['sos', 'so']:
        cube.units = '0.001'

    if name in ['no3', 'o2', 'po4', 'si', 'dissic', 'talk']:
        new_units = 'mmol m-3'

    if name in ['chl', ]:
        new_units = 'mg m-3'

    if name in ['intpp', ]:
        new_units = 'mol m-2 yr-1'

    if name in ['fgco2', ]:
        new_units = 'g m-2 d-1'

    if name in ['spco2', 'dpco2', ]:
        new_units = 'uatm'

    if name in ['mfo', 'amoc', 'msftmyz']:
        # sverdrup are 1000000 m3.s-1, but mfo is kg s-1.
        new_units = 'Tg s-1'

    if name in ['zos', 'zostoga']:
        new_units = 'mm'

    if new_units != '':
        logger.info(' '.join(
            ["Changing units from",
             str(cube.units), 'to', new_units]))
        cube.convert_units(new_units)

    return cube


def match_model_to_key(
        model_type,
        cfg_dict,
        input_files_dict,
):
    """
    Match up model or observations dataset dictionairies from config file.

    This function checks that the control_model, exper_model and
    observational_dataset dictionairies from the recipe are matched with the
    input file dictionary in the cfg metadata.

    Arguments
    ---------
    model_type: str
        The string model_type to match (only used in debugging).
    cfg_dict: dict
        the config dictionary item for this model type, parsed directly from
        the diagnostics/ scripts, part of the recipe.
    input_files_dict: dict
        The input file dictionary, loaded directly from the get_input_files()
         function, in diagnostics_tools.py.

    Returns
    ---------
    dict
        A dictionary of the input files and their linked details.
    """
    for input_file, intput_dict in input_files_dict.items():
        intersect_keys = intput_dict.keys() & cfg_dict.keys()
        match = True
        for key in intersect_keys:
            if intput_dict[key] == cfg_dict[key]:
                continue
            match = False
        if match:
            return input_file
    logger.warning("Unable to match model: %s", model_type)
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


def cube_time_to_float(cube):
    """
    Convert from time coordinate into decimal time.

    Takes an iris time coordinate and returns a list of floats.
    Parameters
    ----------
    cube: iris.cube.Cube
        the opened dataset as a cube.

    Returns
    -------
    list
        List of floats showing the time coordinate in decimal time.

    """
    times = cube.coord('time')
    datetime = guess_calendar_datetime(cube)

    dtimes = times.units.num2date(times.points)
    floattimes = []
    for dtime in dtimes:
        # TODO: it would be better to have a calendar dependent value
        # for daysperyear, as this is not accurate for 360 day calendars.
        daysperyear = 365.25

        try:
            dayofyr = dtime.dayofyr
        except AttributeError:
            time = datetime(dtime.year, dtime.month, dtime.day)
            time0 = datetime(dtime.year, 1, 1, 0, 0)
            dayofyr = (time - time0).days

        floattime = dtime.year + dayofyr / daysperyear + dtime.hour / (
            24. * daysperyear)
        if dtime.hour:
            floattime += dtime.hour / (24. * daysperyear)
        if dtime.minute:
            floattime += dtime.minute / (24. * 60. * daysperyear)
        floattimes.append(floattime)
    return floattimes


def guess_calendar_datetime(cube):
    """
    Guess the cftime.datetime form to create datetimes.

    Parameters
    ----------
    cube: iris.cube.Cube
        the opened dataset as a cube.

    Returns
    -------
    cftime.datetime
        A datetime creator function from cftime, based on the cube's calendar.
    """
    time_coord = cube.coord('time')

    if time_coord.units.calendar in [
            '360_day',
    ]:
        datetime = cftime.Datetime360Day
    elif time_coord.units.calendar in ['365_day', 'noleap']:
        datetime = cftime.DatetimeNoLeap
    elif time_coord.units.calendar in [
            'julian',
    ]:
        datetime = cftime.DatetimeJulian
    elif time_coord.units.calendar in [
            'gregorian',
    ]:
        datetime = cftime.DatetimeGregorian
    elif time_coord.units.calendar in [
            'proleptic_gregorian',
    ]:
        datetime = cftime.DatetimeProlepticGregorian
    else:
        logger.warning('Calendar set to Gregorian, instead of %s',
                       time_coord.units.calendar)
        datetime = cftime.DatetimeGregorian
    return datetime


def get_decade(coord, value):
    """
    Determine the decade.

    Called by iris.coord_categorisation.add_categorised_coord.
    """
    date = coord.units.num2date(value)
    return date.year - date.year % 10


def decadal_average(cube):
    """
    Calculate the decadal_average.

    Parameters
    ----------
    cube: iris.cube.Cube
        The input cube

    Returns
    -------
    iris.cube
    """
    iris.coord_categorisation.add_categorised_coord(cube, 'decade', 'time',
                                                    get_decade)
    return cube.aggregated_by('decade', iris.analysis.MEAN)


def load_thresholds(cfg, metadata):
    """
    Load the thresholds for contour plots from the config files.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.
    metadata: dict
        the metadata dictionary

    Returns
    -------
    list:
        List of thresholds
    """
    thresholds = set()

    if 'threshold' in cfg:
        thresholds.update(float(cfg['threshold']))

    if 'threshold' in metadata:
        thresholds.update(float(metadata['threshold']))

    if 'thresholds' in cfg:
        thresholds.update([float(thres) for thres in cfg['thresholds']])

    if 'thresholds' in metadata:
        thresholds.update([float(thres) for thres in metadata['thresholds']])

    return sorted(list(thresholds))


def get_colour_from_cmap(number, total, cmap='jet'):
    """
    Get a colour `number` of `total` from a cmap.

    This function is used when several lines are created evenly along a
    colour map.

    Parameters
    ----------
    number: int, float
        The
    total: int

    cmap: string,  plt.cm
        A colour map, either by name (string) or from matplotlib
    """
    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)

    if number > total:
        raise ValueError(f'The cannot be larger than the total length '
                         f' of the list ie: {number} > {total}')

    if total > 1:
        colour = cmap(float(number) / float(total - 1.))
    else:
        colour = cmap(0.)
    return colour


def add_legend_outside_right(plot_details, ax1, column_width=0.1, loc='right'):
    """
    Add a legend outside the plot, to the right.

    plot_details is a 2 level dict,
    where the first level is some key (which is hidden)
    and the 2nd level contains the keys:
    'c': color
    'lw': line width
    'label': label for the legend.
    ax1 is the axis where the plot was drawn.

    Parameters
    ----------
    plot_details: dict
        A dictionary of the plot details (color, linestyle, linewidth, label)
    ax1: matplotlib.pyplot.axes
        The pyplot axes to add the
    column_width: float
        The width of the legend column. This is used to adjust for longer words
        in the legends
    loc: string
       Location of the legend. Options are "right" and "below".

    Returns
    -------
    cftime.datetime
        A datetime creator function from cftime, based on the cube's calendar.

    """
    # ####
    # Create dummy axes:
    legend_size = len(plot_details) + 1
    box = ax1.get_position()
    if loc.lower() == 'right':
        nrows = 25
        ncols = int(legend_size / nrows) + 1
        ax1.set_position([
            box.x0, box.y0, box.width * (1. - column_width * ncols), box.height
        ])

    if loc.lower() == 'below':
        ncols = 4
        nrows = int(legend_size / ncols) + 1
        ax1.set_position([
            box.x0, box.y0 + (nrows * column_width), box.width,
            box.height - (nrows * column_width)
        ])

    # Add emply plots to dummy axis.
    for index in sorted(plot_details):
        colour = plot_details[index]['c']

        linewidth = plot_details[index].get('lw', 1)

        linestyle = plot_details[index].get('ls', '-')

        label = plot_details[index].get('label', str(index))

        plt.plot([], [], c=colour, lw=linewidth, ls=linestyle, label=label)

    if loc.lower() == 'right':
        legd = ax1.legend(loc='center left',
                          ncol=ncols,
                          prop={'size': 10},
                          bbox_to_anchor=(1., 0.5))
    if loc.lower() == 'below':
        legd = ax1.legend(loc='upper center',
                          ncol=ncols,
                          prop={'size': 10},
                          bbox_to_anchor=(0.5, -2. * column_width))
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

    Arguments
    ---------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.

    Returns
    ---------
    str
        The image format extention.
    """
    image_extention = default

    # Load format from config.yml and set it as default
    if 'output_file_type' in cfg:
        image_extention = cfg['output_file_type']

    # Load format from config.yml and set it as default
    if 'image_format' in cfg:
        image_extention = cfg['image_format']

    matplotlib_image_formats = plt.gcf().canvas.get_supported_filetypes()
    if image_extention not in matplotlib_image_formats:
        logger.warning(' '.join([
            'Image format ', image_extention, 'not in matplot:',
            ', '.join(matplotlib_image_formats)
        ]))

    image_extention = '.' + image_extention
    image_extention = image_extention.replace('..', '.')
    return image_extention


def get_image_path(
        cfg,
        metadata,
        prefix='diag',
        suffix='image',
        metadata_id_list='default',
):
    """
    Produce a path to the final location of the image.

    The cfg is the opened global config,
    metadata is the metadata dictionary (for the individual dataset file)

    Arguments
    ---------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.
    metadata: dict
        The metadata dictionary for a specific model.
    prefix: str
        A string to prepend to the image basename.
    suffix: str
        A string to append to the image basename
    metadata_id_list: list
        A list of strings to add to the file path. It loads these from the cfg.

    Returns
    ---------
    str
        The ultimate image path

    """
    #####
    if metadata_id_list == 'default':
        metadata_id_list = [
            'project',
            'dataset',
            'mip',
            'exp',
            'ensemble',
            'field',
            'short_name',
            'preprocessor',
            'diagnostic',
            'start_year',
            'end_year',
        ]

    path = folder(cfg['plot_dir'])
    if prefix:
        path += prefix + '_'
    # Check that the keys are in the dict.
    intersection = [va for va in metadata_id_list if va in metadata]
    path += '_'.join([str(metadata[b]) for b in intersection])
    if suffix:
        path += '_' + suffix

    image_extention = get_image_format(cfg)

    if path.find(image_extention) == -1:
        path += image_extention

    path = path.replace(' ', '_')

    logger.info("Image path will be: %s", path)
    return path


def make_cube_layer_dict(cube):
    """
    Take a cube and return a dictionary layer:cube

    Each item in the dict is a layer with a separate cube for each layer.
    ie: cubes[depth] = cube from specific layer

    Cubes with no depth component are returned as dict, where the dict key
    is a blank empty string, and the value is the cube.

    Parameters
    ----------
    cube: iris.cube.Cube
        the opened dataset as a cube.

    Returns
    ---------
    dict
        A dictionary of layer name : layer cube.
    """
    #####
    # Check layering:
    coords = cube.coords()
    layers = []
    for coord in coords:
        if coord.standard_name in ['depth', 'region']:
            layers.append(coord)

    cubes = {}
    if not layers:
        cubes[''] = cube
        return cubes

    # if len(layers) > 1:
    #     # This field has a strange number of layer dimensions.
    #     # depth and regions?
    #     print(cube)
    #     raise ValueError('This cube has both `depth` & `region` coordinates:'
    #                      ' %s', layers)

    # iris stores coords as a list with one entry:
    layer_dim = layers[0]
    if len(layer_dim.points) in [
            1,
    ]:
        cubes[''] = cube
        return cubes

    if layer_dim.standard_name == 'depth':
        coord_dim = cube.coord_dims('depth')[0]
        for layer_index, layer in enumerate(layer_dim.points):
            slices = [slice(None) for index in cube.shape]
            slices[coord_dim] = layer_index
            cubes[layer] = cube[tuple(slices)]

    if layer_dim.standard_name == 'region':
        coord_dim = cube.coord_dims('region')[0]
        for layer_index, layer in enumerate(layer_dim.points):
            slices = [slice(None) for index in cube.shape]
            slices[coord_dim] = layer_index
            layer = layer.replace('_', ' ').title()
            cubes[layer] = cube[tuple(slices)]
    return cubes


def get_cube_range(cubes):
    """
    Determinue the minimum and maximum values of a list of cubes.

    Parameters
    ----------
    cubes: list of iris.cube.Cube
        A list of cubes.

    Returns
    ----------
    list:
        A list of two values: the overall minumum and maximum values of the
        list of cubes.

    """
    mins = []
    maxs = []
    for cube in cubes:
        mins.append(cube.data.min())
        maxs.append(cube.data.max())
    return [
        np.min(mins),
        np.max(maxs),
    ]


def get_cube_range_diff(cubes):
    """
    Determinue the largest deviation from zero in an list of cubes.

    Parameters
    ----------
    cubes: list of iris.cube.Cube
        A list of cubes.

    Returns
    ----------
    list:
        A list of two values: the maximum deviation from zero and its opposite.
    """
    ranges = []
    for cube in cubes:
        ranges.append(np.abs(cube.data.min()))
        ranges.append(np.abs(cube.data.max()))
    return [-1. * np.max(ranges), np.max(ranges)]


def get_array_range(arrays):
    """
    Determinue the minimum and maximum values of a list of arrays..

    Parameters
    ----------
    arrays: list of numpy.array
        A list of numpy.array.

    Returns
    ----------
    list:
        A list of two values, the overall minumum and maximum values of the
        list of cubes.
    """
    mins = []
    maxs = []
    for arr in arrays:
        mins.append(arr.min())
        maxs.append(arr.max())
    logger.info('get_array_range: %s, %s', np.min(mins), np.max(maxs))
    return [np.min(mins), np.max(maxs), ]


def make_mean_of_cube_list(cube_list, operation='mean'):
    """
    Takes the mean of a list of cubes (not an iris.cube.CubeList).
    Assumes all the cubes are the same shape.
    """
    # Fix empty times
    full_times = {}
    times = []
    for cube in cube_list:
        # make time coords uniform:

        cube.coord('time').long_name='Time axis'
        cube.coord('time').attributes={'time_origin': '1950-01-01 00:00:00'}
        times.append(cube.coord('time').points)

        for time in cube.coord('time').points:
            print(cube.name, time, cube.coord('time').units)
            try:
                full_times[time] += 1
            except:
                full_times[time] = 1

    for t, v in sorted(full_times.items()):
        if v != len(cube_list):
            print('FAIL', t, v, '!=', len(cube_list),'\nfull times:',  full_times)
            assert 0

    cube_mean=cube_list[0]
    #try: iris.coord_categorisation.add_year(cube_mean, 'time')
    #except: pass
    #try: iris.coord_categorisation.add_month(cube_mean, 'time')
    #except: pass

    try: cube_mean.remove_coord('year')
    except: pass
    #cube.remove_coord('Year')
    try: model_name = cube_mean.metadata[4]['source_id']
    except: model_name = ''
    print(model_name,  cube_mean.coord('time'))

    for i, cube in enumerate(cube_list[1:]):
        #try: iris.coord_categorisation.add_year(cube, 'time')
        #except: pass
        #try: iris.coord_categorisation.add_month(cube, 'time')
        #except: pass
        try: cube.remove_coord('year')
        except: pass
        #cube.remove_coord('Year')
        try: model_name = cube_mean.metadata[4]['source_id']
        except: model_name = ''
        print(i, model_name, cube.coord('time'))
        cube_mean+=cube
        #print(cube_mean.coord('time'), cube.coord('time'))
    cube_mean = cube_mean/ float(len(cube_list))
    return cube_mean


def make_mean_of_cube_list_notime(cube_list):
    """
    Takes the mean of a list of cubes (not an iris.cube.CubeList).
    Assumes all the cubes are the same shape.
    """
    # Fix empty times
    cube_mean=cube_list[0]
    #try: iris.coord_categorisation.add_year(cube_mean, 'time')
    #except: pass
    #try: iris.coord_categorisation.add_month(cube_mean, 'time')
    #except: pass

    try: cube_mean.remove_coord('year')
    except: pass
    #cube.remove_coord('Year')
    try: model_name = cube_mean.metadata[4]['source_id']
    except: model_name = ''

    #cube_mean = fix_depth(cube_mean)

    for i, cube in enumerate(cube_list[1:]):
        #try: iris.coord_categorisation.add_year(cube, 'time')
        #except: pass
        #try: iris.coord_categorisation.add_month(cube, 'time')
        #except: pass
        #cube = fix_depth(cube)

        try: cube.remove_coord('year')
        except: pass
        #cube.remove_coord('Year')
        try: model_name = cube_mean.metadata[4]['source_id']
        except: model_name = ''
        print(i, model_name)
        cube_mean+=cube
        #print(cube_mean.coord('time'), cube.coord('time'))
    cube_mean = cube_mean/ float(len(cube_list))
    return cube_mean


def prepare_provenance_record(cfg, **provenance_record):
    """
    Prepare informations to feed provenance

    Arguments
    ---------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.
    provenance_record: dict
        dictionary for a specific diagnostic provenance details.
    """
    recipe_path = Path(cfg['run_dir']).parents[1] / cfg['recipe']
    with recipe_path.open() as recipe_file:
        recipe = yaml.safe_load(recipe_file)

    doc = recipe['documentation']
    authors = doc.get('authors', [])
    authors += [
        maintainer for maintainer in doc.get('maintainer', [])
        if maintainer not in authors
    ]
    provenance_record['authors'] = authors
    for key in ['title', 'description', 'projects']:
        val = doc[key]
        if val:
            provenance_record[key] = val
    for key in ['realms', 'themes']:
        val = cfg.get(key)
        if val:
            provenance_record[key] = val
    return provenance_record
