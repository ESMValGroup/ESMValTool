#! /usr/local/sci/bin/python

# This file contains extra routines needed for the validation note
# software but may also be useful elsewhere

import csv
import os
import pdb
import pickle
import string
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

import cf_units
import iris
from iris.analysis.maths import intersection_of_cubes as ioc

import globalvar
import rms

# Import any routines in general/cma_python/cubeinout.py for use in the
#  equation
import cubeinout

# Local variables
# A dictionary to hold regridders (speeding up regridding operations)
regridded_dict = {}
# Accepetable chars in equations
acceptable_char_list = ['+', '-', '*', '/', '^', '(', ')', ',']


def read_control_file(controlfile, ctldir=None):
    """
    Reads the control files and returns the data held in them in
    two nested dictionaries
    """

    # Initialise variables
    in_sec = 0
    large_dict = {}

    # Open the file
    if ctldir:
        controlfile = os.path.join(ctldir, controlfile)
    with open(controlfile, 'r') as f:

        # Read each line in turn
        for line in f:

            # Split it up into substrings and remove whitespace
            line_split = string.split(line, ':')
            line_split = [s.replace('\n', '') for s in line_split]
            line_split_nw = [s.replace(' ', '') for s in line_split]

            if len(line_split) > 1:

                # Ignore any lines beginning with hash
                if line_split_nw[0][0] != '#':

                    # Work out if we are at the start or end of a section
                    if line_split_nw[0] == 'start':
                        in_sec = 1
                        section = line_split[1].strip()
                        setting_dict = {}

                    elif line_split_nw[0] == 'end':
                        in_sec = 0
                        # Add to the large dictionary
                        large_dict[section] = setting_dict

                    else:
                        if in_sec:
                            # Add to the small dictionary
                            setting_dict[line_split_nw[0]] = line_split[1:]

    return large_dict


def read_info_file(infofile):
    """
    Reads the info file and returns the data held in it as a single dictionary
    """

    # Initialise variables
    info_dict = {}

    # Read the valorder file as a csv file
    with open(infofile, 'r') as inf:
        f = csv.reader(inf, delimiter=':')
        for line in f:
            # Ignore any blank lines or those beginning with a hash
            if len(line) > 0:
                if line[0][0] != '#':
                    info_dict[line[0].strip()] = line[1].strip()

    return info_dict


def write_info_file(infofile, info_dict):
    """
    Writes information to an info file.

    infofile=(string) a filename to hold this information. The final file will
             have two columns of data: keys on the left and items on the right,
             separated by a colon.
    info_dict=(dictionary) the information to store in the file
    """

    # Open a file ready for reading
    with open(infofile, 'w') as f:
        # Loop over the info_dict
        for info_key, info_item in info_dict.iteritems():
            f.write(info_key+':'+info_item+'\n')


def key_valid(info_key):
    """Checks that info_key is a key that points to a file name"""
    acceptable_key_list = ['djfm', 'mamm', 'jjam', 'sonm', 'annm',
                           'djfs', 'mams', 'jjas', 'sons', 'anns']
    return info_key in acceptable_key_list


def load_pickle(pickle_file):
    """Loads from a pickle file"""
    with open(pickle_file) as f:
        cube_list = pickle.load(f)

    return cube_list


def save_pickle(cube_list, pickle_file):
    """Saves to a pickle file"""
    with open(pickle_file, 'wb') as f:
        pickle.dump(cube_list, f, protocol=-1)


def split_equation(equation):
    """Reads in an equation and returns the component parts"""

    cubeinout_contents = dir(cubeinout)

    # Initialise arrays
    components = []
    this_comp = ''

    # Loop over each charactor and check for maths symbols. Append the
    # strings to the components array.
    for char in equation:
        if char in acceptable_char_list:

            # If you have a decent length component append it to a list of
            #  components
            this_comp_no_space = this_comp.replace(' ', '')
            if this_comp_no_space != '':

                # Components can't be purely numerical
                try:
                    this_comp_float = float(this_comp_no_space)
                except ValueError:

                    # Components can't be in cubeinout
                    if this_comp_no_space not in cubeinout_contents:
                        components.append(this_comp_no_space)
            this_comp = ''
        else:
            this_comp = this_comp+char

    this_comp_no_space = this_comp.replace(' ', '')
    if this_comp_no_space != '':

        # Components can't be purely numerical
        try:
            this_comp_float = float(this_comp_no_space)
        except ValueError:

            # Components can't be in cubeinout
            if this_comp_no_space not in cubeinout_contents:
                components.append(this_comp_no_space)

    return components


def apply_heavyside(cube, heavyside_cube):
    """
    Uses the heavyside to mask out points where the field is below ground level

    cube = the input field
    heavyside_cube = the heavyside cube that contains 0 for data below ground,
                     1 for data above ground and fractions for partial
                     situations.
    """

    mask = heavyside_cube.data <= 0.1
    cube.data = np.ma.array(cube.data, mask=mask)
    cube.data /= heavyside_cube.data
    return cube


def expand_equation(equation, component_list, dict_name='data_regrid_dict'):
    """
    So far we have an equation but in order to execute it we need
    data_regrid_dict in front of each component in the equation. This routine
    adds those data_dict bits.

    equation = (str) the equation from the plots_file
    components = a list of components where each item is a string of the
                 form source_season_field
    dict_name = (str) the name of the dictionary that contains the
                 variables that you want in the equation. The default (as used
                 by the validation notes) is data_regrid_dict
    """

    cubeinout_contents = dir(cubeinout)

    # long_eqn will form the new equation with data_regrid_dict embedded
    #  within it
    long_eqn = equation

    # Loop backwards over the equation
    len_equation = len(equation)
    in_component = False
    for equation_index in range(len_equation-1, -1, -1):
        char = equation[equation_index]
        if char in acceptable_char_list or char == ' ':
            if in_component:
                start_component = equation_index
                in_component = False

                # Found a component. Now check that it is on the list.
                # If so put a wrapper round it.
                if end_component - start_component > 1:
                    component = equation[start_component+1:end_component+1]
                    if component in component_list:
                        long_eqn = string.join([long_eqn[0:start_component+1],
                                                dict_name,
                                                "['", component, "']",
                                                long_eqn[end_component+1:]],
                                               '')

                    # If it is part of cubeinout put cubeinout in front
                    elif component in cubeinout_contents:
                        long_eqn = string.join([long_eqn[0:start_component+1],
                                                'cubeinout.',
                                                component,
                                                long_eqn[end_component+1:]],
                                               '')
        else:
            if not in_component:
                end_component = equation_index
            in_component = True

    # Do the final component as it is not done in the loop
    if in_component:
        start_component = equation_index-1

        # Found a component. Now check that it is on the list.
        # If so put a wrapper round it.
        if end_component - start_component > 1:
            component = equation[start_component+1:end_component+1]
            if component in component_list:
                long_eqn = string.join([long_eqn[0:start_component+1],
                                        dict_name,
                                        "['", component, "']",
                                        long_eqn[end_component+1:]],
                                       '')

            # If it is part of cubeinout put cubeinout in front
            elif component in cubeinout_contents:
                long_eqn = string.join([long_eqn[0:start_component+1],
                                        'cubeinout.',
                                        component,
                                        long_eqn[end_component+1:]],
                                       '')

    return long_eqn


def get_cube_ready(cube):
    """
    Sets up cubes ready for regridding, in the equation and for
    generating means
    """

    # Remove unnecessary coordinates in cube
    coords_list = ['forecast_reference_time', 'forecast_period', 'source',
                   'season', 'time']
    for coord in coords_list:
        try:
            cube.remove_coord(coord)
        except iris.exceptions.CoordinateNotFoundError:
            pass

    # Also make sure we are not using grid_longitude and grid_latitude
    fix_coords(cube)

    # Test to see if the grid does not have a coordinate system (CERES-EBAF).
    # If so fix it.
    if cube.coord_system() is None:
        cube = fix_cube(cube)
        if globalvar.error:
            return 1

    # Test to see if the grid is a rotated coordinate system (HadSLP2).
    # If so fix it.
    geogcs = iris.coord_systems.GeogCS(6371229.0)
    rot_cs = iris.coord_systems.RotatedGeogCS(0.0, 0.0, ellipsoid=geogcs)
    # TODO: Should this be an 'is' instead of '=='
    if cube.coord_system() == rot_cs:
        cube = fix_cube(cube)
        if globalvar.error:
            return 1

    # Roll the cube but only if you have a longitude dimension and it is global
    coords_list = [this_coord.name() for this_coord in cube.coords()]
    if 'longitude' in coords_list or 'grid_longitude' in coords_list:
        points = cube.coord(axis='x').points
        if points[-1] - points[0] > 350:
            cube = roll_cube(cube)

    # Remove redundant multi-dimensional arrays
    if 'surface_altitude' in coords_list:
        cube.remove_coord('surface_altitude')

    # Guess the bounds
    if not cube.coord(axis='x').has_bounds():
        cube.coord(axis='x').guess_bounds()
    if not cube.coord(axis='y').has_bounds():
        cube.coord(axis='y').guess_bounds()

    return cube


def same_units(cube_dict):
    """Make sure all the cubes in this dictionary have the same units"""
    component_set = cube_dict.keys()
    for comp in component_set:
        cube_dict[comp].units = ''


def get_component_set(equation):
    """Get a set of components used in this equation"""

    component_list = []
    for comp in split_equation(equation):
        if comp != '':
            component_list.append(comp)
    component_set = set(component_list)

    return component_set


def fix_coords(incube):
    """
    Makes so that all grid_longitude axes are renamed longitude
    and all grid_latitude axes are renamed latitude.

    incube = input cube to fix
    returned value = fixed cube
    """

    dims = [cor.standard_name for cor in incube.coords()]
    if 'grid_longitude' in dims:
        incube.coord('grid_longitude').standard_name = 'longitude'
    if 'grid_latitude' in dims:
        incube.coord('grid_latitude').standard_name = 'latitude'
    if u'longitude' in dims:
        incube.coord(u'longitude').standard_name = 'longitude'
    if u'latitude' in dims:
        incube.coord(u'latitude').standard_name = 'latitude'


def fix_cube(incube):
    """Fixes this cube by putting it on a GeogCS coordinate system
    with standard latitudes and longitudes.

    incube = input cube to fix
    returned value = fixed cube
    """

    print 'Incompatible coordinate system detected. Fixing cube.'

    geogcs = iris.coord_systems.GeogCS(6371229.0)

    # Generate new latitudes and longitudes on a GeogCS coordinate system
    old_lats_dim = incube.coord(axis='y')
    new_lats_dim = iris.coords.DimCoord(old_lats_dim.points,
                                        standard_name='latitude',
                                        units=cf_units.Unit('degrees'),
                                        coord_system=geogcs)
    old_lons_dim = incube.coord(axis='x')
    new_lons_dim = iris.coords.DimCoord(old_lons_dim.points,
                                        standard_name='longitude',
                                        units=cf_units.Unit('degrees'),
                                        coord_system=geogcs)

    # Reshape the data to remove any extra, single length, dimensions.
    reshaped_data = np.reshape(incube.data, (new_lats_dim.shape[0],
                                             new_lons_dim.shape[0]))

    # Generate a new cube from the raw data and the new latitudes and
    #  longitudes
    outcube = iris.cube.Cube(reshaped_data,
                             dim_coords_and_dims=[(new_lats_dim, 0),
                                                  (new_lons_dim, 1)])

    return outcube


def regrid_cube(incube, newgrid, scheme_str='AreaWeighted', comp='none'):
    """Regrids a cube to a new grid, retaining missing data:

    incube = input cube to regrid
    rewgrid = new grid cube
    scheme_str = (str) The scheme by which to do the regridding. One of:
       'AreaWeighted' = Use the area weighted scheme (default)
       'Linear' = Use a linear interpolation scheme
    comp = (str) the component name. If not equal to 'none' then it will
           look for existing regridded versions of this component and use that.
    """

    # Assign the area weighted scheme for regridding
    if scheme_str == 'AreaWeighted':
        scheme = iris.analysis.AreaWeighted(mdtol=0.5)
    elif scheme_str == 'Linear':
        scheme = iris.analysis.Linear()
    else:
        print 'ERROR in valmod.regrid_cube: scheme_str value of ', \
            scheme_str, ' not supported'

    # Assign a regridder name
    incube_str = comp + scheme_str + str(incube.shape) + \
        str(incube.coords(axis='x')[0].points[0]) + \
        str(incube.coords(axis='y')[0].points[0])
    newgrid_str = comp + scheme_str + str(newgrid.shape) + \
        str(newgrid.coords(axis='x')[0].points[0]) + \
        str(newgrid.coords(axis='y')[0].points[0])
    if incube_str != newgrid_str:

        regrided_cube = None
        if comp != 'none':

            # See if this regridded cube already exists in memory
            if newgrid_str in regridded_dict.keys():
                print 'Using already regridded version of ', newgrid_str
                regrided_cube = regridded_dict[newgrid_str]

        if regrided_cube is None:
            print 'Regridding cube ', incube_str, ' to ', newgrid_str
            try:
                regrided_cube = incube.regrid(newgrid, scheme)
            except:
                print "ERROR: Cubes are not compatible for regridding."
                print "ERROR: ", sys.exc_info()
                if globalvar.debug:
                    pdb.set_trace()
                else:
                    globalvar.error = True
                    return 1

            # Save this regridded cube for later
            if comp != 'none':
                regridded_dict[newgrid_str] = regrided_cube

    else:

        # Cubes have the same grid. Simply copy across
        print 'Keeping cubes in the same grid: '+incube_str+', '+newgrid_str
        regrided_cube = incube

    return regrided_cube


def roll_cube(cube, start_longitude=-180.0):
    """
    Rolls the cube so that it starts from 180W and goes to 180E
    (instead of from 0E to 360E).

    cube = input cube that you wish to alter.
    start_longitude = destination starting longitude. Default is -180 (180W).
    """

    # Work out how much to shift the longitude axis by
    start_longitude_old = cube.coord(axis='x').points[0]
    degrees_to_shift = start_longitude - start_longitude_old
    fraction_to_shift = degrees_to_shift/360.0

    # However in reality we can only shift by a whole number of grid boxes.
    # Work out how many grid boxes to shift by
    long_axis = cube.coord_dims(cube.coord(axis='x'))[0]
    number_boxes_to_shift = int(cube.shape[long_axis]*fraction_to_shift)

    # Correct the degrees_to_shift to match number_boxes_to_shift
    degrees_to_shift = \
        360.0 * number_boxes_to_shift/float(cube.shape[long_axis])

    # Roll the data
    cube.data = np.roll(cube.data, number_boxes_to_shift, axis=long_axis)

    # Roll the coordinates
    cube.coord(axis='x').points = cube.coord(axis='x').points+degrees_to_shift

    # Roll the bounds
    if cube.coord(axis='x').has_bounds():
        cube.coord(axis='x').bounds = \
            cube.coord(axis='x').bounds+degrees_to_shift

    return cube


def area_avg(cube):
    """
    Perform an area average of a cube using weights to account for
    changes in latitude.
    """

    # What dimension is latitude
    thiscoords = cube.coords()
    lat_dim = -99
    for i in range(len(thiscoords)):
        if thiscoords[i].standard_name == 'latitude':
            lat_dim = i
    if lat_dim == -99:
        raise ValError('Could not find latitude coordinate in cube')

    # Convert latitudes to radians
    lats_radians_array = \
        cube.coord(standard_name='latitude').points/360.0 * 2.0 * np.pi

    # Generate the cos of these latitudes
    cos_lats_array = np.cos(lats_radians_array)

    # Make a 2d array of these cos values
    cos_2d_array = np.zeros(cube.data.shape)
    # cos_2d_array = np.empty_like(cube.data)
    for i in range(cos_2d_array.shape[0]):
        for j in range(cos_2d_array.shape[1]):
            got_data = True
            location = [i, j]
            lat_index = location[lat_dim]
            if hasattr(cube.data, 'mask'):
                if cube.data.mask.size == cube.data.size:
                    got_data = not cube.data.mask[i, j]
            if got_data:
                cos_2d_array[i, j] = cos_lats_array[lat_index]

    # Do some area averaging
    area_average = np.average(cube.data, weights=cos_2d_array)

    return area_average


def add_colorbar(mappable, units=None):
    # adds a colorbar to a plot with some standard keywords set
    cbar_kwargs = {'pad': 0.15, 'orientation': 'horizontal', 'extend': 'both'}

    cax, cbar_kwargs_filter = matplotlib.colorbar.make_axes(plt.gca(),
                                                            **cbar_kwargs)
    cb = plt.colorbar(mappable, cax=cax, **cbar_kwargs_filter)
    cb.ax.tick_params(length=0, labelsize='small')

    # Add the units
    if not globalvar.pub and 'units' in locals():
        cb.set_label(units, size='small')

    return cb


# TODO: This has already been done in Cartopy?
class LatsFormatter(matplotlib.ticker.Formatter):
    """
    This is a custom formatter that converts the native unit of degrees
    into North/South degrees.
    """
    def __call__(self, x, pos=None):
        letter = 'S' if x < 0 else 'N'
        return u"%d%s" % (abs(x), letter)


# TODO: This has already been done in Cartopy?
class LonsFormatter(matplotlib.ticker.Formatter):
    """
    This is a custom formatter that converts the native unit of degrees
    into West/East degrees.
    """
    def __call__(self, x, pos=None):
        letter = 'W' if x < 0 else 'E'
        return u"%d%s" % (abs(x), letter)


class ValError(Exception):
    """
    This is the default exception for errors picked up by the validation
    note software
    """

    def __init__(self, value):
        self.value = value
        print value

    def __str__(self):
        return repr(self.value)


def make_cmap(red_list, green_list, blue_list):
    """
    This will make a colour map given arrays of red, green and blue colours
    """

    # Check that the lengths of the red, green and blue arrays are the same
    if len(red_list) != len(green_list) or len(red_list) != len(blue_list):
        msg = 'The lengths of the red, green and blue colour tables are ' + \
            'not equal'
        raise ValError(msg)
    n_col = len(red_list)

    # Make a blank colour dictionary
    cdict = {}

    # Loop over the colours
    for colour in ['red', 'green', 'blue']:

        # Set up a generic list for each of the colours
        colour_list = eval(colour+'_list')

        # Check that there are no values greater than one
        if max(colour_list) > 1.0:
            msg = 'There is a colour greater than 1.0 in the ' + colour + \
                ' colour list'
            raise ValError(msg)

        # Make a blank colour list just for this colour
        one_cmap_list = []

        for c in range(n_col):

            # What is the array index (from 0 to 1)
            array_index = c/float(n_col-1)

            # Populate the one colour list
            one_cmap_list.append((array_index, colour_list[c], colour_list[c]))

        # Add this to the colour dictionary
        cdict[colour] = one_cmap_list

    # Make the colour map
    cmap = matplotlib.colors.LinearSegmentedColormap('default', cdict)

    return cmap


def file_name(page_num, page_title):
    """
    Make a file name derived from the page number and page title.
    This is done by replacing all spaces and in the page title with
    underscores and removing all brackets.
    """

    # Replace spaces with underscores in page_title
    page_title = page_title.replace(' ', '_')

    # Remove all brackets in page_title
    page_title = page_title.replace('(', '')
    page_title = page_title.replace(')', '')

    # Join page_num with page_title
    if globalvar.pub:
        filename = page_num+'_'+page_title+'.eps'
        # Note: pcolormesh does not work in eps format. Use pdf format and
        #       convert to eps using pdf2ps.
    else:
        filename = page_num+'_'+page_title+'.png'
    return filename


# TODO: Get rid of this routine!!!!
def mk_temp():
    """
    Make a temporary directory on your local disk to store pkl files and
    page titles.

    returned_value = full path to temporary directory.
    """

    localdata = os.getenv('LOCALDATA')
    temp_dir = os.path.join(localdata, 'valnote_temp')
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
    return(temp_dir)


def extract_no_diurnal(all_data_cube, all_constraints=None):
    """
    Extracts one cube from all_data_cube based on the constraints provided.
    It will remove any diurnal cycle diagnostics.

    all_data_cube = a cube list containing all the data
    all_constraints = a list of constraints needed to extract the data you
                      want. This routine will attempt to extract data from
                      each constraint sequentially until data is successfully
                      returned.
    """

    # Extract the data
    if all_constraints:

        # Loop over each constraint
        for this_constraint in all_constraints:

            # Attempt to extract the data
            print 'Extracting using constraint: '
            print this_constraint

            this_data_cube = all_data_cube.extract(this_constraint)

            # Test to see if you managed to get some data.
            # If so exit this for loop.
            cube_length = 1
            if isinstance(this_data_cube, iris.cube.CubeList):
                cube_length = len(this_data_cube)
            if cube_length != 0:
                break

    else:
        this_data_cube = all_data_cube

    # If this is a cube list then get it down to just one cube
    if isinstance(this_data_cube, iris.cube.CubeList):
        cube_length = len(this_data_cube)

        # Check to make sure we have some data
        if cube_length == 0:
            print 'ERROR: No cubes returned from all_data_cube.extract ' + \
                'command in extract_no_diurnal'
            if globalvar.debug:
                # Stop to query why we have no data..
                pdb.set_trace()
                globalvar.error = True
                return 1
            else:
                globalvar.error = True
                return 1

        # If there is just one cube in this list extract that
        if cube_length == 1:
            this_data_cube = this_data_cube[0]

        if cube_length > 1:
            print 'WARNING: detecting cube list in extract_no_diurnal.'
            print '         Using the first cube (but this may be wrong).'
            this_data_cube = this_data_cube[0]

    # Get a list of dimensions
    dims_name_list = []
    for this_dim in this_data_cube.coords():
        if type(this_dim.standard_name) in (str, unicode):
            dims_name_list.append(this_dim.standard_name)
        elif type(this_dim.long_name) in (str, unicode):
            dims_name_list.append(this_dim.long_name)
        else:
            print 'ERROR: No standard name or long name detected for ' + \
                'coordinate in extract_no_diurnal:'
            print this_dim
            if globalvar.debug:
                pdb.set_trace()
            else:
                globalvar.error = True
                return 1

    # Get the longitude, latitude and vertical dimension names
    lat_name = 'x'
    for dim in dims_name_list:
        if 'lat' in dim:
            lat_name = dim
    lon_name = 'x'
    for dim in dims_name_list:
        if 'lon' in dim:
            lon_name = dim
    vert_name = 'x'
    for dim in dims_name_list:
        if 'pressure' in dim or 'vertical' in dim or \
           'height' in dim or 'level' in dim:
            vert_name_temp = dim
            if this_data_cube.coords(vert_name_temp)[0].shape[0] > 1:
                vert_name = vert_name_temp
    pseudo_name = 'x'
    for dim in dims_name_list:
        if 'pseudo_level' in dim:
            pseudo_name = dim

    # Make sure we have found longitudes and latitudes
    if lat_name == 'x' or lon_name == 'x':
        print 'ERROR: could not detect either longitude or latitude axes'
        if globalvar.debug:
            pdb.set_trace()
        else:
            globalvar.error = True
            return 1

    # Initialise some variables for future use
    greatest_bounds_zero = 0.0
    smallest_bounds_two = 1.0e9

    # Assess whether this a two or three dimensional and assign a
    #  list of dimensions accordingly
    if pseudo_name != 'x':
        this_dims_list = [pseudo_name, vert_name, lat_name, lon_name]
    elif vert_name != 'x':
        this_dims_list = [vert_name, lat_name, lon_name]
    else:
        this_dims_list = [lat_name, lon_name]

    # Make a longitude/latitude slice and assign it as the output
    for subset_cube in this_data_cube.slices(this_dims_list):
        if 'output_cube' not in locals():
            output_cube = subset_cube

        # For model output data
        if len(subset_cube.aux_coords) == 3:
            if subset_cube.aux_coords[0].bounds[0][0] \
               >= greatest_bounds_zero \
               and subset_cube.aux_coords[2].bounds[0][0] \
               <= smallest_bounds_two:
                output_cube = subset_cube
                greatest_bounds_zero = subset_cube.aux_coords[0].bounds[0][0]
                smallest_bounds_two = subset_cube.aux_coords[2].bounds[0][0]

    # If the output cube doesn't exist then stop with an error
    if 'output_cube' not in locals():
        print 'ERROR: could not detect time mean diagnostic'
        if globalvar.debug:
            pdb.set_trace()
        else:
            globalvar.error = True
            return 1

    return output_cube


def load_extra_data(extras_dict):
    """Loads in extra data, such as land-sea masks, into a dictionary"""

    # Initialise variable
    globalvar.extra_data_dict = {}

    for key in extras_dict.iterkeys():
        globalvar.extra_data_dict[key] = iris.load_cube(extras_dict[key])


def plot_type_test(incube):
    """
    Detects what type of plot you want to do to display this cube.
    Returns one of:

    'lat_lon'    = Standard latitude and logitude contour plot
    'y_pressure' = Uses pressure as a y axis in contour plot
    'y_height'   = Uses height as a y axis in contour plot
    'x_longitude'= Uses longitude on the x axis for line plot
    'x_latitude' = Uses latitude on the x axis for line plot
    """

    # Make a list of coordinates with longer than 1 dimension.
    # Ignore auxillary coordinates.
    coords_list = []
    dimension_number_list = []
    for this_coord in incube.coords():
        name = this_coord.name()
        if len(this_coord.points) > 1:
            if incube.coord_dims(name)[0] \
               not in dimension_number_list:
                coords_list.append(name)
                dimension_number_list.append(incube.coord_dims(name)[0])

    # Return error if no dimensions
    if len(coords_list) == 0:
        print 'ERROR: zero dimension data detected'
        if globalvar.debug:
            pdb.set_trace()
        else:
            globalvar.error = True
            return 'error'

    # If one dimension then return one of the line plot names
    elif len(coords_list) == 1:
        plot_type = 'error'
        if 'longitude' in coords_list:
            plot_type = 'x_longitude'
        if 'latitude' in coords_list:
            plot_type = 'x_latitude'

    # If two dimensions then return one of the contour plot names
    elif len(coords_list) == 2:
        plot_type = 'lat_lon'
        if 'model_level_number' in coords_list:
            plot_type = 'y_height'
        if 'pressure' in coords_list:
            plot_type = 'y_pressure'

    else:
        print 'ERROR: can not plot cubes with more than 2 dimensions'
        if globalvar.debug:
            pdb.set_trace()
        else:
            globalvar.error = True
            return 'error'

    return plot_type


def perform_equation(equation, data_dict, page_title, key, rms_list=None,
                     filename=None, calc_rms=False, same_levels=True):
    """
    Perform the equations provided in plots_file.dat

    equation = the equation to perform
    data_dict = a dictionary containing the data needed for those equations
    page_title = the title of the page from valorder.dat
    key = the key (a, b, c or d) for the plot
    rms_list = an rms object for storing root mean square values of the final
               field (see rms.py for object structure)
    calc_rms = switch to turn on calculation of rms values
    filename = (str) the filename of the png plot
    same_levels = If multi level data make sure all your levels are the same
                  before you perform the equation (default = True)
    """

    print 'Performing equation '+equation

    # Generate a list of components used in the equation
    component_set = get_component_set(equation)

    # Are any of these components observations.
    using_obs = 0
    for comp in component_set:
        input_info = string.split(comp, '_')
        if input_info[0][0:5] != 'exper' and input_info[0][0:7] != 'control':
            using_obs = 1

        # Check that the component exists in data_dict
        if comp not in data_dict.keys():
            print 'ERROR: Not all variables in equation ' + equation + \
                ' have entries in data_dict. data_dict.keys()='
            print data_dict.keys()
            if globalvar.debug:
                pdb.set_trace()
            else:
                globalvar.error = True
                return 1, globalvar.missing_data

    # Make sure all the fields have correct units and axes
    for comp in component_set:
        data_dict[comp] = get_cube_ready(data_dict[comp])
        if globalvar.error:
            return 1, globalvar.missing_data

    # Get a list of longitude sizes
    n_lons_list = []
    for comp in component_set:
        n_lons_list.append(data_dict[comp].coords(axis='x')[0].shape[0])

    # Find the minimum resolution
    # The lowest resolution
    min_resn = min(n_lons_list)
    # Find out the component that corresponds to that lowest resolution
    min_comp = list(component_set)[n_lons_list.index(min_resn)]

    # If using observations in the equation then regrid to a coarse grid
    # (either N96 or the obs grid if less than N96)
    if using_obs:
        if min_resn < 192:
            regrid_grid = data_dict[min_comp]
        else:
            regrid_grid = globalvar.extra_data_dict['n96_land_frac']
            regrid_grid = get_cube_ready(regrid_grid)
            if globalvar.error:
                return 1, globalvar.missing_data
    else:

        # If no obs are present then use the coarsest resolution for the regrid
        regrid_grid = data_dict[min_comp]

    # Make an array of regridded data
    data_regrid_dict = {}
    for comp in component_set:
        if data_dict[comp] != regrid_grid:
            data_regrid_dict[comp] = regrid_cube(data_dict[comp], regrid_grid,
                                                 comp=comp)
        else:
            data_regrid_dict[comp] = data_dict[comp]
    if globalvar.error:
        return 1, globalvar.missing_data

    # Get a list of components
    component_list = list(component_set)
    n_components = len(component_list)

    # Check to see if any of the input fields have multiple levels.
    if same_levels:
        multi_levels = False
        vertical_level_name_list = ['model_level_number', 'pressure']
        for component in component_list:
            coords_list = [this_coord.name() for this_coord in
                           data_regrid_dict[component].coords()]
            for vertical_level_name in vertical_level_name_list:
                if vertical_level_name in coords_list:
                    multi_levels = True

        # If multipe levels then make sure the levels are the same
        if multi_levels:
            # Count up, comparing vertical grids and removing levels not
            #  the same
            if n_components > 1:
                for c in range(n_components-1):
                    (data_regrid_dict[component_list[c]],
                     data_regrid_dict[component_list[c+1]]) = \
                        ioc(data_regrid_dict[component_list[c]],
                            data_regrid_dict[component_list[c+1]])

            # Count down, comparing vertical grids and removing levels not
            #  the same
            if n_components > 2:
                for c in range(n_components-2, -1, -1):
                    (data_regrid_dict[component_list[c]],
                     data_regrid_dict[component_list[c+1]]) = \
                        ioc(data_regrid_dict[component_list[c]],
                            data_regrid_dict[component_list[c+1]])

    # Make sure they all have the same units
    same_units(data_regrid_dict)

    # Include data_regrid_dict in the equation
    long_eqn = expand_equation(equation, component_list)

    # Perform the equation
    try:
        toplot_cube = eval(long_eqn)
    except:
        print 'ERROR: failed to execute equation ' + long_eqn + \
            ' in perform_equation'
        print "ERROR: ", sys.exc_info()
        if globalvar.debug:
            pdb.set_trace()
        else:
            globalvar.error = True
            return 1, globalvar.missing_data

    if calc_rms:
        rms_float = rms.calc_all(rms_list, toplot_cube, key, page_title,
                                 filename=filename)
    else:
        rms_float = globalvar.missing_data

    return toplot_cube, rms_float


def quiver(u_cube, v_cube, scale=None):
    """
    Generate a vector plot with suitably trimmed down numbers of vectors.
    A wrapper for:
    http://matplotlib.org/1.3.1/api/pyplot_api.html?highlight=quiver#matplotlib.pyplot.quiver

    u_cube = A cube of x components of the vectors
    v_cube = A cube of y components of the vectors
    scale = The magnitude the arrow needs to be for it to be an inch long.
            Try to put in a value about 10 times your largest vector length.
    """

    # Determine the size of the new arrays
    len_x = min([len(u_cube.coords(axis='x')[0].points), 32])
    len_y = min([len(u_cube.coords(axis='y')[0].points), 24])

    # Determine the markup (the difference in the array dimensions]
    markup_x = len(u_cube.coords(axis='x')[0].points) / float(len_x)
    markup_y = len(u_cube.coords(axis='y')[0].points) / float(len_y)

    # Initialise some new arrays
    u_shrink = np.zeros((len_x, len_y))
    v_shrink = np.zeros((len_x, len_y))
    x_shrink = np.zeros(len_x)
    y_shrink = np.zeros(len_y)

    # Copy data across to the new arrays
    for i in range(len_x):
        i_ind = int(round(i*markup_x))
        for j in range(len_y):
            j_ind = int(round(j*markup_y))
            u_shrink[i, j] = u_cube.data[j_ind, i_ind]
            v_shrink[i, j] = v_cube.data[j_ind, i_ind]
            if i == 0:
                y_shrink[j] = u_cube.coords(axis='y')[0].points[j_ind]
        x_shrink[i] = u_cube.coords(axis='x')[0].points[i_ind]

    # Transpose the data as the y axis comes first in numpy arrays
    u_trans = np.transpose(u_shrink)
    v_trans = np.transpose(v_shrink)

    plt.quiver(x_shrink, y_shrink, u_trans, v_trans,
               scale_units='inches', scale=scale, pivot='tail')


def get_land_frac(incube):
    """
    Get the land fraction on the model grid

    incube = model cube (for the grid to use)
    """

    # This code only works if we have the land fraction cubes
    #  already been loaded into extra data
    if len(globalvar.extra_data_dict) == 0:
        print 'ERROR in get_land_frac in valmod.py'
        print '(possibly called from ocean_only in cubeinout.py)'
        print 'You need to run the following before calling these routines:'
        print 'import valmod as vm'
        print "extras_dict=vm.read_info_file(os.path.join(general_dir,"
        print "                                           'control',"
        print "                                           'extras_file.dat'))"
        print "vm.load_extra_data(extras_dict)"
        raise Exception

    # What is the correct land fraction field to use
    resolution = incube.coords(axis='x')[0].shape[0]/2.0
    n_lats = incube.coords(axis='y')[0].shape[0]
    if resolution < 72:
        land_frac = globalvar.extra_data_dict['n48_land_frac']
    if resolution >= 72 and resolution < 156:
        land_frac = globalvar.extra_data_dict['n96_land_frac']
    if resolution >= 156 and resolution < 384:
        land_frac = globalvar.extra_data_dict['n216_land_frac']
    if resolution >= 384:
        land_frac = globalvar.extra_data_dict['n512_land_frac']

    # Regrid it if needed
    if n_lats != land_frac.coords(axis='y')[0].shape[0]:
        print 'Regridding land fraction for '+merged_key
        land_frac = vm.regrid_cube(land_frac, incube)

    return land_frac
