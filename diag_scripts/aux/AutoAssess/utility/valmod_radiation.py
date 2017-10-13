#! /usr/local/sci/bin/python

### THIS IS AN OLDER VERSION OF VALMOD THAT IS REQUIRED BY THE RADIATION CODE ###

# This file contains extra routines needed for the validation note
# software but may also be useful elsewhere 

# Import files
import string, pdb, csv, iris, numpy, sys, os
import matplotlib.ticker
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
import globalvar
import numpy as np
import rms
import cf_units

# Import any routines in general/cma_python/cma_python.py for use in the equation
#from cma_python import *
import cubeinout

# **************** READ_CONTROL_FILE **********************

def read_control_file(controlfile):
    """Reads the control files and returns the data held in them in two nested dictionaries
    """

    # Initialise variables
    in_sec=0
    large_dict={}

    # Open the file
    f = open(controlfile, 'r')

    # Read each line in turn
    for line in f:

        # Split it up into substrings and remove whitespace
        line_split = string.split(line,':')
        line_split = [s.replace('\n','') for s in line_split]
        line_split_nw = [s.replace(' ','') for s in line_split]
   
        if len(line_split) > 1:
        
            # Ignore any lines beginning with hash
            if line_split_nw[0][0] != '#':
    
                # Work out if we are at the start or end of a section
                if line_split_nw[0] == 'start':
                    in_sec=1
                    section=line_split[1]
                    setting_dict={}
       
                elif line_split_nw[0] == 'end':
                    in_sec=0
                    # Add to the large dictionary
                    large_dict[section]=setting_dict
       
                else:
                    if in_sec:
      
                        # Add to the small dictionary
                        setting_dict[line_split_nw[0]]=line_split[1:]
       
    f.close()
    return large_dict

# **************** READ_INFO_FILE **********************
   
def read_info_file(infofile):
    """Reads the info file and returns the data held in it as a single dictionary
    """

    # Initialise variables
    info_dict={}
   
    # Read the valorder file as a csv file
    f = csv.reader(open(infofile,'r'), delimiter=':')
    for line in f:
        
        # Ignore any blank lines or those beginning with a hash
        if len(line) > 0:
            if line[0][0] != '#':
                info_dict[line[0]]=line[1]
      
    return info_dict

# **************** WRITE_INFO_FILE *****************

def write_info_file(infofile, info_dict):
    """Writes information to an info file.
    
    infofile=(string) a filename to hold this information. The final file will have
             two columns of data: keys on the left and items on the right, separated
             by a colon.
    info_dict=(dictionary) the information to store in the file
    """
    
    # Open a file ready for reading
    f = open(infofile, 'w')
    
    # Loop over the info_dict
    for info_key, info_item in info_dict.iteritems():
        f.write(info_key+':'+info_item+'\n')
        
    f.close()
    

# **************** KEY_VALID **********************
    
def key_valid(info_key):
    """Checks that info_key is a key that points to a file name"""
    
    acceptable_key_list=['djfm','mamm','jjam','sonm','annm','djfs','mams','jjas','sons','anns']
    return info_key in acceptable_key_list

# **************** LOAD_PICKLE **********************
    
def load_pickle(pickle_file):
    """Loads from a pickle file"""
    f=open(pickle_file)
    cube_list=pickle.load(f)
    f.close()
    
    return cube_list

# **************** SAVE_PICKLE **********************
    
def save_pickle(cube_list, pickle_file):
    """Saves to a pickle file"""
    f=open(pickle_file, 'wb')
    pickle.dump(cube_list, f)
    f.close()
   
    
# **************** SPLIT_EQUATION **********************
   
def split_equation(equation, n_underscores=0):
    """Reads in an equation and returns the component parts
    n_underscores = (int) number of underscores that the component needs to have in
                    order for it to be a component.
                    default = 0 (ignore this)
    """
      
    # Initialise arrays
    components=[]   
    this_comp=''
   
    # Loop over each charactor and check for maths symbols. Append the strings to the components array.
    for char in equation:
        if char == '+' or char == '-' or char == '*' or char == '/' or char == '^' or char == '(' or char == ')':
        
            # If you have a decent length component append it to a list of components
            if this_comp.replace(' ','') != '':
        
                # Components can't be purely numerical
                try:
                    this_comp_float=float(this_comp)
                except ValueError:
                    if n_underscores > 0:
                        if this_comp.count('_') == n_underscores:
                            components.append(this_comp.replace(' ',''))
                    else:
                        components.append(this_comp.replace(' ',''))
            this_comp=''
        else:
            this_comp=this_comp+char
            
    if this_comp.replace(' ','') != '':

        # Components can't be purely numerical
        try:
            this_comp_float=float(this_comp)
        except ValueError:    
            components.append(this_comp.replace(' ',''))   
   
    return components

# **************** APPLY_HEAVYSIDE **********************

def apply_heavyside(cube, heavyside_cube):
    """Uses the heavyside to mask out points where the field is below ground level
    
    cube = the input field
    heavyside_cube = the heavyside cube that contains 0 for data below ground, 1
                     for data above ground and fractions for partial situations.
    """
    
    mask = heavyside_cube.data <= 0.1
    cube.data = numpy.ma.array(cube.data, mask=mask)
    cube.data /= heavyside_cube.data
    return cube
    
# **************** EXPAND_EQUATION **********************

def expand_equation(equation, component_list,dict_name='data_regrid_dict'):
    """So far we have an equation but in order to execute it we need
    data_regrid_dict in front of each component in the equation. This routine
    adds those data_dict bits.
    
    equation = (string) the equation from the plots_file
    components = a list of components where each item is a string of the
                 form source_season_field
    dict_name = (string) the name of the dictionary that contains the
                 variables that you want in the equation. The default (as used
                 by the validation notes) is data_regrid_dict             
    """
    
    long_equation=equation
    
    # Loop over components
    for comp in component_list:
        long_equation=long_equation.replace(comp,dict_name+"['"+comp+"']",1)
        
    return long_equation

# **************** GET_CUBE_READY **********************

def get_cube_ready(cube):
    """Sets up cubes ready for regridding, in the equation and for generating means"""
    
    # Remove unnecessary coordinates in cube
    coords_list=['forecast_reference_time','forecast_period','source','season','time']     
    for coord in coords_list:
        try:
            cube.remove_coord(coord)
        except iris.exceptions.CoordinateNotFoundError:
            pass
            
    # Also make sure we are not using grid_longitude and grid_latitude
    fix_coords(cube)
    
    # Test to see if the grid does not have a coordinate system (CERES-EBAF). If so fix it.
    if cube.coord_system() == None:
        cube = fix_cube(cube)
        if globalvar.error: return 1
    
    # Test to see if the grid is a rotated coordinate system (HadSLP2). If so fix it.
    if cube.coord_system() == iris.coord_systems.RotatedGeogCS(0.0, 0.0, ellipsoid=iris.coord_systems.GeogCS(6371229.0)):
        cube = fix_cube(cube)
        if globalvar.error: return 1

    # Roll the cube
    coords_list = [this_coord.name() for this_coord in cube.coords()]
    if 'longitude' in coords_list or 'grid_longitude' in coords_list:
        cube = roll_cube(cube)
        
    # Guess the bounds
    if not cube.coord(axis='x').has_bounds(): cube.coord(axis='x').guess_bounds()
    if not cube.coord(axis='y').has_bounds(): cube.coord(axis='y').guess_bounds() 
    
    return cube

# **************** SAME_UNITS **************************

def same_units(cube_dict):
    """Make sure all the cubes in this dectionary have the same coordinates"""
    
    component_set=cube_dict.keys()
    
    first_unit=cube_dict[component_set[0]].units
    
    for comp in component_set: cube_dict[comp].units = first_unit
        

# **************** GET_COMPONENT_SET **********************

def get_component_set(equation, n_underscores=0):
    """Get a set of components used in this equation"""
    
    component_list=[]
    for comp in split_equation(equation, n_underscores=n_underscores):
        if comp != '': component_list.append(comp)
    component_set = set(component_list)
    
    return component_set

# ************* FIX_COORDS ************

def fix_coords(incube):
    """Makes so that all grid_longitude axes are renamed longitude
    and all grid_latitude axes are renamed latitude.
    
    incube = input cube to fix
    returned value = fixed cube
    """
    
    dims=[cor.standard_name for cor in incube.coords()]
    if 'grid_longitude' in dims: incube.coord('grid_longitude').standard_name = 'longitude'
    if 'grid_latitude' in dims: incube.coord('grid_latitude').standard_name = 'latitude'
    if u'longitude' in dims: incube.coord(u'longitude').standard_name = 'longitude'
    if u'latitude' in dims: incube.coord(u'latitude').standard_name = 'latitude'    
    
    
# **************** FIX_CUBE **********************

def fix_cube(incube):
    """Fixes this cube by putting it on a GeogCS coordinate system
    with standard latitudes and longitudes.
    
    incube = input cube to fix
    returned value = fixed cube
    """

    print 'Incompatible coordinate system detected. Fixing cube.'
    
    # Generate new latitudes and longitudes on a GeogCS coordinate system
    old_lats_dim = incube.coord(axis='y')
    new_lats_dim = iris.coords.DimCoord(old_lats_dim.points, standard_name='latitude', units=cf_units.Unit('degrees'), coord_system=iris.coord_systems.GeogCS(6371229.0))
    old_lons_dim = incube.coord(axis='x')
    new_lons_dim = iris.coords.DimCoord(old_lons_dim.points, standard_name='longitude', units=cf_units.Unit('degrees'), coord_system=iris.coord_systems.GeogCS(6371229.0))

    # Reshape the data to remove any extra, single length, dimensions.
    reshaped_data = np.reshape(incube.data, (new_lats_dim.shape[0], new_lons_dim.shape[0]))

    # Generate a new cube from the raw data and the new latitudes and longitudes
    outcube = iris.cube.Cube(reshaped_data, dim_coords_and_dims=[(new_lats_dim, 0), (new_lons_dim, 1)])
    
    # Make sure we have the correct coordinate system
    # outcube.coord(axis='x').coord_system = iris.coord_systems.GeogCS(6371229.0)
    # outcube.coord(axis='y').coord_system = iris.coord_systems.GeogCS(6371229.0)
    
    return outcube

    
# **************** REGRID_CUBE **********************

def regrid_cube(incube, newgrid, scheme_str='AreaWeighted'):
    """Regrids a cube to a new grid, retaining missing data:
    
    incube = input cube to regrid
    rewgrid = new grid cube
    scheme_str = (string) The scheme by which to do the regridding. One of:
       'AreaWeighted' = Use the area weighted scheme (default)
       'Linear' = Use a linear interpolation scheme
    """    
                
    # Assign the area weighted scheme for regridding
    if scheme_str == 'AreaWeighted': scheme = iris.analysis.AreaWeighted(mdtol=0.5)
    elif scheme_str == 'Linear': scheme = iris.analysis.Linear()
    else: print 'ERROR in valmod.regrid_cube: scheme_str value of ',scheme_str,' not supported'
    
    # Regrid the cube

    try:
        regrided_cube = incube.regrid(newgrid, scheme)
    except:
        print "ERROR: Cubes are not compatible for regridding." 
        print "ERROR: ", sys.exc_info()
        if globalvar.debug:
            pdb.set_trace()     
        else:
            globalvar.error=True
            return 1

    return regrided_cube
    
# **************** ROLL_CUBE *************************

def roll_cube(cube, start_longitude=-180.0):
    """Rolls the cube so that it starts from 180W and goes to 180E (instead of from 0E to 360E).
    
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
    degrees_to_shift = 360.0 * number_boxes_to_shift/float(cube.shape[long_axis])
    
    # Roll the data
    cube.data = numpy.roll(cube.data, number_boxes_to_shift, axis=long_axis)
    
    # Roll the coordinates
    cube.coord(axis='x').points = cube.coord(axis='x').points+degrees_to_shift
    
    # Roll the bounds
    if cube.coord(axis='x').has_bounds(): cube.coord(axis='x').bounds = cube.coord(axis='x').bounds+degrees_to_shift
    
    return cube

# **************** AREA_AVG **************************

def area_avg(cube):
    """Perform an area average of a cube using weights to account for
    changes in latitude."""
    
    # What dimension is latitude
    thiscoords = cube.coords()
    lat_dim=-99
    for i in range(len(thiscoords)):
        if thiscoords[i].standard_name == 'latitude': lat_dim=i
    if lat_dim == -99: raise ValError('Could not find latitude coordinate in cube')
    
    # Convert latitudes to radians
    lats_radians_array = cube.coord(standard_name='latitude').points/360.0 * 2.0 * numpy.pi

    # Generate the cos of these latitudes
    cos_lats_array = numpy.cos(lats_radians_array)

    # Make a 2d array of these cos values
    cos_2d_array = numpy.zeros(cube.data.shape)
    #cos_2d_array = numpy.empty_like(cube.data)
    for i in range(cos_2d_array.shape[0]):
        for j in range(cos_2d_array.shape[1]):
            got_data=True
            location=[i,j]
            lat_index=location[lat_dim]
            if hasattr(cube.data, 'mask'):
                if cube.data.mask.size == cube.data.size:
                    got_data = not cube.data.mask[i,j]
            if got_data: cos_2d_array[i,j] = cos_lats_array[lat_index]

    # Do some area averaging
    area_average = numpy.average(cube.data, weights=cos_2d_array)
    
    return area_average
   
# **************** ADD_COLORBAR **********************
    
def add_colorbar(mappable,units=None):
    # adds a colorbar to a plot with some standard keywords set
    cbar_kwargs = {'pad': 0.1, 'orientation': 'horizontal', 'extend': 'both'}
    
    cax, cbar_kwargs_filter = matplotlib.colorbar.make_axes(plt.gca(), **cbar_kwargs)
    cb = plt.colorbar(mappable, cax=cax, **cbar_kwargs_filter)
    cb.ax.tick_params(length=0)
    
    # Add the units
    if not globalvar.pub and 'units' in locals(): cb.set_label(units)

    return cb
    
# **************** LatsFormatter **********************

class LatsFormatter(matplotlib.ticker.Formatter):
    """
    This is a custom formatter that converts the native unit of degrees into North/South degrees.
    """
    def __call__(self, x, pos=None):
        letter = 'S' if x < 0 else 'N' 
        return u"%d%s" % (abs(x), letter)
        
# **************** LonsFormatter **********************

class LonsFormatter(matplotlib.ticker.Formatter):
    """
    This is a custom formatter that converts the native unit of degrees into West/East degrees.
    """
    def __call__(self, x, pos=None):
        letter = 'W' if x < 0 else 'E' 
        return u"%d%s" % (abs(x), letter)
        
# **************** ValError ************************
                    
class ValError(Exception):
    """
    This is the default exception for errors picked up by the validation note software
    """
    def __init__(self, value):
        self.value = value
        print value
    def __str__(self):
        return repr(self.value)
        
# **************** MAKE_CMAP ***************************

def make_cmap(red_list, green_list, blue_list):
    """
    This will make a colour map given arrays of red, green and blue colours
    """
    
    # Check that the lengths of the red, green and blue arrays are the same
    if len(red_list) != len(green_list) or len(red_list) != len(blue_list): raise ValError('The lengths of the red, green and blue colour tables are not equal')
    n_col = len(red_list)
    
    # Make a blank colour dictionary
    cdict={}
    
    # Loop over the colours
    for colour in ['red','green','blue']:
    
        # Set up a generic list for each of the colours
        colour_list=eval(colour+'_list')
        
        # Check that there are no values greater than one
        if max(colour_list) > 1.0: raise ValError('There is a colour greater than 1.0 in the '+colour+' colour list')
        
        # Make a blank colour list just for this colour
        one_cmap_list=[]
        
        for c in range(n_col):
        
            # What is the array index (from 0 to 1)
            array_index=c/float(n_col-1)
        
            # Populate the one colour list
            one_cmap_list.append((array_index,colour_list[c],colour_list[c]))
            
        # Add this to the colour dictionary
        cdict[colour]=one_cmap_list
            
    # Make the colour map 
    cmap=matplotlib.colors.LinearSegmentedColormap('default', cdict)
    
    return cmap

# ************** FILE_NAME ************************

def file_name(page_num,page_title):
    """ Make a file name derived from the page number and page title.
    This is done by replacing all spaces and in the page title with underscores and
    removing all brackets."""
    
    # Replace spaces with underscores in page_title
    page_title = page_title.replace(' ','_')
    
    # Remove all brackets in page_title
    page_title = page_title.replace('(','')
    page_title = page_title.replace(')','')
    
    # Join page_num with page_title
    if globalvar.pub:
        filename=page_num+'_'+page_title+'.eps'
    else:
        filename=page_num+'_'+page_title+'.png'
    return filename

# ************* MK_TEMP ******************

def mk_temp():
    """ Make a temporary directory on your local disk to store pkl files and
    page titles. 
    
    returned_value = full path to temporary directory.
    """

    localdata=os.getenv('LOCALDATA')
    temp_dir=os.path.join(localdata,'valnote_temp')
    if not os.path.exists(temp_dir): os.mkdir(temp_dir)
    return(temp_dir)

# ************* EXTRACT_NO_DIURNAL ******************

def extract_no_diurnal(all_data_cube, all_constraints=None):
    """ Extracts one cube from all_data_cube based on the constraints provided.
    It will remove any diurnal cycle diagnostics.
    
    all_data_cube = a cube list containing all the data
    all_constraints = a list of constraints needed to extract the data you want. This routine will
                      attempt to extract data from each constraint sequentially until data is
                      successfully returned.
    """

    # Extract the data
    if all_constraints:
    
        # Loop over each constraint      
        for this_constraint in all_constraints:
        
            # Attempt to extract the data
            print 'Extracting using constraint: '
            print this_constraint
            
            this_data_cube=all_data_cube.extract(this_constraint)
        
            # Test to see if you managed to get some data. If so exit this for loop.
            cube_length=1
            if isinstance(this_data_cube, iris.cube.CubeList): cube_length=len(this_data_cube)
            if cube_length != 0:
                break
              
    else:
        this_data_cube=all_data_cube
    
    # If this is a cube list then get it down to just one cube
    if isinstance(this_data_cube, iris.cube.CubeList):
        cube_length=len(this_data_cube)
        
        # Check to make sure we have some data
        if cube_length == 0:
	    print 'ERROR: No cubes returned from all_data_cube.extract command in extract_no_diurnal'
	    if globalvar.debug:
	        pdb.set_trace()     
	    else:
	        globalvar.error=True
	        return 1
 
        # If there is just one cube in this list extract that
        if cube_length == 1:
            this_data_cube = this_data_cube[0]
            
        if cube_length > 1:
            print 'WARNING: detecting cube list in extract_no_diurnal. Using the first cube (but this may be wrong).'
            this_data_cube = this_data_cube[0]
           
    # Get a list of dimensions
    dims_name_list=[]
    for this_dim in this_data_cube.coords():
        if type(this_dim.standard_name) in (str, unicode):
            dims_name_list.append(this_dim.standard_name)
        elif type(this_dim.long_name) in (str, unicode):
            dims_name_list.append(this_dim.long_name)
        else:
            print 'ERROR: No standard name or long name detected for coordinate in extract_no_diurnal:'
            print this_dim
            if globalvar.debug:
                pdb.set_trace()      
            else:
                globalvar.error=True
                return 1        
    
    # Get the longitude, latitude and vertical dimension names
    lat_name='x'
    for dim in dims_name_list:
        if 'lat' in dim: lat_name=dim
    lon_name='x'
    for dim in dims_name_list:
        if 'lon' in dim: lon_name=dim
    vert_name='x'
    for dim in dims_name_list:
        if 'pressure' in dim or 'vertical' in dim or 'height' in dim or 'level' in dim: 
            vert_name_temp=dim
            if this_data_cube.coords(vert_name_temp)[0].shape[0] > 1: vert_name=vert_name_temp
    
    # Make sure we have found longitudes and latitudes   
    if lat_name == 'x' or lon_name == 'x':
        print 'ERROR: could not detect either longitude or latitude axes'
        if globalvar.debug:
            pdb.set_trace()      
        else:
            globalvar.error=True
            return 1

    # Initialise some variables for future use
    greatest_bounds_zero = 0.0
    smallest_bounds_two = 1.0e9
    
    # Assess whether this a two or three dimensional and assign a list of dimensions accordingly
    if vert_name != 'x':
       this_dims_list=[vert_name, lat_name, lon_name]
    else:
       this_dims_list=[lat_name, lon_name]
    
    # Make a longitude/latitude slice and assign it as the output
    for subset_cube in this_data_cube.slices(this_dims_list): 
        if 'output_cube' not in locals(): output_cube=subset_cube
        
        # For model output data
        if len(subset_cube.aux_coords) == 3:
            if subset_cube.aux_coords[0].bounds[0][0] >= greatest_bounds_zero and subset_cube.aux_coords[2].bounds[0][0] <= smallest_bounds_two:
                output_cube=subset_cube
                greatest_bounds_zero = subset_cube.aux_coords[0].bounds[0][0]
                smallest_bounds_two = subset_cube.aux_coords[2].bounds[0][0]
    
    # If the output cube doesn't exist then stop with an error
    if 'output_cube' not in locals():
        print 'ERROR: could not detect time mean diagnostic'
        if globalvar.debug:
            pdb.set_trace()      
        else:
            globalvar.error=True
            return 1        
    
    return output_cube

# **************** LOAD_EXTRA_DATA **********************
    
def load_extra_data(extras_dict):
    """Loads in extra data, such as land-sea masks, into a dictionary"""
    
    # Initialise variable
    globalvar.extra_data_dict={}
    print 'The line of globalvar.extra_data_dict={}', globalvar.extra_data_dict
    for key in extras_dict.iterkeys():

        globalvar.extra_data_dict[key] = iris.load_cube(extras_dict[key])
        globalvar.extra_data_dict[key] = get_cube_ready(globalvar.extra_data_dict[key])
# **************** PERFORM_EQUATION **********************

def perform_equation(equation,data_dict,page_title,key,rms_list=None,filename=None,calc_rms=False,n_underscores=0):
    """Perform the equations provided in plots_file.dat
    
    equation = the equation to perform
    data_dict = a dictionary containing the data needed for those equations
    page_title = the title of the page from valorder.dat
    key = the key (a, b, c or d) for the plot
    rms_list = an rms object for storing root mean square values of the final field (see rms.py for object structure)
    calc_rms = switch to turn on calculation of rms values
    filename = (string) the filename of the png plot
    n_underscores = (int) number of underscores that the component needs to have in
                    order for it to be a component.
                    default = 0 (ignore this)
    """
    
    print 'Performing equation '+equation
           
    # Generate a list of components used in the equation
    component_set=get_component_set(equation, n_underscores=n_underscores)
    
    # Are any of these components observations.
    using_obs=0
    for comp in component_set:
        input_info = string.split(comp,'_')
        if input_info[0] != 'exper' and input_info[0] != 'control':
            using_obs=1
            
        # Check that the component exists in data_dict
        if comp not in data_dict.keys():
            print 'ERROR: Not all variables in equation '+equation+' have entries in data_dict. data_dict.keys()='
            print data_dict.keys()
            if globalvar.debug:
                pdb.set_trace()
            else:
                globalvar.error=True
                return 1, -999

    # Make sure all the fields have correct units and axes
    for comp in component_set:
        data_dict[comp] = get_cube_ready(data_dict[comp])
        if globalvar.error: return 1, -999
    same_units(data_dict)    

    # Get a list of longitude sizes
    n_lons_list=[]
    for comp in component_set:
        n_lons_list.append(data_dict[comp].coords(axis='x')[0].shape[0])
        #print 'comp', comp, n_lons_list, data_dict[comp].coords(axis='x')[0].shape
    # Find the minimum resolution
    min_resn = min(n_lons_list)                                             # The lowest resolution
    min_comp = list(component_set)[n_lons_list.index(min_resn)]             # Find out the component that corresponds to that lowest resolution    
    
    # If using observations in the equation then regrid to a coarse grid (either N96 or the obs grid if less than N96) 
    #if using_obs:    
    #    if min_resn < 192:
    #        regrid_grid = data_dict[min_comp]
    #    else:
    #        regrid_grid = globalvar.extra_data_dict['n96_land_frac']
    #        regrid_grid = get_cube_ready(regrid_grid) 
    #        if globalvar.error: return 1, -999   
    #else:
        
        # If no obs are present then use the coarsest resolution for the regrid
    #    regrid_grid = data_dict[min_comp]
    
    regrid_grid = globalvar.extra_data_dict['1deg_1deg_landsea_frac']
    
    # Make an array of regridded data
    data_regrid_dict={}
    for comp in component_set:
        if data_dict[comp] != regrid_grid:
            data_regrid_dict[comp] = regrid_cube(data_dict[comp],regrid_grid)
        else:
            data_regrid_dict[comp] = data_dict[comp]
    if globalvar.error: return 1, -999  
        
    # Get a list of components
    component_list = list(component_set)
    n_components = len(component_list)
    
    # If this is a zonal mean then make sure all your fields have the same levels
    if globalvar.plot_type ==  'zonal_mean':
        
        # Count up, comparing vertical grids and removing levels not the same
        if n_components > 1:
            for c in range(n_components-1):
                data_regrid_dict[component_list[c]], data_regrid_dict[component_list[c+1]]=iris.analysis.maths.intersection_of_cubes(data_regrid_dict[component_list[c]], data_regrid_dict[component_list[c+1]])
            
        # Count down, comparing vertical grids and removing levels not the same
        if n_components > 2:
            for c in range(n_components-2,0,-1):
                data_regrid_dict[component_list[c]], data_regrid_dict[component_list[c+1]]=iris.analysis.maths.intersection_of_cubes(data_regrid_dict[component_list[c]], data_regrid_dict[component_list[c+1]])         
    
    # Include data_regrid_dict in the equation
    long_equation = expand_equation(equation, component_list) 
    
    # Perform pre equation processing. e.g. zonal meaning or vertical meaning. This is specified by
    # the proc input in item_file.dat  
    if globalvar.valnote:
        for comp in component_set:
            input_info = string.split(comp,'_')
            if globalvar.item_dict[input_info[2]].has_key('proc'):            
                if globalvar.item_dict[input_info[2]]['proc'][0] == 'zonal_mean':
                    data_regrid_dict[comp] = data_regrid_dict[comp].collapsed('longitude', iris.analysis.MEAN)
                if globalvar.item_dict[input_info[2]]['proc'][0] == 'vertical_mean': 
                    data_regrid_dict[comp] = data_regrid_dict[comp].collapsed('pressure', iris.analysis.MEAN)   

    # Perform the equation
    try:
        toplot_cube=eval(long_equation)
    except:
        print 'ERROR: failed to execute equation '+long_equation+' in perform_equation'
        print "ERROR: ", sys.exc_info()
        if globalvar.debug:
            pdb.set_trace()      
        else:
            globalvar.error=True
            return 1, -999   
    
    # Perform post equation processing. e.g. zonal meaning or vertical meaning. This is specified by
    # the proc input in levels_file.dat
    if globalvar.valnote:
        level_type = globalvar.plots_dict[page_title][key][1]
        coords_list = [this_coord.name() for this_coord in toplot_cube.coords()]
        if globalvar.levels_dict[level_type].has_key('proc'):
            if globalvar.levels_dict[level_type]['proc'][0] == 'zonal_mean' and 'longitude' in coords_list:
                toplot_cube = toplot_cube.collapsed('longitude', iris.analysis.MEAN)
            if globalvar.levels_dict[level_type]['proc'][0] == 'vertical_mean':
                toplot_cube = toplot_cube.collapsed('pressure', iris.analysis.MEAN)
                
    if calc_rms:
        rms_float = rms.calc_all(rms_list, toplot_cube, key, page_title, filename=filename)  
    else:
        rms_float = -999                            
    
    return toplot_cube, rms_float
