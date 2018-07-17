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
import numpy as np
import cf_units



def get_cube_ready(cube):
    """Sets up cubes ready for regridding, in the equation and for generating means"""
    
    # Remove unnecessary coordinates in cube
    to_remove_list=['forecast_reference_time',
                    'forecast_period','source','season','time']     
    for coord in cube.coords():
        if coord.name() in to_remove_list:
            cube.remove_coord(coord)
        
    # Guess the bounds
    if not cube.coord(axis='x').has_bounds():
        cube.coord(axis='x').guess_bounds()
    if not cube.coord(axis='y').has_bounds():
        cube.coord(axis='y').guess_bounds() 
    
    return cube
 

def area_avg(cube, coord1=None, coord2=None):
    """Perform an area average of a cube using weights to account for
    changes in latitude."""
    import iris.analysis.cartography
    for coord in (coord1, coord2):
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(cube)
    result = cube.collapsed(
        [coord1, coord2], iris.analysis.MEAN, weights=grid_areas)
    return result


def add_colorbar(mappable,units=None):
    # adds a colorbar to a plot with some standard keywords set
    cbar_kwargs = {'pad': 0.1, 'orientation': 'horizontal', 'extend': 'both'}
    
    cax, cbar_kwargs_filter = matplotlib.colorbar.make_axes(plt.gca(), **cbar_kwargs)
    cb = plt.colorbar(mappable, cax=cax, **cbar_kwargs_filter)
    cb.ax.tick_params(length=0)
    
    # Add the units
    if not globalvar.pub and 'units' in locals(): cb.set_label(units)

    return cb


class LatsFormatter(matplotlib.ticker.Formatter):
    """
    This is a custom formatter that converts the native unit of degrees into North/South degrees.
    """
    def __call__(self, x, pos=None):
        letter = 'S' if x < 0 else 'N' 
        return u"%d%s" % (abs(x), letter)


class LonsFormatter(matplotlib.ticker.Formatter):
    """
    This is a custom formatter that converts the native unit of degrees into West/East degrees.
    """
    def __call__(self, x, pos=None):
        letter = 'W' if x < 0 else 'E' 
        return u"%d%s" % (abs(x), letter)

                    
class ValError(Exception):
    """
    This is the default exception for errors picked up by the validation note software
    """
    def __init__(self, value):
        self.value = value
        print(value)
    def __str__(self):
        return repr(self.value)


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
            print('Extracting using constraint: ')
            print(this_constraint)
            
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
            print('ERROR: No cubes returned from all_data_cube.extract command in extract_no_diurnal')
            if globalvar.debug:
                pdb.set_trace()     
            else:
                globalvar.error=True
                return 1
 
        # If there is just one cube in this list extract that
        if cube_length == 1:
            this_data_cube = this_data_cube[0]
            
        if cube_length > 1:
            print('WARNING: detecting cube list in extract_no_diurnal. Using the first cube (but this may be wrong).')
            this_data_cube = this_data_cube[0]
           
    # Get a list of dimensions
    dims_name_list=[]
    for this_dim in this_data_cube.coords():
        if type(this_dim.standard_name) in (str, unicode):
            dims_name_list.append(this_dim.standard_name)
        elif type(this_dim.long_name) in (str, unicode):
            dims_name_list.append(this_dim.long_name)
        else:
            print('ERROR: No standard name or long name detected for coordinate in extract_no_diurnal:')
            print(this_dim)
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
        print('ERROR: could not detect either longitude or latitude axes')
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
        print('ERROR: could not detect time mean diagnostic')
        if globalvar.debug:
            pdb.set_trace()      
        else:
            globalvar.error=True
            return 1        
    
    return output_cube


def perform_equation(dataset_1, dataset_2, analysis_type):
    """
    Perform simple cube subtraction
    
    analysis_type = type of analysis (zonal_mean, vertical_mean,...)
    """
    # Make sure all the fields have correct units
    dataset_1_ready = get_cube_ready(dataset_1)
    dataset_2_ready = get_cube_ready(dataset_2)

    if  analysis_type == 'zonal_mean':
        dataset_1_mean = dataset_1_ready.collapsed('longitude',
                                                   iris.analysis.MEAN)
        dataset_2_mean = dataset_2_ready.collapsed('longitude',
                                                   iris.analysis.MEAN)

    elif analysis_type == 'vertical_mean':
        dataset_1_mean = dataset_1_ready.collapsed('pressure',
                                                   iris.analysis.MEAN)
        dataset_2_mean = dataset_2_ready.collapsed('pressure',
                                                   iris.analysis.MEAN)
    elif analysis_type == 'lat_lon':
        dataset_1_mean = dataset_1_ready
        dataset_2_mean = dataset_2_ready

    # Perform simple difference
    toplot_cube = dataset_1_mean - dataset_2_mean
    
    return toplot_cube
