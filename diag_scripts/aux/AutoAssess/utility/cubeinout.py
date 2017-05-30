'''
File that contains code that modifies a cube in some form or another and
returns another (modified) cube.

This code is used by the validation notes and the validation notes expect
the output to be a cube. The dimensions of the output cube don't have to
be the same as the dimensions of the input cube so therefore operations like
zonal meaning are allowed here.
'''

import math

import numpy as np

import iris

import valmod as vm

# Set up global constants
global radius_earth, gravity_earth
radius_earth = 6.37123e+06
gravity_earth = 9.80665

# Set up dictionaries
land_frac_dict = {}


def zy_stream(vfield):
    '''
    Function to calculate the vertical integral of vfield from the surface to the
    top of the atmosphere and so calculate the meridional streamfunction. (vfield
    is the v component of the wind zonally averaged)

    As in pp package vstrm routine, the streamfunction is given by:-

              /
          S= \ 2.pi.a.cos(lat).v dp
              /

    vfield - input vfield must be an iris cube of v winds varying only in the
    latitude and pressure dimensions.
    '''

    # Set up variables needed for looping over pressure coordinates
    pressure_coord = vfield.coord('pressure').points.tolist()
    pressure_coord.reverse()
    pressure_coord_len = len(pressure_coord)

    # Set up latitude variables
    latitude_coord_arr = vfield.coord('latitude').points
    latitude_coord_len = len(latitude_coord_arr)

    # Generate numpy arrays for holding the integrals and streamfunction data
    integral = np.zeros((latitude_coord_len))
    strm_arr = np.zeros((pressure_coord_len, latitude_coord_len))

    # Loop over pressure levels from bottom level (1000 hPa) to top
    for p_index in range(pressure_coord_len):
        if p_index > 0:

            # Work out the pressure difference between this and the last level
            dp = pressure_coord[p_index] - pressure_coord[p_index-1]

            # Integral of mean velocity in layer wrt pressure in mb, hence x50 ( x100 x0.5)
            vfield_this = vfield.extract(iris.Constraint(pressure=pressure_coord[p_index]))
            vfield_last = vfield.extract(iris.Constraint(pressure=pressure_coord[p_index-1]))
            integral += (vfield_this.data+vfield_last.data)*50.0*dp

            # Set values for factor to multiply v by in integral
            cos_latitude = np.cos(latitude_coord_arr*math.pi/180.0)
            factor = 2.0*math.pi*radius_earth*cos_latitude/gravity_earth

            # Add this to the final stream function array.
            # Note use reverse indexing as the order in the pressure dimension is reversed
            strm_arr[-1*p_index-1, :] = integral * factor

    # Make a cube from the resulting streamfunction
    strm_cube = iris.cube.Cube(strm_arr, aux_coords_and_dims=[(vfield.coord('pressure'), 0),
                                                              (vfield.coord('latitude'), 1)])

    return strm_cube


def log10(incube):
    '''
    Function to return the log (to base 10) of an input cube.
    '''

    # Do the maths
    log_array = np.log10(incube.data)

    # Put that data into a new cube
    outcube = incube.copy()
    outcube.data = log_array

    return outcube


def zonal_mean(incube):
    '''
    Function to generate a zonal mean of a cube.
    '''

    coords_list = [this_coord.name() for this_coord in incube.coords()]

    # If surface altitude is in the coordinates then remove it
    if 'surface_altitude' in coords_list:
        incube.remove_coord('surface_altitude')

    # Do the zonal mean
    if 'longitude' in coords_list:
        outcube = incube.collapsed('longitude', iris.analysis.MEAN)
    else:
        print 'ERROR (cubeinout.zonal_mean): longitude not a coordinate'

    return outcube


def meridional_mean(incube):
    '''
    Function to generate a meridional mean of a cube.
    '''

    coords_list = [this_coord.name() for this_coord in incube.coords()]
    if 'latitude' in coords_list:
        grid_areas = iris.analysis.cartography.area_weights(incube)
        outcube = incube.collapsed('latitude', iris.analysis.MEAN,
                                   weights=grid_areas)
    else:
        print 'ERROR (cubeinout.meridional_mean): latitude not a coordinate'

    return outcube


def vertical_mean(incube, zcoord_name=None):
    '''
    Function to generate a vertical mean of a cube.
    '''

    # Determine what your cube's vertical coordinate is
    if zcoord_name is None:
        coords_list = [this_coord.name() for this_coord in incube.coords()]
        if 'pressure' in coords_list:
            zcoord_name = 'pressure'
        elif 'model_level_number' in coords_list:
            zcoord_name = 'model_level_number'
        else:
            print "ERROR: vertical_mean can't determine what vertical " + \
                "coordinate your cube has"

    # Do the vertical meaning
    outcube = incube.collapsed(zcoord_name, iris.analysis.MEAN)

    return outcube


def ocean_only(incube):
    '''
    Function to mask out the land points
    '''

    # Get the land fraction field on the correct grid
    land_frac = vm.get_land_frac(incube)

    # Make and apply a mask
    mask = land_frac.data >= 0.5
    incube.data = np.ma.array(incube.data, mask=mask)
    return incube


def process_isccp(isccp_hist, isccp_weights, cloud_type):
    '''
    Function to process the isccp data and generate a cloud amount for a particular cloud type.
    isccp_hist = 4D cube with dimenstions: pseudo_level, pressure, latitude, longitude
                 The ISCCP histograms
    isccp_weights = (cube) Weighting function per lat-lon grid point
    cloud_type = (integer) one of:
                 0 = 'Optically thin clouds too thin to be seen by ISCCP'
                 1 = 'Low-top Thin Cloud'
                 2 = 'Low-top Medium Cloud'
                 3 = 'Low-top Thick Cloud'
                 4 = 'Mid-level-top Thin Cloud'
                 5 = 'Mid-level-top Medium Cloud'
                 6 = 'Mid-level-top Thick Cloud'
                 7 = 'High-top Thin Cloud'
                 8 = 'High-top Medium Cloud'
                 9 = 'High-top Thick Cloud'
    '''

    # Define the matrix. Each point on the matrix has a number which refers
    # to the cloud_type number.
    d2_type_matrix = np.array([[0, 7, 7, 8, 8, 9, 9],
                               [0, 7, 7, 8, 8, 9, 9],
                               [0, 7, 7, 8, 8, 9, 9],
                               [0, 4, 4, 5, 5, 6, 6],
                               [0, 4, 4, 5, 5, 6, 6],
                               [0, 1, 1, 2, 2, 3, 3],
                               [0, 1, 1, 2, 2, 3, 3]])

    # Generate a list of cubes to add together
    isccp_cubelist = iris.cube.CubeList()
    for i in range(d2_type_matrix.shape[0]):
        for j in range(d2_type_matrix.shape[1]):
            if d2_type_matrix[j, i] == cloud_type:
                isccp_cubelist.append(isccp_hist[i, j, :, :] / isccp_weights)

    # Add them together
    if len(isccp_cubelist) > 0:
        isccp_merge = isccp_cubelist.merge()[0]
        isccp_sub_total = isccp_merge.collapsed('pseudo_level',
                                                iris.analysis.SUM)
        isccp_total = isccp_sub_total.collapsed('pressure', iris.analysis.SUM)
    else:
        print 'ERROR could not find any cloud_type in the d2_type_matrix ' + \
            'in process_isccp (cubeinout.py)'

    return isccp_total


def anomaly_wrt_zonal_mean(incube):
    '''
    Function to calculate departures from the zonal mean = field - [field]
    '''

    zonal_mean_cube = zonal_mean(incube)
    anomaly_wrt_zonal_mean = incube - zonal_mean_cube
    return anomaly_wrt_zonal_mean
