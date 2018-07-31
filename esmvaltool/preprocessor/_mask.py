"""
_mask.py

module that performs missing values masking
and geographical area eslection
"""

from __future__ import print_function

import os
import logging

import iris
import numpy as np
from iris.analysis import Aggregator
from iris.util import rolling_window

logger = logging.getLogger(__name__)


def _check_dims(cube, mask_cube):
    """Check for same dims for mask and data"""
    x_dim = cube.coord('longitude').points.ndim
    y_dim = cube.coord('latitude').points.ndim
    mx_dim = mask_cube.coord('longitude').points.ndim
    my_dim = mask_cube.coord('latitude').points.ndim
    len_x = len(cube.coord('longitude').points)
    len_y = len(cube.coord('latitude').points)
    len_mx = len(mask_cube.coord('longitude').points)
    len_my = len(mask_cube.coord('latitude').points)
    if x_dim == mx_dim and y_dim == my_dim:
        if len_x == len_mx and len_y == len_my:
            logger.debug('Data cube and fx mask have same dims')
            return True
        else:
            logger.error('Data cube and fx mask have different grids!')
    else:
        logger.error('Data cube and fx mask differ in dims\n')
        logger.error('x=(%i, %i), y=(%i, %i)', x_dim, mx_dim, y_dim, my_dim)


def _get_fx_mask(fx_data, fx_option, mask_type):
    """Build a 50 percent land or sea mask"""
    inmask = np.zeros_like(fx_data, bool)
    if mask_type == 'sftlf':
        if fx_option == 'land':
            # Mask land out
            inmask[fx_data > 50.] = True
        elif fx_option == 'sea':
            # Mask sea out
            inmask[fx_data <= 50.] = True
    elif mask_type == 'sftof':
        if fx_option == 'land':
            # Mask land out
            inmask[fx_data < 50.] = True
        elif fx_option == 'sea':
            # Mask sea out
            inmask[fx_data >= 50.] = True

    return inmask


def _apply_fx_mask(fx_mask, var_data):
    """Apply the fx mask"""
    # Broadcast mask
    var_mask = np.zeros_like(var_data, bool)
    var_mask = np.broadcast_to(fx_mask, var_mask.shape).copy()

    # Aplly mask accross
    if np.ma.is_masked(var_data):
        var_mask |= var_data.mask

    # Build the new masked data
    var_data = np.ma.array(var_data, mask=var_mask, fill_value=1e+20)

    return var_data


def mask_landsea(cube, fx_file, mask_out):
    """Apply a land/sea mask"""
    # mask_out: is either 'land' or 'sea'
    # Dict to store the Natural Earth masks
    cwd = os.path.dirname(__file__)
    # ne_10m_land is fast; ne_10m_ocean is very slow
    shapefiles = {
        'land': os.path.join(cwd, 'ne_masks/ne_10m_land.shp'),
        'sea': os.path.join(cwd, 'ne_masks/ne_50m_ocean.shp')
    }

    if fx_file:
        # Try loading; some fx files may be broken (example bcc)
        try:
            fx_cube = iris.load_cube(fx_file)
            # preserve importance order: stflf then sftof
            if os.path.basename(fx_file).split('_')[0] == 'sftlf':
                if _check_dims(cube, fx_cube):
                    landsea_mask = _get_fx_mask(fx_cube.data, mask_out,
                                                'sftlf')
                    cube.data = _apply_fx_mask(landsea_mask, cube.data)
            elif os.path.basename(fx_file).split('_')[0] == 'sftof':
                if _check_dims(cube, fx_cube):
                    landsea_mask = _get_fx_mask(fx_cube.data, mask_out,
                                                'sftof')
                    cube.data = _apply_fx_mask(landsea_mask, cube.data)
            else:
                logger.warning('Masking with %s file', shapefiles[mask_out])
                cube = _mask_with_shp(cube, shapefiles[mask_out])
        except iris.exceptions.TranslationError as msg:
            logger.warning('Could not load fx file !')
            logger.warning(msg)
            logger.warning('Will use Natural Earth mask instead')
            cube = _mask_with_shp(cube, shapefiles[mask_out])
    else:
        # Mask with Natural Earth (NE) files
        logger.warning('Masking with %s file', shapefiles[mask_out])
        cube = _mask_with_shp(cube, shapefiles[mask_out])

    return cube


def masked_cube_simple(mycube, slicevar, v_min, v_max, threshold):
    """
    Mask function 1 -- simple cube cropping

    masking for a specific variable slicevar (string)
    arguments: cube, variable, min value, max value, threshold

    """
    import numpy.ma as ma
    coord_names = [coord.name() for coord in mycube.coords()]
    if slicevar in coord_names:
        coord = mycube.coord(slicevar)
        print('Masking on variable: %s' % coord.standard_name)
        cubeslice = mycube.extract(
            iris.Constraint(coord_values={
                coord.standard_name:
                lambda cell: v_min <= cell.point <= v_max
            }))
        if cubeslice is not None:
            masked_cubeslice = cubeslice.copy()
            masked_cubeslice.data = ma.masked_greater(cubeslice.data,
                                                      threshold)
            return masked_cubeslice
        else:
            logger.info('NOT masking the cube')
            return mycube
    else:
        logger.info('Var is not a cube dimension, leaving cube untouched')
        return mycube


def masked_cube_lonlat(mycube, lonlat_list, threshold):
    """
    Mask function 2 -- simple cube cropping on (min,max) lon,lat

    Builds a box and keeps only the values inside the box
    args: cube, min value, max value, where value=(lon, lat), threshold

    """
    import numpy.ma as ma
    lon1, lon2, lat1, lat2 = lonlat_list
    cubeslice = mycube.extract(
        iris.Constraint(
            longitude=lambda v: lon1 <= v.point <= lon2,
            latitude=lambda v: lat1 <= v.point <= lat2))
    if cubeslice is not None:
        masked_cubeslice = cubeslice.copy()
        masked_cubeslice.data = ma.masked_greater(cubeslice.data, threshold)
        print('Masking cube on lon-lat')
        return masked_cubeslice
    else:
        print('NOT masking the cube')
        return mycube


def _get_geometry_from_shp(shapefilename):
    """Get the mask geometry out from a shapefile"""
    import cartopy.io.shapereader as shpreader
    reader = shpreader.Reader(shapefilename)
    # Index 0 grabs the lowest resolution mask (no zoom)
    main_geom = [contour for contour in reader.geometries()][0]
    return main_geom


def _mask_with_shp(cube, shapefilename):
    """Apply a Natural Earth land/sea mask"""
    import shapely.vectorized as shp_vect

    # Create the region
    region = _get_geometry_from_shp(shapefilename)

    # Create a mask for the data
    mask = np.zeros(cube.shape, dtype=bool)

    # Create a set of x,y points from the cube
    # 1D regular grids
    if cube.coord('longitude').points.ndim < 2:
        x_p, y_p = np.meshgrid(
            cube.coord(axis='X').points, cube.coord(axis='Y').points)
    # 2D irregular grids; spit an error for now
    else:
        logger.error('No fx-files found (sftlf or sftof)!\n \
                     2D grids are suboptimally masked with\n \
                     Natural Earth masks. Exiting.')

    # Wrap around longitude coordinate to match data
    x_p_180 = np.where(x_p >= 180., x_p - 360., x_p)
    # the NE mask has no points at x = -180 and y = +/-90
    # so we will fool it and apply the mask at (-179, -89, 89) instead
    x_p_180 = np.where(x_p_180 == -180., x_p_180 + 1., x_p_180)
    y_p_0 = np.where(y_p == -90., y_p + 1., y_p)
    y_p_90 = np.where(y_p_0 == 90., y_p_0 - 1., y_p_0)

    # Build mask with vectorization
    if len(cube.data.shape) == 3:
        mask[:] = shp_vect.contains(region, x_p_180, y_p_90)
    elif len(cube.data.shape) == 4:
        mask[:, :] = shp_vect.contains(region, x_p_180, y_p_90)

    # Then apply the mask
    if isinstance(cube.data, np.ma.MaskedArray):
        cube.data.mask |= mask
    else:
        cube.data = np.ma.masked_array(cube.data, mask)

    return cube


def polygon_shape(xlist, ylist):
    """
    Make a polygon

    Function that takes a list of x-coordinates and a list of y-coordinates
    and returns a polygon and its (x,y) points on the polygon's border
    """
    from shapely.geometry import Polygon
    poly = Polygon(xlist, ylist)
    x_p, y_p = poly.exterior.coords.xy
    return poly, x_p, y_p


"""
Calculating a custom statistic
==============================

This example shows how to define and use a custom
:class:`iris.analysis.Aggregator`, that provides a new statistical operator for
use with cube aggregation functions such as :meth:`~iris.cube.Cube.collapsed`,
:meth:`~iris.cube.Cube.aggregated_by` or
:meth:`~iris.cube.Cube.rolling_window`.

In this case, we have a time sequence of measurements (time unit dt), and we
want to calculate how many times N the measurements exceed a certain threshold
R over a sliding window dT (multiple of dt). The threshold could be 0 for any
unwanted value for instance.
"""


# Define a function to perform the custom statistical operation.
# Note: in order to meet the requirements of iris.analysis.Aggregator, it must
# do the calculation over an arbitrary (given) data axis.
def count_spells(data, threshold, axis, spell_length):
    """
    Count data occurences

    Function to calculate the number of points in a sequence where the value
    has exceeded a threshold value for at least a certain number of timepoints.

    Generalised to operate on multiple time sequences arranged on a specific
    axis of a multidimensional array.

    Args:

    * data (array):
        raw data to be compared with value threshold.

    * threshold (float):
        threshold point for 'significant' datapoints.

    * axis (int):
        number of the array dimension mapping the time sequences.
        (Can also be negative, e.g. '-1' means last dimension)

    * spell_length (int):
        number of consecutive times at which value > threshold to "count".

    """
    if axis < 0:
        # just cope with negative axis numbers
        axis += data.ndim
    # Threshold the data to find the 'significant' points.
    data_hits = data > threshold
    # Make an array with data values "windowed" along the time axis.
    ###############################################################
    # WARNING: default step is = window size i.e. no overlapping
    # if you want overlapping windows set the step to be m*spell_length
    # where m is a float
    ###############################################################
    hit_windows = rolling_window(
        data_hits, window=spell_length, step=spell_length, axis=axis)
    # Find the windows "full of True-s" (along the added 'window axis').
    full_windows = np.all(hit_windows, axis=axis + 1)
    # Count points fulfilling the condition (along the time axis).
    spell_point_counts = np.sum(full_windows, axis=axis, dtype=int)
    return spell_point_counts


def window_counts(mycube, value_threshold, window_size, pctile):
    """
    Find data counts in a time window

    Function that returns a flat array containing
    the number of data points within a time window `window_size'
    per grid point that satify a condition
    value > value_threshold.
    It also returns statistical measures for the flat array
    window_counts[0] = array
    window_counts[1] = mean(array)
    window_counts[2] = std(array)
    window_counts[3] = percentile(array, pctile)
    """
    # Make an aggregator from the user function.
    spell_count = Aggregator(
        'spell_count', count_spells, units_func=lambda units: 1)

    # Calculate the statistic.
    counts_windowed_cube = mycube.collapsed(
        'time',
        spell_count,
        threshold=value_threshold,
        spell_length=window_size)

    # if one wants to print the whole array
    # np.set_printoptions(threshold=np.nan)
    r_p = counts_windowed_cube.data.flatten()
    meanr = np.mean(r_p)
    stdr = np.std(r_p)
    prcr = np.percentile(r_p, pctile)
    return r_p, meanr, stdr, prcr


def mask_cube_counts(mycube, value_threshold, counts_threshold, window_size):
    """Build the counts mask"""
    # Make an aggregator from the user function.
    spell_count = Aggregator(
        'spell_count', count_spells, units_func=lambda units: 1)

    # Calculate the statistic.
    counts_windowed_cube = mycube.collapsed(
        'time',
        spell_count,
        threshold=value_threshold,
        spell_length=window_size)

    mask = counts_windowed_cube.data >= counts_threshold
    mask.astype(np.int)
    # preserving the original cube metadata
    dummyar = np.ones(mycube.data.shape, dtype=mycube.data.dtype)
    newmask = dummyar * mask
    newmask[newmask == 0] = 1e+20  # np.nan
    masked_cube = mycube.copy()
    # masked_cube.data = masked_cube.data * newmask
    masked_cube.data = newmask
    return counts_windowed_cube, newmask, masked_cube


def mask_threshold(mycube, threshold):
    """
    Mask with threshold

    Takes a MINIMUM value `threshold'
    and removes by masking off anything that's below it in the cube data
    """
    import numpy.ma as ma
    mcube = mycube.copy()
    # apply masking for threshold of MINIMUM value threshold
    mcube.data = ma.masked_less(mycube.data, threshold)
    return mcube


def mask_fillvalues(cubes, threshold_fraction, min_value=-1.e10,
                    time_window=1):
    """Get the final fillvalues mask"""
    # function idea copied from preprocess.py

    # Ensure all cubes have masked arrays
    for cube in cubes:
        cube.data = np.ma.fix_invalid(cube.data, copy=False)

    # Get the fillvalue masks from all datasets
    masks = (_get_fillvalues_mask(cube, threshold_fraction, min_value,
                                  time_window) for cube in cubes)

    # Combine all fillvalue masks
    combined_mask = None
    for mask in masks:
        if combined_mask is None:
            combined_mask = np.zeros_like(mask)
        # Select only valid (not all masked) pressure levels
        n_dims = len(mask.shape)
        if n_dims == 2:
            valid = ~np.all(mask)
            if valid:
                combined_mask |= mask
        elif n_dims == 3:
            valid = ~np.all(mask, axis=(1, 2))
            combined_mask[valid] |= mask[valid]
        else:
            raise NotImplementedError("Unable to handle {} dimensional data"
                                      .format(n_dims))

    if np.any(combined_mask):
        # Apply masks
        logger.debug("Applying fillvalues mask")
        for cube in cubes:
            cube.data.mask |= combined_mask

    return cubes


def _get_fillvalues_mask(cube, threshold_fraction, min_value, time_window):
    # function idea copied from preprocess.py
    logger.debug("Creating fillvalues mask")

    # basic checks
    if threshold_fraction < 0 or threshold_fraction > 1.0:
        raise ValueError(
            "Fraction of missing values {} should be between 0 and 1.0"
            .format(threshold_fraction))
    nr_time_points = len(cube.coord('time').points)
    if time_window > nr_time_points:
        logger.warning("Time window (in time units) larger "
                       "than total time span")

    max_counts_per_time_window = nr_time_points / time_window
    # round to lower integer
    counts_threshold = int(max_counts_per_time_window * threshold_fraction)

    # Make an aggregator
    spell_count = Aggregator(
        'spell_count', count_spells, units_func=lambda units: 1)

    # Calculate the statistic.
    counts_windowed_cube = cube.collapsed(
        'time', spell_count, threshold=min_value, spell_length=time_window)

    # Create mask
    mask = counts_windowed_cube.data < counts_threshold
    if np.ma.isMaskedArray(mask):
        mask = mask.data | mask.mask

    return mask
