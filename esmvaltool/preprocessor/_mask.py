"""
_mask.py

module that performs missing values masking
and geographical area eslection
"""

import logging
import os

import cartopy.io.shapereader as shpreader
import iris
import numpy as np
import shapely.vectorized as shp_vect
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
    if (x_dim == mx_dim and y_dim == my_dim and len_x == len_mx
            and len_y == len_my):
        logger.debug('Data cube and fx mask have same dims')
        return True

    logger.debug(
        'Data cube and fx mask differ in dims: '
        'cube: ((%i, %i), grid=(%i, %i)), mask: ((%i, %i), grid=(%i, %i))',
        x_dim, y_dim, len_x, len_y, mx_dim, my_dim, len_mx, len_my)
    return False


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
    elif mask_type == 'sftgif':
        if fx_option == 'ice':
            # Mask ice out
            inmask[fx_data > 50.] = True
        elif fx_option == 'landsea':
            # Mask landsea out
            inmask[fx_data <= 50.] = True

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


def mask_landsea(cube, fx_files, mask_out):
    """
    Mask out either land or sea

    Function that masks out either land mass or seas (oceans, seas and lakes)

    It uses dedicated fx files (sftlf or sftof) or, in their absence, it
    applies a Natural Earth mask (land or ocean contours). Not that the
    Natural Earth masks have different resolutions: 10m for land, and 50m
    for seas; these are more than enough for ESMValTool puprpose.

    Parameters
    ----------

    * cube (iris.Cube.cube instance):
        data cube to be masked.

    * fx_files (list):
        list holding the full paths to fx files.

    * mask_out (string):
        either "land" to mask out land mass or "sea" to mask out seas.

    Returns
    -------
    masked iris cube

    """
    # Dict to store the Natural Earth masks
    cwd = os.path.dirname(__file__)
    # ne_10m_land is fast; ne_10m_ocean is very slow
    shapefiles = {
        'land': os.path.join(cwd, 'ne_masks/ne_10m_land.shp'),
        'sea': os.path.join(cwd, 'ne_masks/ne_50m_ocean.shp')
    }

    if fx_files:
        fx_cubes = {}
        for fx_file in fx_files:
            fx_root = os.path.basename(fx_file).split('_')[0]
            fx_cubes[fx_root] = iris.load_cube(fx_file)

        # preserve importance order: try stflf first then sftof
        if ('sftlf' in fx_cubes.keys()
                and _check_dims(cube, fx_cubes['sftlf'])):
            landsea_mask = _get_fx_mask(fx_cubes['sftlf'].data, mask_out,
                                        'sftlf')
            cube.data = _apply_fx_mask(landsea_mask, cube.data)
            logger.debug("Applying land-sea mask: sftlf")
        elif ('sftof' in fx_cubes.keys()
              and _check_dims(cube, fx_cubes['sftof'])):
            landsea_mask = _get_fx_mask(fx_cubes['sftof'].data, mask_out,
                                        'sftof')
            cube.data = _apply_fx_mask(landsea_mask, cube.data)
            logger.debug("Applying land-sea mask: sftof")
        else:
            if cube.coord('longitude').points.ndim < 2:
                cube = _mask_with_shp(cube, shapefiles[mask_out])
                logger.debug(
                    "Applying land-sea mask from Natural Earth"
                    " shapefile: \n%s", shapefiles[mask_out])
            else:
                logger.error("Use of shapefiles with irregular grids not "
                             "yet implemented, land-sea mask not applied")
    else:
        if cube.coord('longitude').points.ndim < 2:
            cube = _mask_with_shp(cube, shapefiles[mask_out])
            logger.debug(
                "Applying land-sea mask from Natural Earth"
                " shapefile: \n%s", shapefiles[mask_out])
        else:
            logger.error("Use of shapefiles with irregular grids not "
                         "yet implemented, land-sea mask not applied")

    return cube


def mask_landseaice(cube, fx_files, mask_out):
    """
    Mask out either landsea (combined) or ice

    Function that masks out either landsea (land and seas) or ice (Antarctica
    and Greenland and some wee glaciers). It uses dedicated fx files (sftgif).

    Parameters
    ----------

    * cube (iris.Cube.cube instance):
        data cube to be masked.

    * fx_files (list):
        list holding the full paths to fx files.

    * mask_out (string):
        either "landsea" to mask out landsea or "ice" to mask out ice.

    Returns
    -------
    masked iris cube

    """
    # sftgif is the only one so far
    if fx_files:
        for fx_file in fx_files:
            fx_cube = iris.load_cube(fx_file)

            if _check_dims(cube, fx_cube):
                landice_mask = _get_fx_mask(fx_cube.data, mask_out, 'sftgif')
                cube.data = _apply_fx_mask(landice_mask, cube.data)
                logger.debug("Applying landsea-ice mask: sftgif")
    else:
        logger.warning("Landsea-ice mask could not be found ")

    return cube


def _get_geometry_from_shp(shapefilename):
    """Get the mask geometry out from a shapefile"""
    reader = shpreader.Reader(shapefilename)
    # Index 0 grabs the lowest resolution mask (no zoom)
    main_geom = [contour for contour in reader.geometries()][0]
    return main_geom


def _mask_with_shp(cube, shapefilename):
    """Apply a Natural Earth land/sea mask"""
    # Create the region
    region = _get_geometry_from_shp(shapefilename)

    # Create a mask for the data
    mask = np.zeros(cube.shape, dtype=bool)

    # Create a set of x,y points from the cube
    # 1D regular grids
    if cube.coord('longitude').points.ndim < 2:
        x_p, y_p = np.meshgrid(
            cube.coord(axis='X').points,
            cube.coord(axis='Y').points)
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


def mask_above_threshold(cube, threshold):
    """
    Mask above a specific threshold value.

    Takes a value 'threshold' and masks off anything that is above
    it in the cube data. Values equal to the threshold are not masked.
    """
    cube.data = np.ma.masked_where(cube.data > threshold, cube.data)
    return cube


def mask_below_threshold(cube, threshold):
    """
    Mask below a specific threshold value.

    Takes a value 'threshold' and masks off anything that is below
    it in the cube data. Values equal to the threshold are not masked.
    """
    cube.data = np.ma.masked_where(cube.data < threshold, cube.data)
    return cube


def mask_inside_range(cube, minimum, maximum):
    """
    Mask inside a specific threshold range.

    Takes a MINIMUM and a MAXIMUM value for the range, and masks off anything
    that's between the two in the cube data.
    """
    cube.data = np.ma.masked_inside(cube.data, minimum, maximum)
    return cube


def mask_outside_range(cube, minimum, maximum):
    """
    Mask outside a specific threshold range.

    Takes a MINIMUM and a MAXIMUM value for the range, and masks off anything
    that's outside the two in the cube data.
    """
    cube.data = np.ma.masked_outside(cube.data, minimum, maximum)
    return cube


def mask_fillvalues(products,
                    threshold_fraction,
                    min_value=-1.e10,
                    time_window=1):
    """Compute and apply a multi-dataset fillvalues mask"""
    combined_mask = None

    logger.debug("Creating fillvalues mask")
    used = set()
    for product in products:
        for cube in product.cubes:
            cube.data = np.ma.fix_invalid(cube.data, copy=False)
            mask = _get_fillvalues_mask(cube, threshold_fraction, min_value,
                                        time_window)
            if combined_mask is None:
                combined_mask = np.zeros_like(mask)
            # Select only valid (not all masked) pressure levels
            n_dims = len(mask.shape)
            if n_dims == 2:
                valid = ~np.all(mask)
                if valid:
                    combined_mask |= mask
                    used.add(product)
            elif n_dims == 3:
                valid = ~np.all(mask, axis=(1, 2))
                combined_mask[valid] |= mask[valid]
                if np.any(valid):
                    used.add(product)
            else:
                raise NotImplementedError(
                    "Unable to handle {} dimensional data".format(n_dims))

    if np.any(combined_mask):
        logger.debug("Applying fillvalues mask")
        used = {p.copy_provenance() for p in used}
        for product in products:
            for cube in product.cubes:
                cube.data.mask |= combined_mask
            for other in used:
                if other.filename != product.filename:
                    product.wasderivedfrom(other)

    return products


def _get_fillvalues_mask(cube, threshold_fraction, min_value, time_window):

    # basic checks
    if threshold_fraction < 0 or threshold_fraction > 1.0:
        raise ValueError(
            "Fraction of missing values {} should be between 0 and 1.0".format(
                threshold_fraction))
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
