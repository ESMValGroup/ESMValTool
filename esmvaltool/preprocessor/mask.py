########################################################################
# MASKING
# V.Predoi, University of Reading, May 2017
########################################################################
import iris

def masked_cube_simple(mycube, slicevar, v1, v2, threshold):
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
    	cubeslice = mycube.extract(iris.Constraint(coord_values = {coord.standard_name:lambda cell: v1 <= cell.point <= v2}))
    	if cubeslice is not None:
        	masked_cubeslice = cubeslice.copy()
        	masked_cubeslice.data = ma.masked_greater(cubeslice.data, threshold)
        	print('Masking cube keeping only what is in between %f and %f'% (v1, v2))
        	return masked_cubeslice
    	else:
        	print('NOT masking the cube')
        	return mycube
    else:
        print('Variable is not a cube dimension, leaving cube untouched')
        return mycube

def masked_cube_lonlat(mycube, lon1, lon2, lat1, lat2, threshold):
    """

    Mask function 2 -- simple cube cropping on (min,max) lon,lat
    Buikds a box and keeps only the values inside the box
    arguments: cube, min value, max value, where value=(lon, lat), threshold

    """
    import numpy.ma as ma
    cubeslice = mycube.extract(iris.Constraint(longitude = lambda v: lon1 <= v.point <= lon2, latitude = lambda v: lat1 <= v.point <= lat2))
    if cubeslice is not None:
        masked_cubeslice = cubeslice.copy()
        masked_cubeslice.data = ma.masked_greater(cubeslice.data, threshold)
        print('Masking cube on lon-lat')
        return masked_cubeslice
    else:
        print('NOT masking the cube')
        return mycube

def cube_shape(mycube):
    """

    Function that converts a cube into a shapely MultiPoint geometry

    """
    import shapely.geometry as sg
    lon = mycube.coord('longitude')
    lat = mycube.coord('latitude')
    region = sg.MultiPoint(zip(lon.points.flat, lat.points.flat))
    return region

def maskgeometry(shapefilename, att, argv):

    """
    This function takes in a shapefile shapefilename
    and creates a specific geometry based on a set of conditions
    on contour attributes att is argv e.g.
    contour.attributes['name'] == 'land_mass'

    """

    import cartopy.io.shapereader as shpreader
    reader = shpreader.Reader(shapefilename)
    contours = reader.records()
    contour_polygons, = [contour.geometry for contour in contours if
                         contour.attributes[att] == argv]
    main_geom = sorted(contour_polygons.geoms, key=lambda geom: geom.area)[-1]
    return main_geom

def mask_2d(mycube, geom):
    """
    This function masks off any given 2D geometry geom
    and keeps only the values that fall in geom, nulling
    everything else in mycube
    WARNING: as of now this function works with cubes that have time as coord
    it is adusted to save time and loop over only lon-lat points

    """
    import numpy as np
    from shapely.geometry import Point
    ccube = mycube.collapsed('time', iris.analysis.MEAN)
    mask = np.ones(ccube.data.shape)
    p = -1
    for i in np.ndindex(ccube.data.shape):
        if i[0] != p:
            print i[0],
            p = i[0]
        this_cube = ccube[i]
        this_lat = this_cube.coord('latitude').points[0]
        this_lon = this_cube.coord('longitude').points[0] - 360
        this_point = Point(this_lon, this_lat)
        mask[i] = this_point.within(geom)

    mycube.data = mycube.data*mask
    return mycube

def polygon_shape(xlist, ylist):
    """

    Function that takes a list of x-coordinates and a list of y-coordinates
    and returns a polygon and its (x,y) points on the polygon's border

    """
    from shapely.geometry import Polygon
    poly = Polygon(xlist, ylist)
    x, y = poly.exterior.coords.xy
    return poly, x, y


"""
Calculating a custom statistic
==============================

This example shows how to define and use a custom
:class:`iris.analysis.Aggregator`, that provides a new statistical operator for
use with cube aggregation functions such as :meth:`~iris.cube.Cube.collapsed`,
:meth:`~iris.cube.Cube.aggregated_by` or
:meth:`~iris.cube.Cube.rolling_window`.

In this case, we have a time sequence of measurements (time unit dt), and we want to calculate how many times N
the measurements exceed a certain threshold R over a sliding window dT (multiple of dt). The threshold could be 0
for any unwanted value for instance.

"""
import numpy as np
from iris.analysis import Aggregator
from iris.util import rolling_window


# Define a function to perform the custom statistical operation.
# Note: in order to meet the requirements of iris.analysis.Aggregator, it must
# do the calculation over an arbitrary (given) data axis.
def count_spells(data, threshold, axis, spell_length):
    """
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
    data_hits = data>threshold
    # Make an array with data values "windowed" along the time axis.
    ###############################################################
    # WARNING: default step is = window size i.e. no overlapping
    # if you want overlapping windows set the step to be m*spell_length
    # where m is a float
    ###############################################################
    hit_windows = rolling_window(data_hits, window=spell_length, step=spell_length, axis=axis)
    # Find the windows "full of True-s" (along the added 'window axis').
    full_windows = np.all(hit_windows, axis=axis+1)
    # Count points fulfilling the condition (along the time axis).
    spell_point_counts = np.sum(full_windows, axis=axis, dtype=int)
    return spell_point_counts

def window_counts(mycube, value_threshold, window_size, pctile):
    """
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
    SPELL_COUNT = Aggregator('spell_count',
                             count_spells,
                             units_func=lambda units: 1)

    # Calculate the statistic.
    counts_windowed_cube = mycube.collapsed('time', SPELL_COUNT,
                                          threshold=value_threshold,
                                          spell_length=window_size)

    #if one wants to print the whole array
    #np.set_printoptions(threshold=np.nan)
    r = counts_windowed_cube.data.flatten()
    meanr = np.mean(r)
    stdr = np.std(r)
    prcr = np.percentile(r, pctile)
    return r,meanr,stdr,prcr

def mask_cube_counts(mycube, value_threshold, counts_threshold, window_size):
    # Make an aggregator from the user function.
    SPELL_COUNT = Aggregator('spell_count',
                             count_spells,
                             units_func=lambda units: 1)

    # Calculate the statistic.
    counts_windowed_cube = mycube.collapsed('time', SPELL_COUNT,
                                          threshold=value_threshold,
                                          spell_length=window_size)

    mask = counts_windowed_cube.data > counts_threshold
    mask.astype(np.int)
    # preserving the original cube metadata
    masked_cube = mycube.copy()
    masked_cube.data = mycube.data*mask
    return counts_windowed_cube, mask, masked_cube

def mask_threshold(mycube, threshold):
    """
    Takes a MINIMUM value `threshold'
    and removes by masking off anything that's below it in the cube data
    """
    import numpy.ma as ma
    mcube = mycube.copy()
    # apply masking for threshold of MINIMUM value threshold
    mcube.data = ma.masked_less(mycube.data, threshold)
    return mcube
