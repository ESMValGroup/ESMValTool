########################################################################
# PREPROCESSING TOOLS
########################################################################
import iris
import numpy as np
from auxiliary import info
from iris.analysis import Aggregator
from iris.util import rolling_window
import logging

logger = logging.getLogger(__name__)


#########################################################################
# FILE OPERATIONS
#########################################################################
# a couple functions needed by cmor reformatting (the new python one)
def get_attr_from_field_coord(ncfield, coord_name, attr):
    if coord_name is not None:
        attrs = ncfield.cf_group[coord_name].cf_attrs()
        attr_val = [value for (key, value) in attrs if key == attr]
        if attr_val:
            return attr_val[0]
    return None


# Use this callback to fix anything Iris tries to break!
# noinspection PyUnusedLocal
# function also used in preprocess.py
def merge_callback(raw_cube, field, filename):
    # Remove attributes that cause issues with merging and concatenation
    for attr in ['creation_date', 'tracking_id', 'history']:
        if attr in raw_cube.attributes:
            del raw_cube.attributes[attr]
    for coord in raw_cube.coords():
        # Iris chooses to change longitude and latitude units to degrees
        #  regardless of value in file, so reinstating file value
        if coord.standard_name in ['longitude', 'latitude']:
            units = get_attr_from_field_coord(field, coord.var_name, 'units')
            if units is not None:
                coord.units = units


# merge multiple files assigned to a same diagnostic and variable
def glob(file_list, varname):
    """
    Function that takes a list of nc files and globs them into a single one
    """
    # there may be the case where the nc file contains multiple cubes
    # and/or the time calendars are crooked
    # -> these exceptional cases are nicely solved by applying Javier's nice
    # iris fixing (merge_callback and get_attr_from_field_coord are also
    # in preprocess.py but keeping them here in case we will have to
    # change things)

    var_name = varname

    def cube_var_name(raw_cube):
        return raw_cube.var_name == var_name

    var_cons = iris.Constraint(cube_func=cube_var_name)
    # force single cube; this function defaults a list of cubes
    cl = [
        iris.load(a, var_cons, callback=merge_callback)[0] for a in file_list
    ]

    c = iris.cube.CubeList(cl)

    try:
        concatenated = c.concatenate()
        try:
            logger.info("Successfully concatenated cubes")
            return concatenated[0]
        except (OSError, iris.exceptions.IrisError) as exc:
            logger.warning("Could not save concatenated cube, keeping a "
                           "list of files - %s ", exc)
            return 0
    except iris.exceptions.ConcatenateError as exc:
        error_message = "Problem trying to concatenate the following cubes:\n"
        for cube in cl:
            error_message += cube.summary(shorten=True) + '\n'
        logger.warning(
            "Could not concatenate cubes, keeping a list of files - %s",
            error_message)
        return 0


############################################################################
# MASKING
############################################################################
def fx_mask(mycube, fx):
    masked_cube = mycube.copy()
    masked_cube.data = mycube.data * fx.data / 100.
    return masked_cube


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
        cubeslice = mycube.extract(
            iris.Constraint(coord_values={
                coord.standard_name:
                lambda cell: v1 <= cell.point <= v2
            }))
        if cubeslice is not None:
            masked_cubeslice = cubeslice.copy()
            masked_cubeslice.data = ma.masked_greater(cubeslice.data,
                                                      threshold)
            print('Masking cube keeping only what is in between %f and %f' %
                  (v1, v2))
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
    contour_polygons, = [
        contour.geometry for contour in contours
        if contour.attributes[att] == argv
    ]
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

    mycube.data = mycube.data * mask
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
    SPELL_COUNT = Aggregator(
        'spell_count', count_spells, units_func=lambda units: 1)

    # Calculate the statistic.
    counts_windowed_cube = mycube.collapsed(
        'time',
        SPELL_COUNT,
        threshold=value_threshold,
        spell_length=window_size)

    # if one wants to print the whole array
    # np.set_printoptions(threshold=np.nan)
    r = counts_windowed_cube.data.flatten()
    meanr = np.mean(r)
    stdr = np.std(r)
    prcr = np.percentile(r, pctile)
    return r, meanr, stdr, prcr


def mask_cube_counts(mycube, value_threshold, counts_threshold, window_size):

    # Make an aggregator from the user function.
    SPELL_COUNT = Aggregator(
        'spell_count', count_spells, units_func=lambda units: 1)

    # Calculate the statistic.
    counts_windowed_cube = mycube.collapsed(
        'time',
        SPELL_COUNT,
        threshold=value_threshold,
        spell_length=window_size)

    mask = counts_windowed_cube.data > counts_threshold
    mask.astype(np.int)
    # preserving the original cube metadata
    masked_cube = mycube.copy()
    masked_cube.data = mycube.data * mask
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


########################################################################
# TIME AND AREA OPERATIONS
########################################################################


# slice cube over a restricted time period
def time_slice(mycube, yr1, mo1, d1, yr2, mo2, d2):
    """
    Function that returns a subset of the original cube (slice)
    given two dates of interest date1 and date2
    date1 and date2 should be given in a yr,mo,d (int)format e.g.
    time_slice(cube,2006,2,2,2010,1,1) or
    time_slice(cube,'2006','2','2','2010','1','1');

    Returns a cube
    """
    import datetime
    import iris.unit
    myDate1 = datetime.datetime(int(yr1), int(mo1), int(d1))
    myDate2 = datetime.datetime(int(yr2), int(mo2), int(d2))
    t1 = mycube.coord('time').units.date2num(myDate1)
    t2 = mycube.coord('time').units.date2num(myDate2)
    myConstraint = iris.Constraint(time=lambda t: (
        t1 < mycube.coord('time').units.date2num(t.point) and
        t2 > mycube.coord('time').units.date2num(t.point)))
    cubeslice = mycube.extract(myConstraint)
    return cubeslice


# slice cube over a restricted area (box)
def area_slice(mycube, long1, long2, lat1, lat2):
    """
    Function that subsets a cube on a box (long1,long2,lat1,lat2)
    This function is a restriction of masked_cube_lonlat();
    Returns a cube
    """
    sublon = iris.Constraint(
        longitude=lambda cell: float(long1) <= cell <= float(long2))
    sublat = iris.Constraint(
        latitude=lambda cell: float(lat1) <= cell <= float(lat2))
    region_subset = mycube.extract(sublon & sublat)
    return region_subset


# get the time average
def time_average(mycube):
    """
    Function to get the time average over MEAN;
    Returns a cube
    """
    var_mean = mycube.collapsed('time', iris.analysis.MEAN)
    return var_mean


# get the probability a value is greater than a threshold
def proportion_greater(mycube, coord1, threshold):
    """
    Function that returns the probability
    that a cetain variable coord1 (string) is greater than
    a threshold threshold (float or string), across a cube mycube;
    Returns a cube
    """
    thr = float(threshold)
    result = mycube.collapsed(
        coord1, iris.analysis.PROPORTION, function=lambda values: values > thr)
    return result


# get zonal means
def zonal_means(mycube, coord1, mean_type):
    """
    Function that returns zonal means along a coordinate coord1;
    the type of mean is controlled by mean_type variable (string):
        'mean' -> MEAN
        'stdev' -> STD_DEV
        'variance' -> VARIANCE

    Returns a cube
    """
    if mean_type == 'mean':
        result = mycube.collapsed(coord1, iris.analysis.MEAN)
    elif mean_type == 'stdev':
        result = mycube.collapsed(coord1, iris.analysis.STD_DEV)
    elif mean_type == 'variance':
        result = mycube.collapsed(coord1, iris.analysis.VARIANCE)
    return result


# get the area average
def area_average(mycube, coord1, coord2):
    """
    Function that determines the area average
    Can be used with coord1 and coord2 (strings,
    usually 'longitude' and 'latitude' but depends on the cube);
    Returns a cube
    """
    import iris.analysis.cartography
    mycube.coord(coord1).guess_bounds()
    mycube.coord(coord2).guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(mycube)
    result = mycube.collapsed(
        [coord1, coord2], iris.analysis.MEAN, weights=grid_areas)
    return result


# get the seasonal mean
def seasonal_mean(mycube):
    """
    Function to compute seasonal means with MEAN
    Chunks time in 3-month periods and computes means over them;
    Returns a cube
    """
    import iris.coord_categorisation
    iris.coord_categorisation.add_season(mycube, 'time', name='clim_season')
    iris.coord_categorisation.add_season_year(
        mycube, 'time', name='season_year')
    annual_seasonal_mean = mycube.aggregated_by(['clim_season', 'season_year'],
                                                iris.analysis.MEAN)

    def spans_three_months(time):
        return (time.bound[1] - time.bound[0]) == 2160

    three_months_bound = iris.Constraint(time=spans_three_months)
    resc = annual_seasonal_mean.extract(three_months_bound)
    print(resc)
    return resc


# operate along a trajectory line
def trajectory_cube(mycube, long1, long2, lat1, lat2, plong1, plong2, plat1,
                    plat2, samplecounts):
    """
    Function that subsets a cube on a box (long1,long2,lat1,lat2)
    then creates a trajectory with waypoints (plong1,plong2,plat1, plat2),
    populates it with samplecounts number of points
    and subsets the cube along the trajectory
    """
    from iris.analysis import trajectory
    sublon = iris.Constraint(
        longitude=lambda cell: float(long1) <= cell <= float(long2))
    sublat = iris.Constraint(
        latitude=lambda cell: float(lat1) <= cell <= float(lat2))
    wspd_subset = mycube.extract(sublon & sublat)
    pnts = [{
        'longitude': float(plong1),
        'latitude': float(plat1)
    }, {
        'longitude': float(plong2),
        'latitude': float(plat2)
    }]
    traj = trajectory.Trajectory(pnts, sample_count=int(samplecounts))
    lon = [d['longitude'] for d in traj.sampled_points]
    lat = [d['latitude'] for d in traj.sampled_points]
    sampled_points = [('longitude', lon), ('latitude', lat)]
    section = trajectory.interpolate(wspd_subset, sampled_points)
    lon, lat = wspd_subset.coord('longitude').points, wspd_subset.coord(
        'latitude').points
    return section, lon, lat
