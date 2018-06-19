"""
Time and area operations on data cubes

Allows for selecting data subsets using certain time bounds;
selecting geographical regions; constructing seasonal and area
averages; checks on data time frequencies (daily, monthly etc)
"""
from datetime import timedelta
import iris


# slice cube over a restricted time period
def time_slice(mycube, yr1, mo1, d1, yr2, mo2, d2):
    """
    Slice cube on time

    Function that returns a subset of the original cube (slice)
    given two dates of interest date1 and date2
    date1 and date2 should be given in a yr,mo,d (int)format e.g.
    time_slice(cube,2006,2,2,2010,1,1) or
    time_slice(cube,'2006','2','2','2010','1','1');

    Returns a cube
    """
    import datetime
    time_units = mycube.coord('time').units
    if time_units.calendar == '360_day':
        if d1 > 30:
            d1 = 30
        if d2 > 30:
            d2 = 30
    my_date1 = datetime.datetime(int(yr1), int(mo1), int(d1))
    my_date2 = datetime.datetime(int(yr2), int(mo2), int(d2))

    t1 = time_units.date2num(my_date1)
    t2 = time_units.date2num(my_date2)
    # TODO replace the block below for when using iris 2.0
    # my_constraint = iris.Constraint(time=lambda t: (
    #     t1 < time_units.date2num(t.point) < t2))
    my_constraint = iris.Constraint(time=lambda t: (t1 < t.point < t2))
    cube_slice = mycube.extract(my_constraint)
    return cube_slice


# slice cube over a restricted area (box)
def area_slice(mycube, long1, long2, lat1, lat2):
    """
    Subset a cube on area

    Function that subsets a cube on a box (long1,long2,lat1,lat2)
    This function is a restriction of masked_cube_lonlat();
    Returns a cube
    """

    # Converts Negative longitudes to 0 -> 360. standard
    if long1 < 0.:
        long1 += 360.
    if long2 < 0.:
        long2 += 360.

    if long2 < long1:
        # if you want to look at a region both sides of
        # the zero longitude ie, such as the Atlantic Ocean!
        sublon = iris.Constraint(
            longitude=lambda cell: not ((float(long1) >= cell)
                                        * (cell >= float(long2))))
    else:
        sublon = iris.Constraint(
            longitude=lambda cell: float(long1) <= cell <= float(long2))
    sublat = iris.Constraint(
        latitude=lambda cell: float(lat1) <= cell <= float(lat2))
    region_subset = mycube.extract(sublon & sublat)
    return region_subset


# slice cube over a restricted area (box)
def volume_slice(mycube, long1, long2, lat1, lat2, z1, z2):
    """
    Subset a cube on volume

    Function that subsets a cube on a box (long1,long2,lat1,lat2,z1,z2)
    This function is a restriction of masked_cube_lonlat();
    Returns a cube
    """

    # Converts Negative longitudes to 0 -> 360. standard
    if long1 < 0.:
        long1 += 360.
    if long2 < 0.:
        long2 += 360.

    if long2 < long1:
        # if you want to look at a region both sides of
        # the zero longitude ie, such as the Atlantic Ocean!
        sublon = iris.Constraint(
            longitude=lambda cell:
                not ((float(long1) >= cell) * (cell >= float(long2))))
    else:
        sublon = iris.Constraint(
            longitude=lambda cell: float(long1) <= cell <= float(long2))
    sublat = iris.Constraint(
        latitude=lambda cell: float(lat1) <= cell <= float(lat2))

    subz = iris.Constraint(depth=lambda cell: float(z1) <= cell <= float(z2))

    region_subset = mycube.extract(sublon & sublat & subz)
    return region_subset


# get the time average
def time_average(mycube):
    """Get the time average over MEAN; returns a cube"""
    var_mean = mycube.collapsed('time', iris.analysis.MEAN)
    return var_mean


# get the probability a value is greater than a threshold
def proportion_greater(mycube, coord1, threshold):
    """
    Proportion greater

    Return the probability that a cetain variable coord1 (string)
    is greater than a threshold threshold (float or string),
    across a cube mycube; returns a cube
    """
    thr = float(threshold)
    result = mycube.collapsed(
        coord1, iris.analysis.PROPORTION, function=lambda values: values > thr)
    return result


# get zonal means
def zonal_means(mycube, coord1, mean_type):
    """
    Get zonal means

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
    Determine the area average.

    Can be used with coord1 and coord2 (strings,
    usually 'longitude' and 'latitude' but depends on the cube);
    Returns a cube
    """

    import iris.analysis.cartography
    # CMOR ised data should already have bounds?
    #    mycube.coord(coord1).guess_bounds()
    #    mycube.coord(coord2).guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(mycube)
    result = mycube.collapsed(
        [coord1, coord2], iris.analysis.MEAN, weights=grid_areas)
    return result


# get the volume average
def volume_average(
        mycube,
        coordz,
        coord1,
        coord2,
):
    """
    Determine the area average.

    Can be used with coord1 and coord2 (strings,
    usually 'longitude' and 'latitude' but depends on the cube);

    Returns a cube
    """

    import iris.analysis.cartography
    # CMOR ised data should already have bounds?
    #    mycube.coord(coord1).guess_bounds()
    #    mycube.coord(coord2).guess_bounds()

    depth = mycube.coord(coordz)
    thickness = depth.bounds[:, 1] - depth.bounds[:, 0]  # 1D

    area = iris.analysis.cartography.area_weights(mycube)
    if area.ndim == 4 and thickness.ndim == 1:
        grid_volume = area * thickness[None, :, None, None]
    else:
        grid_volume = area * thickness
    result = mycube.collapsed(
        [coordz, coord1, coord2], iris.analysis.MEAN, weights=grid_volume)

    return result


# get the depth integration
def depth_integration(mycube, coordz):
    """
    Determine the total sum over the vertical component.

    Requires a 3D cube, and the name of the z coordinate.
    Returns a cube 2D.
    """

    depth = mycube.coord(coordz)
    thickness = depth.bounds[..., 1] - depth.bounds[..., 0]  # 1D

    if depth.ndim == 1:
        slices = [None for i in mycube.shape]
        coord_dim = mycube.coord_dims(coordz)[0]
        slices[coord_dim] = slice(None)
        thickness = thickness[tuple(slices)]

    result = mycube * thickness
    result = mycube.collapsed(
        [
            coordz,
        ],
        iris.analysis.SUM,
    )
    result.rename('Depth_integrated_' + str(mycube.name()))

    # This doesn't work:
    # Not able to change units on cube.
    # Waiting for news from iris community.
    # result.units = Unit('m') * result.units
    return result


# extract along a constand latitude or longitude
def extract_slice(mycube, latitude=None, longitude=None):
    """
    Extract data along a constant latitude or longitude.
    A range may also be extracted using a minimum and maximum
    value for latitude or longitude.
    """

    second_coord = False
    lats = mycube.coord('latitude')
    lons = mycube.coord('longitude')

    if lats.ndim == 2:
        raise ValueError(
            'extract_slice: Not implemented for irregular arrays!')

    if isinstance(latitude, float) and isinstance(longitude, float):
        raise ValueError(
            'extract_slice: Cant slice along lat and lon at the same time'
        )

    if isinstance(latitude, list) and isinstance(longitude, list):
        raise ValueError(
            'extract_slice: Can\'t reduce lat and lon at the same time'
        )

    #####
    # Look for the first coordinate.
    if isinstance(latitude, float):
        d = lats.nearest_neighbour_index(latitude)
        coord_dim = mycube.coord_dims('latitude')[0]

    if isinstance(longitude, float):
        d = lons.nearest_neighbour_index(longitude)
        coord_dim = mycube.coord_dims('longitude')[0]

    #####
    # Look for the second coordinate.
    if isinstance(latitude, list):
        second_coord = 'latitude'
        if len(latitude) > 2:
            raise ValueError(
                'extract_slice: latitude slice has too many points: {}'.format(
                    latitude))
        d1 = lats.nearest_neighbour_index(latitude[0])
        d2 = lats.nearest_neighbour_index(latitude[1])
        print('extract_slice(second_coord,', second_coord, d1, d2)

    if isinstance(longitude, list):
        second_coord = 'longitude'
        if len(longitude) > 2:
            raise ValueError(
                'extract_slice: longitude slice has too many points: {}'.
                format(longitude))
        d1 = lons.nearest_neighbour_index(longitude[0])
        d2 = lons.nearest_neighbour_index(longitude[1])
        print('extract_slice(second_coord,', second_coord, d1, d2)

    #####
    # Extracting the line of constant longitude/latitude
    slices = [slice(None) for i in mycube.shape]
    slices[coord_dim] = d

    if second_coord:
        coord_dim2 = mycube.coord_dims(second_coord)[0]
        slices[coord_dim2] = slice(d1, d2)

    newcube = mycube[tuple(slices)]
    print('extract_slice(slicess),', slices, newcube.shape,
          ('was', mycube.shape))
    return newcube


# extract along a trajectory
def extract_trajectory(mycube, latitudes, longitudes, n):
    """
    Extract data along a trajectory.

    """
    from iris.analysis.trajectory import interpolate
    import numpy as np

    if len(latitudes) != len(longitudes):
        raise ValueError(
            'Longitude & Latitude coordinates have different lengths'
        )

    if len(latitudes) == len(longitudes) == 2:
        minlat, maxlat = np.min(latitudes), np.max(latitudes)
        minlon, maxlon = np.min(longitudes), np.max(longitudes)

        longitudes = np.linspace(minlat, maxlat, num=n)
        latitudes = np.linspace(minlon, maxlon, num=n)

    points = [('latitude', latitudes), ('longitude', longitudes)]
    interpolated_cube = interpolate(mycube, points)  # Very slow!
    return interpolated_cube


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
        """Check for three months"""
        return (time.bound[1] - time.bound[0]) == 2160

    three_months_bound = iris.Constraint(time=spans_three_months)
    return annual_seasonal_mean.extract(three_months_bound)


# operate along a trajectory line
def trajectory_cube(mycube, long1, long2, lat1, lat2, plong1, plong2, plat1,
                    plat2, samplecounts):
    """
    Build a trajectory

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
    lon = wspd_subset.coord('longitude').points
    lat = wspd_subset.coord('latitude').points
    return section, lon, lat


# set of time axis checks
# funcs that perform checks on the time axis
# of data cubes and validates the type of data:
# daily, monthly, seasonal or yearly
class NoBoundsError(ValueError):
    """OBS files dont have bounds"""

    pass


def is_daily(cube):
    """Test whether the time coordinate contains only daily bound periods."""

    def is_day(bound):
        """Count days"""
        time_span = timedelta(days=(bound[1] - bound[0]))
        return timedelta(days=1) == time_span

    if not cube.coord('time').has_bounds():
        raise NoBoundsError()
    return all([is_day(bound) for bound in cube.coord('time').bounds])


def is_monthly(cube):
    """A month is a period of at least 28 days, up to 31 days."""

    def is_month(bound):
        """Count months"""
        time_span = timedelta(days=(bound[1] - bound[0]))
        return timedelta(days=31) >= time_span >= timedelta(days=28)

    if not cube.coord('time').has_bounds():
        raise NoBoundsError()
    return all([is_month(bound) for bound in cube.coord('time').bounds])


def is_seasonal(cube):
    """
    Check if data is seasonal

    A season is a period of 3 months, i.e.
    at least 89 days, and up to 92 days.
    """

    def is_season(bound):
        """Count seasons"""
        time_span = timedelta(days=(bound[1] - bound[0]))
        is_seas = timedelta(days=31 + 30 + 31) >= time_span >= \
            timedelta(days=28 + 31 + 30)
        return is_seas

    if not cube.coord('time').has_bounds():
        raise NoBoundsError()
    return all([is_season(bound) for bound in cube.coord('time').bounds])


def is_yearly(cube):
    """A year is a period of at least 360 days, up to 366 days."""

    def is_year(bound):
        """Count years"""
        t_s = timedelta(days=(bound[1] - bound[0]))
        return timedelta(days=365) == t_s or timedelta(days=360) == t_s

    if not cube.coord('time').has_bounds():
        raise NoBoundsError()
    return all([is_year(bound) for bound in cube.coord('time').bounds])
