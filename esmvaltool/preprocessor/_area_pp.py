"""
Area operations on data cubes

Allows for selecting data subsets using certain latitude and longitude bounds;
selecting geographical regions; constructing area averages; etc.
"""
import iris


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
    # CMOR ised data should already have bounds?
    #    mycube.coord(coord1).guess_bounds()
    #    mycube.coord(coord2).guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(mycube)
    result = mycube.collapsed(
        [coord1, coord2], iris.analysis.MEAN, weights=grid_areas)
    return result


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
