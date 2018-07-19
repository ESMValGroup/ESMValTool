"""
Area operations on data cubes

Allows for selecting data subsets using certain latitude and longitude bounds;
selecting geographical regions; constructing area averages; etc.
"""
import iris


# slice cube over a restricted area (box)
def area_slice(cube, start_longitude, end_longitude, start_latitude,
               end_latitude):
    """
    Subset a cube on area

    Function that subsets a cube on a box (start_longitude, end_longitude,
    start_latitude, end_latitude)
    This function is a restriction of masked_cube_lonlat();

    Arguments
    ---------
        cube: iris.cube.Cube
            input cube.

        start_longitude: float
            Western boundary longitude.

        end_longitude: float
              Eastern boundary longitude.

        start_latitude: float
             Southern Boundary latitude.

        end_latitude: float
             Northern Boundary Latitude.

    Returns
    -------
    iris.cube.Cube
        smaller cube.
    """
    # Converts Negative longitudes to 0 -> 360. standard
    start_longitude = float(start_longitude)
    end_longitude = float(end_longitude)
    start_latitude = float(start_latitude)
    end_latitude = float(end_latitude)

    if start_longitude < 0.:
        start_longitude += 360.
    if end_longitude < 0.:
        end_longitude += 360.

    if end_longitude < start_longitude:
        # if you want to look at a region both sides of
        # the zero longitude ie, such as the Atlantic Ocean!
        sublon = iris.Constraint(
            longitude=lambda cell: ((start_longitude <= cell <= 360.)
                                    + (0. <= cell <= end_longitude)))
    else:
        sublon = iris.Constraint(
            longitude=lambda cell: start_longitude <= cell <= end_longitude)

    sublat = iris.Constraint(
        latitude=lambda cell: start_latitude <= cell <= end_latitude)

    region_subset = cube.extract(sublon & sublat)
    return region_subset


# get zonal means
def zonal_means(cube, coordinate, mean_type):
    """
    Get zonal means

    Function that returns zonal means along a coordinate `coordinate`;
    the type of mean is controlled by mean_type variable (string):
        'mean' -> MEAN
        'stdev' -> STD_DEV
        'variance' -> VARIANCE
        'min' -> MIN
        'max' -> MAX

    Arguments
    ---------
        cube: iris.cube.Cube
            input cube.

         coordinate: str
             name of coordinate to make mean

         mean_type: str
             Type of analysis to use, from iris.analysis.

    Returns
    -------
    iris.cube.Cube
        Returns a cube
    """
    if mean_type == 'mean':
        result = cube.collapsed(coordinate, iris.analysis.MEAN)
    elif mean_type == 'stdev':
        result = cube.collapsed(coordinate, iris.analysis.STD_DEV)
    elif mean_type == 'variance':
        result = cube.collapsed(coordinate, iris.analysis.VARIANCE)
    elif mean_type.lower() in ['minimum', 'min']:
        result = cube.collapsed(coordinate, iris.analysis.MIN)
    elif mean_type.lower() in ['maximum', 'max']:
        result = cube.collapsed(coordinate, iris.analysis.MAX)
    return result


# get the area average
def area_average(cube, coord1, coord2):
    """
    Determine the area average.

    Can be used with coord1 and coord2 (strings,
    usually 'longitude' and 'latitude' but depends on the cube);

    Arguments
    ---------
        cube: iris.cube.Cube
            input cube.

        coord1: str
            name of first coordinate

        coord2: str
            name of second coordinate

    Returns
    -------
    iris.cube.Cube
        collapsed cube.
    """
    # CMOR ised data should already have bounds?
    #    cube.coord(coord1).guess_bounds()
    #    cube.coord(coord2).guess_bounds()
    for coord in (coord1, coord2):
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(cube)
    result = cube.collapsed(
        [coord1, coord2], iris.analysis.MEAN, weights=grid_areas)
    return result


# operate along a trajectory line
def trajectory_cube(cube, long1, long2, lat1, lat2, plong1, plong2, plat1,
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
    wspd_subset = cube.extract(sublon & sublat)
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
