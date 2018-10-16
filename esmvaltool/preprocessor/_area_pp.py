"""
Area operations on data cubes.

Allows for selecting data subsets using certain latitude and longitude bounds;
selecting geographical regions; constructing area averages; etc.
"""
import logging

import iris


logger = logging.getLogger(__name__)


# guess bounds tool
def _guess_bounds(cube, coords):
    """Guess bounds of a cube, or not."""
    # check for bounds just in case
    for coord in coords:
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    return cube


# slice cube over a restricted area (box)
def area_slice(cube, start_longitude, end_longitude, start_latitude,
               end_latitude):
    """
    Subset a cube on area.

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

    region_subset = cube.intersection(
        longitude=(start_longitude, end_longitude),
        latitude=(start_latitude, end_latitude))
    region_subset = region_subset.intersection(longitude=(0., 360.))

    return region_subset


# get zonal means
def zonal_means(cube, coordinate, mean_type):
    """
    Get zonal means.

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

        coord1, coord2: str, str
            coords to use

    Returns
    -------
    iris.cube.Cube
        collapsed cube.
    """
    # check for bounds just in case
    coords = [coord1, coord2]
    cube = _guess_bounds(cube, coords)
    grid_areas = iris.analysis.cartography.area_weights(cube)
    result = cube.collapsed(coords, iris.analysis.MEAN, weights=grid_areas)
    return result


# to be used by the generalized area average
# def _area_average_general_cecks(coords, aggregator, weighted=False):
#    """Perform a set of pre-checks before area_average_general."""
#    if coords is None:
#        coords = ['latitude', 'longitude']
#
#    # Make sure that aggregator is an Aggregator instance
#    if not isinstance(aggregator, iris.analysis.Aggregator):
#        raise ValueError('Aggregator not an iris aggregator!')
#    # If doing weighted aggregation make sure that aggregator
#    # is a WeightAggregator instance
#    if weighted:
#        if not isinstance(aggregator, iris.analysis.WeightedAggregator):
#            raise ValueError('Aggregator not a weighted iris aggregator!')
#
#    return coords, aggregator


# this is a function that will have to be implemented in the future
# in lieu of the less general area average
# def area_average_general(cube,
#                         weighted=True,
#                         coords=None,
#                         aggregator=iris.analysis.MEAN):
#    """
#    Routine to calculate weighted horizontal area aggregations.
#
#    Routine defaults to longitude and latitude, but can be configured to
#    collapse over any coordinate in the cube.
#
#    Inputs:
#
#    cube = cube to aggregate
#
#    Keywords:
#
#    weighted = perform area weighted aggregation (default: True)
#    aggregator = aggregator for collapsed method (default: iris.analysis.MEAN)
#    coords = list of coordinates to collapse cube over
#             (default: ["latitude", "longitude"])
#
#    Return:
#
#    aggregated cube.
#    """
#    coords, aggregator = _area_average_general_cecks(coords,
#                                                     aggregator, weighted)
#
#    # If doing area-weighted aggregation then calculate area weights
#    if weighted:
#        # Coords need bounding
#        cube = _guess_bounds(cube, coords)
#        weights = iris.analysis.cartography.area_weights(cube)
#
#    return cube.collapsed(coords, aggregator, weights=weights)
