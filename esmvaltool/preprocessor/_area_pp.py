"""
Area operations on data cubes.

Allows for selecting data subsets using certain latitude and longitude bounds;
selecting geographical regions; constructing area averages; etc.
"""
import logging
import numpy as np
import numpy.ma as ma

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


def _area_average_general_cecks(coords, aggregator, weighted=False):
    """Perform a set of pre-checks before area_average_general."""
    if coords is None:
        coords = ['latitude', 'longitude']

    # Make sure that aggregator is an Aggregator instance
    if not isinstance(aggregator, iris.analysis.Aggregator):
        logger.warning('Aggregator %s not an iris aggregator, \
                       performing MEAN instead', aggregator)
        aggregator = iris.analysis.MEAN
    # If doing weighted aggregation make sure that aggregator
    # is a WeightAggregator instance
    if weighted:
        if not isinstance(aggregator, iris.analysis.WeightedAggregator):
            logger.warning('Aggregator %s not a weighted iris aggregator, \
                           performing MEAN instead', aggregator)
            aggregator = iris.analysis.MEAN

    return coords, aggregator


def _apply_area_average_mask(cube, aggkeys, mask,
                             logicmask=False, intkeys=None):
    """Apply the needed type of mask while area averaging."""
    # Extract region of mask to match data
    if intkeys:
        newmask = mask.intersection(ignore_bounds=True, **intkeys)
    else:
        newmask = mask.copy()
    if 'weights' in aggkeys:
        if logicmask:
            aggkeys['weights'] = ma.array(
                data=aggkeys['weights'], mask=newmask.data)
        else:
            aggkeys['weights'] *= newmask.data
    else:
        if logicmask:
            cube.data = ma.array(data=cube.data, mask=newmask.data)
        else:
            cube.data *= newmask.data

    return cube


def area_average_general(cube,
                         weighted=True,
                         mask=None,
                         logicmask=False,
                         coords=None,
                         aggregator=iris.analysis.MEAN,
                         **aggkeys):
    """
    Routine to calculate weighted horizontal area aggregations.

    Routine defaults to longitude and latitude, but can be configured to
    collapse over any coordinate in the cube.

    Inputs:

    cube = cube to aggregate

    Keywords:

    weighted = perform area weighted aggregation (default: True)
    mask = cube containing mask data (default: None)
    logicmask = Does mask contain logical data (default: False)
    aggregator = aggregator for collapsed method (default: iris.analysis.MEAN)
    coords = list of coordinates to collapse cube over
             (default: ["latitude", "longitude"])
    "coord" = (coord_min, coord_max) - range of coordinate to collapse over
    **kwargs = any keywords required for the aggregator

    Return:

    aggregated cube.
    """
    coords, aggregator = _area_average_general_cecks(coords,
                                                     aggregator, weighted)
    # Extract region specification if available
    intkeys = {}
    for coord in coords:
        if coord in aggkeys:
            intkeys[coord] = aggkeys.pop(coord)

    # Extract region if required
    if intkeys:
        newcube = cube.intersection(ignore_bounds=True, **intkeys)
        # For some reason cube.intersection() promotes dtype of coordinate
        # arrays to float64, whereas cube.extract() doesn't. Need to make
        # sure behaviour is identical.
        for coord in intkeys.keys():
            newcube.coord(coord).points = \
                newcube.coord(coord).points.astype(np.float32, copy=False)
    else:
        newcube = cube.copy()

    # If doing area-weighted aggregation then calculate area weights
    if weighted:
        # Coords need bounding
        newcube = _guess_bounds(newcube, coords)
        aggkeys['weights'] = iris.analysis.cartography.area_weights(newcube)

    # Apply mask
    if mask:
        _apply_area_average_mask(cube, aggkeys, mask, logicmask, intkeys)

    return newcube.collapsed(coords, aggregator, **aggkeys)
