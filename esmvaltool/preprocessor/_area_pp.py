"""
Area operations on data cubes.

Allows for selecting data subsets using certain latitude and longitude bounds;
selecting geographical regions; constructing area averages; etc.
"""
import logging

import numpy as np
import iris
from iris.coords import AuxCoord
from iris.aux_factory import AuxCoordFactory


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

    try:
        region_subset = cube.intersection(
            longitude=(start_longitude, end_longitude),
            latitude=(start_latitude, end_latitude))
        region_subset = region_subset.intersection(longitude=(0., 360.))
        return region_subset
    except iris.exceptions.CoordinateMultiDimError:
        dims = set(cube.coord_dims('longitude') + cube.coord_dims('latitude'))
        def in_region(lat, lon):
            if start_longitude >= 0:
                if lon < 0:
                    lon += 360
            elif start_longitude < 0:
                if lon > 180:
                    lon -= 360

            if start_latitude > lat  or lat > end_latitude:
                return 0
            elif start_longitude > lon  or lon > end_longitude:
                return 0
            return 1
        vectorized = np.vectorize(in_region)
        in_region = vectorized(
            cube.coord('latitude').points,
            cube.coord('longitude').points,
        )
        lat_lon_dims = list(set(
            cube.coord_dims('latitude') +cube.coord_dims('longitude')
        ))
        lat_lon_dims.sort()
        for index, dimension in enumerate(lat_lon_dims):
            aux_coord = AuxCoord(
                np.any(
                    in_region,
                    axis=tuple(set(range(len(lat_lon_dims))) - {index})
                ),
                var_name='extraction'
            )
            cube.add_aux_coord(aux_coord, dimension)
            cube = cube.extract(iris.Constraint(extraction=True))
            cube.remove_coord('extraction')
        return cube


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
