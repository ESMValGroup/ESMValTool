"""Handle some coordinate functionality of Iris cubes

Handles:

- fixing and adding additional coordinates to a cube

- creating a target grid for regridding

- extract areas from an iris cube

"""

import logging
import numpy as np
import iris
import iris.coords
import iris.coord_categorisation
from .constraints import CoordConstraint


logger = logging.getLogger(__name__)  # pylint: disable=invalid-name


def fixcoords(cube, realization):
    """Add a number of auxiliary coordinates, and adds bounds for
    longitude and latitude.

    The input cube is changed in-place, as well as returned


    Add the following auxiliary coordinates to the input cube:

    - realization: from `realization` argument

    - year: year

    - season: 'djf', 'mam', 'jja', 'son'

    - season_year: years starting in December, so that complete
      seasons fit in a year, and winter is not chopped in parts.

    - month: months

    """

    realization = iris.coords.AuxCoord(realization, 'realization')
    cube.add_aux_coord(realization)

    iris.coord_categorisation.add_season(cube, 'time')
    iris.coord_categorisation.add_season_year(cube, 'time')
    iris.coord_categorisation.add_year(cube, 'time')
    iris.coord_categorisation.add_month(cube, 'time')

    # Longitude and latitude need bounds for calculating area weights
    # (i.e., latitude corrections)
    for name in ('longitude', 'latitude'):
        try:
            cube.coord(name).guess_bounds()
        except ValueError as exc:
            if str(exc) == "Coord already has bounds. Remove the bounds before guessing new ones.":
                pass
            else:
                raise
    return cube


def create_grid():
    """DUMMY DOCSTRING"""
    # TO DO: specify grid size and offset in configuration
    logger.info("Creating grid")
    longitude = np.linspace(0.5, 359.5, 360)
    latitude = np.linspace(-89.5, 89.5, 180)
    data = np.zeros((len(latitude), len(longitude)), dtype=np.float)
    longitude = iris.coords.DimCoord(longitude, standard_name='longitude', units='degrees')
    longitude.guess_bounds()
    latitude = iris.coords.DimCoord(latitude, standard_name='latitude', units='degrees')
    latitude.guess_bounds()

    targetgrid = iris.cube.Cube(data, dim_coords_and_dims=[(latitude, 0), (longitude, 1)])
    return targetgrid


def parse_area(area):
    """Parse areas as defined in the config file & format"""

    if 'latitude' in area:
        lat = area['latitude']
    elif 'lat' in area:
        lat = area['lat']
    elif 's' in area and 'n' in area:
        lat = [area['s'], area['n']]
    else:
        raise KeyError("missing latitude in area definition")
    if 'longitude' in area:
        lon = area['longitude']
    elif 'lon' in area:
        lon = area['lon']
    elif 'long' in area:
        lon = area['long']
    elif 'e' in area and 'w' in area:
        lon = [area['e'], area['w']]
    else:
        raise KeyError("missing latitude in area definition")

    if isinstance(lon, (int, float)) and isinstance(lat, (int, float)):
        return None, [('latitude', lat), ('longitude', lon)]

    if len(lon) != 2:
        raise ValueError("longitude area definition is the wrong format")
    if len(lat) != 2:
        raise ValueError("latitude area definition is the wrong format")
    lon_constraint = CoordConstraint(area['longitude'][0], area['longitude'][1])
    lat_constraint = CoordConstraint(area['latitude'][0], area['latitude'][1])
    constraint = iris.Constraint(longitude=lon_constraint, latitude=lat_constraint)
    return constraint, None


def extract_areas(cube, areas=None, targetgrid=None, average_area=True, gridscheme='area'):
    """DUMMY DOCSTRING"""
    if areas is None:
        areas = {'global': None}
    if targetgrid is not None:
        if gridscheme == 'area':
            scheme = iris.analysis.AreaWeighted()
        elif gridscheme == 'linear':
            scheme = iris.analysis.Linear()
        else:
            scheme = iris.analysis.Linear()
        gridcube = cube.regrid(targetgrid, scheme)
    else:
        gridcube = cube

    cubes = {}
    for name, area in areas.items():
        excube = gridcube.copy()
        if area is not None:
            # pragma pylint: disable=unsupported-membership-test
            # pragma pylint: disable=unsubscriptable-object
            logger.info("Extracting area %s", name)
            # Account for the various types of input:
            # - iris.Constraint
            # - {'latitude': [<north>, <south>], 'longitude': [<west>, <east>]}
            # - {'latitude': <north>, 'longitude': <east>}
            if isinstance(area, iris.Constraint):
                excube = excube.extract(area)
            elif isinstance(area, dict):
                constraint, coords = parse_area(area)
                if coords and not constraint:
                    excube = excube.interpolate(coords, iris.analysis.Linear())
                else:
                    excube = excube.extract(constraint)

        if excube is None:
            logger.warning("Area extraction failed for cube %r", cube)
            cubes[name] = None
            continue

        # Take care not to attempt to average over an "area" of a single grid point
        if average_area and (len(excube.coord('latitude').points) > 1 or
                             len(excube.coord('longitude').points) > 1):
            weights = iris.analysis.cartography.area_weights(excube)
            logger.info("Averaging area %s", name)
            excube_meanarea = excube.collapsed(['latitude', 'longitude'],
                                               iris.analysis.MEAN, weights=weights)
        else:
            excube_meanarea = excube
        cubes[name] = excube_meanarea

    return cubes
