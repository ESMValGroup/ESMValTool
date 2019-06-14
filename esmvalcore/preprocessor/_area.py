"""
Area operations on data cubes.

Allows for selecting data subsets using certain latitude and longitude bounds;
selecting geographical regions; constructing area averages; etc.
"""
import logging

import iris
from dask import array as da

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
def extract_region(cube, start_longitude, end_longitude, start_latitude,
                   end_latitude):
    """
    Extract a region from a cube.

    Function that subsets a cube on a box (start_longitude, end_longitude,
    start_latitude, end_latitude)
    This function is a restriction of masked_cube_lonlat().

    Parameters
    ----------
    cube: iris.cube.Cube
        input data cube.
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

    if cube.coord('latitude').ndim == 1:
        region_subset = cube.intersection(
            longitude=(start_longitude, end_longitude),
            latitude=(start_latitude, end_latitude))
        region_subset = region_subset.intersection(longitude=(0., 360.))
        return region_subset
    # irregular grids
    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points
    select_lats = start_latitude < lats < end_latitude
    select_lons = start_longitude < lons < end_longitude
    selection = select_lats & select_lons
    data = da.ma.masked_where(~selection, cube.core_data())
    return cube.copy(data)


def get_iris_analysis_operation(operator):
    """
    Determine the iris analysis operator from a string.

    Map string to functional operator.

    Parameters
    ----------
    operator: str
        A named operator.

    Returns
    -------
        function: A function from iris.analysis

    Raises
    ------
    ValueError
        operator not in allowed operators list.
        allowed operators: mean, median, std_dev, variance, min, max
    """
    operators = ['mean', 'median', 'std_dev', 'variance', 'min', 'max']
    operator = operator.lower()
    if operator not in operators:
        raise ValueError("operator {} not recognised. "
                         "Accepted values are: {}."
                         "".format(operator, ', '.join(operators)))
    operation = getattr(iris.analysis, operator.upper())
    return operation


def zonal_means(cube, coordinate, mean_type):
    """
    Get zonal means.

    Function that returns zonal means along a coordinate `coordinate`;
    the type of mean is controlled by mean_type variable (string):
    - 'mean' -> MEAN
    - 'median' -> MEDIAN
    - 'std_dev' -> STD_DEV
    - 'variance' -> VARIANCE
    - 'min' -> MIN
    - 'max' -> MAX

    Parameters
    ----------
    cube: iris.cube.Cube
        input cube.
    coordinate: str
        name of coordinate to make mean.
    mean_type: str
        Type of analysis to use, from iris.analysis.

    Returns
    -------
    iris.cube.Cube
    """
    operation = get_iris_analysis_operation(mean_type)
    return cube.collapsed(coordinate, operation)


def tile_grid_areas(cube, fx_files):
    """
    Tile the grid area data to match the dataset cube.

    Parameters
    ----------
    cube: iris.cube.Cube
        input cube.
    fx_files: dict
        dictionary of field:filename for the fx_files

    Returns
    -------
    iris.cube.Cube
        Freshly tiled grid areas cube.
    """
    grid_areas = None
    if fx_files:
        for key, fx_file in fx_files.items():
            if fx_file is None:
                continue
            logger.info('Attempting to load %s from file: %s', key, fx_file)
            fx_cube = iris.load_cube(fx_file)

            grid_areas = fx_cube.core_data()
            if cube.ndim == 4 and grid_areas.ndim == 2:
                grid_areas = da.tile(grid_areas,
                                     [cube.shape[0], cube.shape[1], 1, 1])
            elif cube.ndim == 4 and grid_areas.ndim == 3:
                grid_areas = da.tile(grid_areas, [cube.shape[0], 1, 1, 1])
            elif cube.ndim == 3 and grid_areas.ndim == 2:
                grid_areas = da.tile(grid_areas, [cube.shape[0], 1, 1])
            else:
                raise ValueError('Grid and dataset number of dimensions not '
                                 'recognised: {} and {}.'
                                 ''.format(cube.ndim, grid_areas.ndim))
    return grid_areas


# get the area average
def area_statistics(cube, operator, fx_files=None):
    """
    Apply a statistical operator in the horizontal direction.

    The average in the horizontal direction. We assume that the
    horizontal directions are ['longitude', 'latutude'].

    This function can be used to apply
    several different operations in the horizonal plane: mean, standard
    deviation, median variance, minimum and maximum. These options are
    specified using the `operator` argument and the following key word
    arguments:

    +------------+--------------------------------------------------+
    | `mean`     | Area weighted mean.                              |
    +------------+--------------------------------------------------+
    | `median`   | Median (not area weighted)                       |
    +------------+--------------------------------------------------+
    | `std_dev`  | Standard Deviation (not area weighted)           |
    +------------+--------------------------------------------------+
    | `variance` | Variance (not area weighted)                     |
    +------------+--------------------------------------------------+
    | `min`:     | Minimum value                                    |
    +------------+--------------------------------------------------+
    | `max`      | Maximum value                                    |
    +------------+--------------------------------------------------+

    Parameters
    ----------
        cube: iris.cube.Cube
            Input cube.
        operator: str
            The operation, options: mean, median, min, max, std_dev, variance
        fx_files: dict
            dictionary of field:filename for the fx_files

    Returns
    -------
    iris.cube.Cube
        collapsed cube.

    Raises
    ------
    iris.exceptions.CoordinateMultiDimError
        Exception for latitude axis with dim > 2.
    ValueError
        if input data cube has different shape than grid area weights
    """
    grid_areas = tile_grid_areas(cube, fx_files)

    if not fx_files and cube.coord('latitude').points.ndim == 2:
        logger.error(
            'fx_file needed to calculate grid cell area for irregular grids.')
        raise iris.exceptions.CoordinateMultiDimError(cube.coord('latitude'))

    coord_names = ['longitude', 'latitude']
    if grid_areas is None or not grid_areas.any():
        cube = _guess_bounds(cube, coord_names)
        grid_areas = iris.analysis.cartography.area_weights(cube)
        logger.info('Calculated grid area shape: %s', grid_areas.shape)

    if cube.shape != grid_areas.shape:
        raise ValueError('Cube shape ({}) doesn`t match grid area shape '
                         '({})'.format(cube.shape, grid_areas.shape))

    operation = get_iris_analysis_operation(operator)

    # TODO: implement weighted stdev, median, s var when available in iris.
    # See iris issue: https://github.com/SciTools/iris/issues/3208

    if operator == 'mean':
        return cube.collapsed(coord_names,
                              operation,
                              weights=grid_areas)

    # Many IRIS analysis functions do not accept weights arguments.
    return cube.collapsed(coord_names, operation)


def extract_named_regions(cube, regions):
    """
    Extract a specific named region.

    The region coordinate exist in certain CMIP datasets.
    This preprocessor allows a specific named regions to be extracted.

    Parameters
    ----------
    cube: iris.cube.Cube
       input cube.
    regions: str, list
        A region or list of regions to extract.

    Returns
    -------
    iris.cube.Cube
        collapsed cube.

    Raises
    ------
    ValueError
        regions is not list or tuple or set.
    ValueError
        region not included in cube.
    """
    # Make sure regions is a list of strings
    if isinstance(regions, str):
        regions = [regions]

    if not isinstance(regions, (list, tuple, set)):
        raise TypeError(
            'Regions "{}" is not an acceptable format.'.format(regions))

    available_regions = set(cube.coord('region').points)
    invalid_regions = set(regions) - available_regions
    if invalid_regions:
        raise ValueError('Region(s) "{}" not in cube region(s): {}'.format(
            invalid_regions, available_regions))

    constraints = iris.Constraint(region=lambda r: r in regions)
    cube = cube.extract(constraint=constraints)
    return cube
