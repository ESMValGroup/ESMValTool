"""
Volume and z coordinate operations on data cubes.

Allows for selecting data subsets using certain volume bounds;
selecting depth or height regions; constructing volumetric averages;
"""
from copy import deepcopy

import logging

import iris
import numpy as np

logger = logging.getLogger(__name__)


def extract_volume(cube, z_min, z_max):
    """
    Subset a cube based on a range of values in the z-coordinate.

    Function that subsets a cube on a box (z_min, z_max)
    This function is a restriction of masked_cube_lonlat();
    Note that this requires the requested z-coordinate range to be the
    same sign as the iris cube. ie, if the cube has z-coordinate as
    negative, then z_min and z_max need to be negative numbers.

    Parameters
    ----------
    cube: iris.cube.Cube
        input cube.
    z_min: float
        minimum depth to extract.
    z_max: float
        maximum depth to extract.

    Returns
    -------
    iris.cube.Cube
        z-coord extracted cube.
    """
    if z_min > z_max:
        # minimum is below maximum, so switch them around
        zmax = float(z_min)
        zmin = float(z_max)
    else:
        zmax = float(z_max)
        zmin = float(z_min)

    z_constraint = iris.Constraint(
        coord_values={
            cube.coord(axis='Z'): lambda cell: zmin < cell.point < zmax})

    return cube.extract(z_constraint)


def _create_cube_time(src_cube, data, times):
    """
    Generate a new cube with the volume averaged data.

    The resultant cube is seeded with `src_cube` metadata and coordinates,
    excluding any source coordinates that span the associated vertical
    dimension. The `times` of interpolation are used along with the
    associated source cube time coordinate metadata to add a new
    time coordinate to the resultant cube.

    Based on the _create_cube method from _regrid.py.

    Parameters
    ----------
    src_cube : cube
        The source cube that was vertically interpolated.
    data : array
        The payload resulting from interpolating the source cube
        over the specified times.
    times : array
        The array of times.

    Returns
    -------
    cube

    .. note::

        If there is only one level of interpolation, the resultant cube
        will be collapsed over the associated vertical dimension, and a
        scalar vertical coordinate will be added.

    """
    # Get the source cube vertical coordinate and associated dimension.
    src_times = src_cube.coord('time')
    t_dim, = src_cube.coord_dims(src_times)

    if data.shape[t_dim] != len(times):
        emsg = ('Mismatch between data and times for data dimension {!r}, '
                'got data shape {!r} with times shape {!r}.')
        raise ValueError(emsg.format(t_dim, data.shape, times.shape))

    # Construct the resultant cube with the interpolated data
    # and the source cube metadata.
    kwargs = deepcopy(src_cube.metadata)._asdict()
    result = iris.cube.Cube(data, **kwargs)

    # Add the appropriate coordinates to the cube, excluding
    # any coordinates that span the z-dimension of interpolation.
    for coord in src_cube.dim_coords:
        [dim] = src_cube.coord_dims(coord)
        if dim != t_dim:
            result.add_dim_coord(coord.copy(), dim)

    for coord in src_cube.aux_coords:
        dims = src_cube.coord_dims(coord)
        if t_dim not in dims:
            result.add_aux_coord(coord.copy(), dims)

    for coord in src_cube.derived_coords:
        dims = src_cube.coord_dims(coord)
        if t_dim not in dims:
            result.add_aux_coord(coord.copy(), dims)

    # Construct the new vertical coordinate for the interpolated
    # z-dimension, using the associated source coordinate metadata.
    kwargs = deepcopy(src_times._as_defn())._asdict()

    try:
        coord = iris.coords.DimCoord(times, **kwargs)
        result.add_dim_coord(coord, t_dim)
    except ValueError:
        coord = iris.coords.AuxCoord(times, **kwargs)
        result.add_aux_coord(coord, t_dim)

    return result


def calculate_volume(cube):
    """
    Calculate volume from a cube.

    This function is used when the volume netcdf fx_files can't be found.

    Parameters
    ----------
    cube: iris.cube.Cube
        input cube.

    Returns
    -------
    float
        grid volume.
    """
    # ####
    # Load depth field and figure out which dim is which.
    depth = cube.coord(axis='z')
    z_dim = cube.coord_dims(cube.coord(axis='z'))[0]

    # ####
    # Load z direction thickness
    thickness = depth.bounds[..., 1] - depth.bounds[..., 0]

    # ####
    # Calculate grid volume:
    area = iris.analysis.cartography.area_weights(cube)
    if thickness.ndim == 1 and z_dim == 1:
        grid_volume = area * thickness[None, :, None, None]
    if thickness.ndim == 4 and z_dim == 1:
        grid_volume = area * thickness[:, :]

    return grid_volume


def volume_statistics(
        cube,
        operator,
        fx_files=None):
    """
    Apply a statistical operation over a volume.

    The volume average is weighted acoording to the cell volume. Cell volume
    is calculated from iris's cartography tool multiplied by the cell
    thickness.

    Parameters
    ----------
        cube: iris.cube.Cube
            Input cube.
        operator: str
            The operation to apply to the cube, options are: 'mean'.
        fx_files: dict
            dictionary of field:filename for the fx_files

    Returns
    -------
    iris.cube.Cube
        collapsed cube.

    Raises
    ------
    ValueError
        if input cube shape differs from grid volume cube shape.
    """
    # TODO: Test sigma coordinates.
    # TODO: Add other operations.

    # ####
    # Load z coordinate field and figure out which dim is which.
    t_dim = cube.coord_dims('time')[0]

    grid_volume_found = False
    grid_volume = None
    if fx_files:
        for key, fx_file in fx_files.items():
            if fx_file is None:
                continue
            logger.info('Attempting to load %s from file: %s', key, fx_file)
            fx_cube = iris.load_cube(fx_file)

            grid_volume = fx_cube.data
            grid_volume_found = True
            cube_shape = cube.data.shape

    if not grid_volume_found:
        grid_volume = calculate_volume(cube)

    # Check whether the dimensions are right.
    if cube.data.ndim == 4 and grid_volume.ndim == 3:
        grid_volume = np.tile(grid_volume,
                              [cube_shape[0], 1, 1, 1])

    if cube.data.shape != grid_volume.shape:
        raise ValueError('Cube shape ({}) doesn`t match grid volume shape '
                         '({})'.format(cube.data.shape, grid_volume.shape))

    # #####
    # Calculate global volume weighted average
    result = []
    # #####
    # iterate over time and z-coordinate dimensions.
    for time_itr in range(cube.shape[t_dim]):
        # ####
        # create empty output arrays
        column = []
        depth_volume = []

        # ####
        # iterate over time and z-coordinate dimensions.
        for z_itr in range(cube.shape[1]):
            # ####
            # Calculate weighted mean for this time and layer
            if operator == 'mean':
                total = cube[time_itr, z_itr].collapsed(
                    [cube.coord(axis='z'),
                     'longitude', 'latitude'],
                    iris.analysis.MEAN,
                    weights=grid_volume[time_itr, z_itr]).data
            else:
                raise ValueError('Volume operator ({}) not '
                                 'recognised.'.format(operator))
            column.append(total)

            try:
                layer_vol = np.ma.masked_where(
                    cube[time_itr, z_itr].data.mask,
                    grid_volume[time_itr, z_itr]).sum()

            except AttributeError:
                # ####
                # No mask in the cube data.
                layer_vol = grid_volume.sum()
            depth_volume.append(layer_vol)
        # ####
        # Calculate weighted mean over the water volumn
        result.append(np.average(column, weights=depth_volume))

    # ####
    # Send time series and dummy cube to cube creating tool.
    times = np.array(cube.coord('time').points.astype(float))
    result = np.array(result)

    # #####
    # Create a small dummy output array for the output cube
    if operator == 'mean':
        src_cube = cube[:2, :2].collapsed([cube.coord(axis='z'),
                                           'longitude', 'latitude'],
                                          iris.analysis.MEAN,
                                          weights=grid_volume[:2, :2], )

    return _create_cube_time(src_cube, result, times)


def depth_integration(cube):
    """
    Determine the total sum over the vertical component.

    Requires a 3D cube. The z-coordinate
    integration is calculated by taking the sum in the z direction of the
    cell contents multiplied by the cell thickness.

    Arguments
    ---------
    cube: iris.cube.Cube
        input cube.

    Returns
    -------
    iris.cube.Cube
        collapsed cube.
    """
    # ####
    depth = cube.coord(axis='z')
    thickness = depth.bounds[..., 1] - depth.bounds[..., 0]

    if depth.ndim == 1:
        slices = [None for i in cube.shape]
        coord_dim = cube.coord_dims(cube.coord(axis='z'))[0]
        slices[coord_dim] = slice(None)
        thickness = np.abs(thickness[tuple(slices)])

    ones = np.ones_like(cube.data)

    weights = thickness * ones

    result = cube.collapsed(cube.coord(axis='z'), iris.analysis.SUM,
                            weights=weights)

    result.rename('Depth_integrated_' + str(cube.name()))
    # result.units = Unit('m') * result.units # This doesn't work:
    # TODO: Change units on cube to reflect 2D concentration (not 3D)
    # Waiting for news from iris community.
    return result


def extract_transect(cube, latitude=None, longitude=None):
    """
    Extract data along a line of constant latitude or longitude.

    Both arguments, latitude and longitude, are treated identically.
    Either argument can be a single float, or a pair of floats, or can be
    left empty.
    The single float indicates the latitude or longitude along which the
    transect should be extracted.
    A pair of floats indicate the range that the transect should be
    extracted along the secondairy axis.

    For instance `'extract_transect(cube, longitude=-28)'` will produce a
    transect along 28 West.

    Also, `'extract_transect(cube, longitude=-28, latitude=[-50, 50])'` will
    produce a transect along 28 West  between 50 south and 50 North.

    This function is not yet implemented for irregular arrays - instead
    try the extract_trajectory function, but note that it is currently
    very slow. Alternatively, use the regrid preprocessor to regrid along
    a regular grid and then extract the transect.

    Parameters
    ----------
    cube: iris.cube.Cube
        input cube.
    latitude: None, float or [float, float], optional
        transect latiude or range.
    longitude:  None, float or [float, float], optional
        transect longitude or range.

    Returns
    -------
    iris.cube.Cube
        collapsed cube.

    Raises
    ------
    ValueError
        slice extraction not implemented for irregular grids.
    ValueError
        latitude and longitute are both floats or lists; not allowed
        to slice on both axes at the same time.
    """
    # ###
    coord_dim2 = False
    second_coord_range = False
    lats = cube.coord('latitude')
    lons = cube.coord('longitude')

    if lats.ndim == 2:
        raise ValueError(
            'extract_slice: Not implemented for irregular arrays!'
            + '\nTry regridding the data first.')

    if isinstance(latitude, float) and isinstance(longitude, float):
        raise ValueError(
            'extract_slice: Cant slice along lat and lon at the same time'
        )

    if isinstance(latitude, list) and isinstance(longitude, list):
        raise ValueError(
            'extract_slice: Can\'t reduce lat and lon at the same time'
        )

    for dim_name, dim_cut, coord in zip(['latitude', 'longitude'],
                                        [latitude, longitude], [lats, lons]):
        # ####
        # Look for the first coordinate.
        if isinstance(dim_cut, float):
            coord_index = coord.nearest_neighbour_index(dim_cut)
            coord_dim = cube.coord_dims(dim_name)[0]

        # ####
        # Look for the second coordinate.
        if isinstance(dim_cut, list):
            coord_dim2 = cube.coord_dims(dim_name)[0]
            second_coord_range = [coord.nearest_neighbour_index(dim_cut[0]),
                                  coord.nearest_neighbour_index(dim_cut[1])]
    # ####
    # Extracting the line of constant longitude/latitude
    slices = [slice(None) for i in cube.shape]
    slices[coord_dim] = coord_index

    if second_coord_range:
        slices[coord_dim2] = slice(second_coord_range[0],
                                   second_coord_range[1])
    return cube[tuple(slices)]


def extract_trajectory(cube, latitudes, longitudes, number_points=2):
    """
    Extract data along a trajectory.

    latitudes and longitudes are the pairs of coordinates for two points.
    number_points is the number of points between the two points.

    This version uses the expensive interpolate method, but it may be
    necceasiry for irregular grids.

    If only two latitude and longitude coordinates are given,
    extract_trajectory will produce a cube will extrapolate along a line
    bewteen those two points, and will add `number_points` points between
    the two corners.

    If more than two points are provided, then
    extract_trajectory will produce a cube which has extrapolated the data
    of the cube to those points, and `number_points` is not needed.

    Parameters
    ----------
    cube: iris.cube.Cube
        input cube.
    latitudes: list
        list of latitude coordinates (floats).
    longitudes: list
        list of longitude coordinates (floats).
    number_points: int
        number of points to extrapolate (optional).

    Returns
    -------
    iris.cube.Cube
        collapsed cube.

    Raises
    ------
    ValueError
        if latitude and longitude have different dimensions.
    """
    from iris.analysis.trajectory import interpolate

    if len(latitudes) != len(longitudes):
        raise ValueError(
            'Longitude & Latitude coordinates have different lengths'
        )

    if len(latitudes) == len(longitudes) == 2:
        minlat, maxlat = np.min(latitudes), np.max(latitudes)
        minlon, maxlon = np.min(longitudes), np.max(longitudes)

        longitudes = np.linspace(minlat, maxlat, num=number_points)
        latitudes = np.linspace(minlon, maxlon, num=number_points)

    points = [('latitude', latitudes), ('longitude', longitudes)]
    interpolated_cube = interpolate(cube, points)  # Very slow!
    return interpolated_cube
