"""
Volume and z coordinate operations on data cubes.

Allows for selecting data subsets using certain volume bounds;
selecting depth or height regions; constructing volumetric averages;
"""
from copy import deepcopy

import iris
import numpy as np


def volume_slice(cube, z_min, z_max):
    """
    Subset a cube on volume

    Function that subsets a cube on a box (z_min, z_max)
    This function is a restriction of masked_cube_lonlat();
    Note that this requires the requested depth range to be the same sign
    as the iris cube. ie, if the cube has depth as negative, then z_min
    and z_max need to be negative numbers.

    Arguments
    ---------
        cube: iris.cube.Cube
            input cube.

        z_min: float
            minimum depth to extract.

        z_max: float
            maximum depth to extract.

    Returns
    -------
    iris.cube.Cube
        extracted cube.
    """
    if z_min > z_max:
        # minimum is below maximum, so switch them around
        zmax = z_min
        zmin = z_max
    else:
        zmax = z_max
        zmin = z_min

    subz = iris.Constraint(
        depth=lambda cell: float(zmin) <= cell <= float(zmax))

    region_subset = cube.extract(subz)
    return region_subset


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
    print(src_cube.coords)
    print(len(data))
    print(times)
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

    # Collapse the z-dimension for the scalar case.
    if times.size == 1:
        slicer = [slice(None)] * result.ndim
        slicer[t_dim] = 0
        result = result[tuple(slicer)]

    return result


def volume_average(cube, coordz, coord1, coord2):
    """
    Determine the volume average.

    The volume average is weighted acoording to the cell volume. Cell volume
    is calculated from iris's cartography tool multiplied by the cell
    thickness.

    Arguments
    ---------
        cube: iris.cube.Cube
            input cube.

        coordz: str
            name of depth coordinate

        coord1: str
            name of first coordinate

        coord2: str
            name of second coordinate

    Returns
    -------
    iris.cube.Cube
        collapsed cube.
    """
    # TODO: Add sigma depth coordinates.
    # TODO: Calculate cell volume, but it may be already included in netcdf.

    # ####
    # Load depth field and figure out which dim is which.
    depth = cube.coord(coordz)
    t_dim = cube.coord_dims('time')[0]
    z_dim = cube.coord_dims(coordz)[0]

    # ####
    # Load z direction thickness
    thickness = depth.bounds[..., 1] - depth.bounds[..., 0]

    if cube.shape[t_dim] < 2:
        # ####
        # Very small cube and should be able to do it in one.
        # Too small to make a dummy cube.
        area = iris.analysis.cartography.area_weights(cube)

        # ####
        # Calculate grid volume:
        if thickness.ndim == 1 and z_dim == 1:
            grid_volume = area * thickness[None, :, None, None]
        if thickness.ndim == 4:
            grid_volume = area * thickness

        return cube.collapsed([coordz, coord1, coord2],
                              iris.analysis.MEAN,
                              weights=grid_volume, )

    # ####
    # Calculate grid volume:
    area = iris.analysis.cartography.area_weights(cube[:2, :2])
    if thickness.ndim == 1 and z_dim == 1:
        grid_volume = area * thickness[None, :2, None, None]
    if thickness.ndim == 4 and z_dim == 1:
        grid_volume = area * thickness[:, :2]

    # #####
    # Create a small dummy output array
    src_cube = cube[:2, :2].collapsed([coordz, coord1, coord2],
                                      iris.analysis.MEAN,
                                      weights=grid_volume, )

    # #####
    # Calculate global volume weighted average
    result = []
    # #####
    # iterate over time and depth dimensions.
    for time_itr in range(cube.shape[t_dim]):
        # ####
        # create empty output arrays
        column = []
        depth_volume = []

        # ####
        # assume cell area is the same thoughout the water column
        area = iris.analysis.cartography.area_weights(cube[time_itr, 0])

        # ####
        # iterate over time and depth dimensions.
        for z_itr in range(cube.shape[1]):
            # ####
            # Calculate grid volume for this time and layer

            if thickness.ndim == 1:
                grid_volume = area * thickness[z_itr]
            if thickness.ndim == 4:
                grid_volume = area * thickness[time_itr, z_itr]

            # ####
            # Calculate weighted mean for this time and layer
            total = cube[time_itr, z_itr].collapsed([coordz, coord1, coord2],
                                                    iris.analysis.MEAN,
                                                    weights=grid_volume).data
            column.append(total)

            try:
                layer_vol = np.ma.masked_where(cube[time_itr, z_itr].data.mask,
                                               grid_volume).sum()
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
    return _create_cube_time(src_cube, result, times)


def depth_integration(cube, coordz):
    """
    Determine the total sum over the vertical component.

    Requires a 3D cube, and the name of the z coordinate. The depth
    integration is calculated by taking the sum in the z direction
    of the cell contents multiplied by the cell thickness.

    Arguments
    ---------
        cube: iris.cube.Cube
            input cube.

        coordz: str
            name of depth coordinate

    Returns
    -------
    iris.cube.Cube
        collapsed cube.
    """
    # ####
    depth = cube.coord(coordz)
    thickness = depth.bounds[..., 1] - depth.bounds[..., 0]

    if depth.ndim == 1:
        slices = [None for i in cube.shape]
        coord_dim = cube.coord_dims(coordz)[0]
        slices[coord_dim] = slice(None)
        thickness = np.abs(thickness[tuple(slices)])

    ones = np.ones_like(cube.data)

    weights = thickness * ones

    result = cube.collapsed(coordz, iris.analysis.SUM,
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

    ie:
      extract_transect(cube, longitude=-28)
        will produce a transect along 28 West.

      extract_transect(cube, longitude=-28, latitude=[-50, 50])
        will produce a transect along 28 West  between 50 south and 50 North.

    This function is not yet implemented for irregular arrays - instead
    try the extract_trajectory function, but note that it is currently
    very slow. Alternatively, use the regrid preprocessor to regrid along
    a regular grid and then extract the transect.

    Arguments
    ---------
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
    """
    # ###
    coord_dim2 = False
    second_coord_range = False
    lats = cube.coord('latitude')
    lons = cube.coord('longitude')

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

    for dim_name, dim_cut, coord in zip(['latitude', 'longitude'],
                                        [latitude, longitude], [lats, lons]):
        # ####
        # Look for the first coordinate.
        if isinstance(dim_cut, float):
            coord_index = lats.nearest_neighbour_index(dim_cut)
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

    Arguments
    ---------
        cube: iris.cube.Cube
            input cube.

        latitudes: list of floats
            list of latitude coordinates.

        longitudes: list of floats
            list of longitude coordinates.

        number_points: int
            number of points to extrapolate (optional).

    Returns
    -------
    iris.cube.Cube
        collapsed cube.
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
