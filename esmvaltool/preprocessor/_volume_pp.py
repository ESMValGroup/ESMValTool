"""
Volume and z coordinate operations on data cubes.

Allows for selecting data subsets using certain volume bounds;
selecting depth or height regions; constructing volumetric averages;
"""
import iris
import numpy as np


# slice cube over a restricted area (box)
def volume_slice(cube, z_min, z_max):
    """
    Subset a cube on volume

    Function that subsets a cube on a box (z_min,z_max)
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
    # CMOR ised data should already have bounds?
    #    cube.coord(coord1).guess_bounds()
    #    cube.coord(coord2).guess_bounds()
    depth = cube.coord(coordz)
    thickness = depth.bounds[..., 1] - depth.bounds[..., 0]

    area = iris.analysis.cartography.area_weights(cube)

    if depth.ndim == 1:
        slices = [None for i in cube.shape]
        coord_dim = cube.coord_dims(coordz)[0]
        slices[coord_dim] = slice(None)
        thickness = np.abs(thickness[tuple(slices)])

    grid_volume = area * thickness

    result = cube.collapsed(
        [coordz, coord1, coord2], iris.analysis.MEAN, weights=grid_volume)

    return result


# get the depth integration
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
    ####
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

      extract_transect(cube, longitude=-28, latitude=[-50,50])
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
    ####
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
        #####
        # Look for the first coordinate.
        if isinstance(dim_cut, float):
            coord_index = lats.nearest_neighbour_index(dim_cut)
            coord_dim = cube.coord_dims(dim_name)[0]

        #####
        # Look for the second coordinate.
        if isinstance(dim_cut, list):
            coord_dim2 = cube.coord_dims(dim_name)[0]
            second_coord_range = [coord.nearest_neighbour_index(dim_cut[0]),
                                  coord.nearest_neighbour_index(dim_cut[1])]
    #####
    # Extracting the line of constant longitude/latitude
    slices = [slice(None) for i in cube.shape]
    slices[coord_dim] = coord_index

    if second_coord_range:
        slices[coord_dim2] = slice(second_coord_range[0],
                                   second_coord_range[1])
    return cube[tuple(slices)]


# extract along a trajectory
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
