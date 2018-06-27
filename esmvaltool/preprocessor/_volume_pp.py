"""
Volume and z coordinate operations on data cubes.

Allows for selecting data subsets using certain volume bounds;
selecting depth or height regions; constructing volumetric averages;
"""
import iris


# slice cube over a restricted area (box)
def volume_slice(mycube, z_min, z_max):
    """
    Subset a cube on volume

    Function that subsets a cube on a box (z_min,z_max)
    This function is a restriction of masked_cube_lonlat();
    Returns a cube
    """
    subz = iris.Constraint(
        depth=lambda cell: float(z_min) <= cell <= float(z_max))

    region_subset = mycube.extract(subz)
    return region_subset


def volume_average(mycube, coordz, coord1, coord2,):
    """
    Determine the area average.

    Can be used with coord1 and coord2 (strings,
    usually 'longitude' and 'latitude' but depends on the cube);

    Returns a cube
    """
    # CMOR ised data should already have bounds?
    #    mycube.coord(coord1).guess_bounds()
    #    mycube.coord(coord2).guess_bounds()

    depth = mycube.coord(coordz)
    thickness = depth.bounds[:, 1] - depth.bounds[:, 0]  # 1D

    area = iris.analysis.cartography.area_weights(mycube)
    if area.ndim == 4 and thickness.ndim == 1:
        grid_volume = area * thickness[None, :, None, None]
    else:
        grid_volume = area * thickness
    result = mycube.collapsed(
        [coordz, coord1, coord2], iris.analysis.MEAN, weights=grid_volume)

    return result


# get the depth integration
def depth_integration(mycube, coordz):
    """
    Determine the total sum over the vertical component.

    Requires a 3D cube, and the name of the z coordinate.
    Returns a cube 2D.
    """
    ####
    depth = mycube.coord(coordz)
    thickness = depth.bounds[..., 1] - depth.bounds[..., 0]  # 1D

    if depth.ndim == 1:
        slices = [None for i in mycube.shape]
        coord_dim = mycube.coord_dims(coordz)[0]
        slices[coord_dim] = slice(None)
        thickness = thickness[tuple(slices)]

    result = mycube * thickness
    result = mycube.collapsed(
        [
            coordz,
        ],
        iris.analysis.SUM,
    )
    result.rename('Depth_integrated_' + str(mycube.name()))

    # result.units = Unit('m') * result.units # This doesn't work:
    # TODO: Change units on cube to reflect 2D concentration (not 3D)
    # Waiting for news from iris community.
    return result


def extract_transect(mycube, latitude=None, longitude=None):
    """
    Extract data along a line of constant latitude or longitude.

    A range may also be extracted using a minimum and maximum
    value for latitude or longitude.

    This function is not yet implemented for irregular arrays - instead
    try the extract_trajectory function, but note that it is currently
    very slow.
    """
    ####
    coord_dim2 = False
    second_coord_range = False
    lats = mycube.coord('latitude')
    lons = mycube.coord('longitude')

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
            coord_dim = mycube.coord_dims(dim_name)[0]

        #####
        # Look for the second coordinate.
        if isinstance(dim_cut, list):
            coord_dim2 = mycube.coord_dims(dim_name)[0]
            second_coord_range = [coord.nearest_neighbour_index(dim_cut[0]),
                                  coord.nearest_neighbour_index(dim_cut[1])]
    #####
    # Extracting the line of constant longitude/latitude
    slices = [slice(None) for i in mycube.shape]
    slices[coord_dim] = coord_index

    if second_coord_range:
        slices[coord_dim2] = slice(second_coord_range[0],
                                   second_coord_range[1])
    return mycube[tuple(slices)]


# extract along a trajectory
def extract_trajectory(mycube, latitudes, longitudes, number_points):
    """
    Extract data along a trajectory.

    latitudes and longitudes are the pairs of coordinates for two points.
    number_points is the number of points between the two points.

    This version uses the expensive interpolate method, but it may be
    necceasiry for irregular grids.
    """
    from iris.analysis.trajectory import interpolate
    import numpy as np

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
    interpolated_cube = interpolate(mycube, points)  # Very slow!
    return interpolated_cube
