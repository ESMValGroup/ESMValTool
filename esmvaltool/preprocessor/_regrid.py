"""Horizontal and vertical regridding module."""

import os
import re
from copy import deepcopy

import iris
import numpy as np
import six
import stratify
from iris.analysis import AreaWeighted, Linear, Nearest, UnstructuredNearest

from ._io import concatenate_callback, load
from ._regrid_esmpy import ESMF_REGRID_METHODS
from ._regrid_esmpy import regrid as esmpy_regrid
from ..cmor.fix import fix_file, fix_metadata
from ..cmor.table import CMOR_TABLES

# Regular expression to parse a "MxN" cell-specification.
_CELL_SPEC = re.compile(
    r'''\A
        \s*(?P<dlon>\d+(\.\d+)?)\s*
        x
        \s*(?P<dlat>\d+(\.\d+)?)\s*
        \Z
     ''', re.IGNORECASE | re.VERBOSE)

# Default fill-value.
_MDI = 1e+20

# Stock cube - global grid extents (degrees).
_LAT_MIN = -90.0
_LAT_MAX = 90.0
_LAT_RANGE = _LAT_MAX - _LAT_MIN
_LON_MIN = 0.0
_LON_MAX = 360.0
_LON_RANGE = _LON_MAX - _LON_MIN

# A cached stock of standard horizontal target grids.
_CACHE = dict()

# Supported horizontal regridding schemes.
HORIZONTAL_SCHEMES = {
    'linear': Linear(extrapolation_mode='mask'),
    'linear_extrapolate': Linear(extrapolation_mode='extrapolate'),
    'nearest': Nearest(extrapolation_mode='mask'),
    'area_weighted': AreaWeighted(),
    'unstructured_nearest': UnstructuredNearest(),
}

# Supported vertical interpolation schemes.
VERTICAL_SCHEMES = ('linear', 'nearest',
                    'linear_horizontal_extrapolate_vertical',
                    'nearest_horizontal_extrapolate_vertical')


def parse_cell_spec(spec):
    """
    Parse an MxN cell specification string.

    Parameters
    ----------
    spec: str

    Returns
    -------
    tuple
        tuple of (float, float) of parsed (lon, lat)

    Raises
    ------
    ValueError
        if the MxN cell specification is malformed.
    ValueError
        invalid longitude and latitude delta in cell specification.
    """
    cell_match = _CELL_SPEC.match(spec)
    if cell_match is None:
        emsg = 'Invalid MxN cell specification for grid, got {!r}.'
        raise ValueError(emsg.format(spec))

    cell_group = cell_match.groupdict()
    dlon = float(cell_group['dlon'])
    dlat = float(cell_group['dlat'])

    if (np.trunc(_LON_RANGE / dlon) * dlon) != _LON_RANGE:
        emsg = ('Invalid longitude delta in MxN cell specification '
                'for grid, got {!r}.')
        raise ValueError(emsg.format(dlon))

    if (np.trunc(_LAT_RANGE / dlat) * dlat) != _LAT_RANGE:
        emsg = ('Invalid latitude delta in MxN cell specification '
                'for grid, got {!r}.')
        raise ValueError(emsg.format(dlat))

    return dlon, dlat


def _stock_cube(spec, lat_offset=True, lon_offset=True):
    """
    Create a stock cube.

    Create a global cube with M degree-east by N degree-north regular grid
    cells.

    The longitude range is from 0 to 360 degrees. The latitude range is from
    -90 to 90 degrees. Each cell grid point is calculated as the mid-point of
    the associated MxN cell.

    Parameters
    ----------
    spec : str
        Specifies the 'MxN' degree cell-specification for the global grid.
    lat_offset : bool
        Offset the grid centers of the latitude coordinate w.r.t. the
        pole by half a grid step. This argument is ignored if `target_grid`
        is a cube or file.
    lon_offset : bool
        Offset the grid centers of the longitude coordinate w.r.t. Greenwich
        meridian by half a grid step.
        This argument is ignored if `target_grid` is a cube or file.

    Returns
    -------
        A :class:`~iris.cube.Cube`.

    """
    dlon, dlat = parse_cell_spec(spec)
    mid_dlon, mid_dlat = dlon / 2, dlat / 2

    # Construct the latitude coordinate, with bounds.
    if lat_offset:
        latdata = np.linspace(_LAT_MIN + mid_dlat, _LAT_MAX - mid_dlat,
                              _LAT_RANGE / dlat)
    else:
        latdata = np.linspace(_LAT_MIN, _LAT_MAX, _LAT_RANGE / dlat + 1)

    # Construct the longitude coordinat, with bounds.
    if lon_offset:
        londata = np.linspace(_LON_MIN + mid_dlon, _LON_MAX - mid_dlon,
                              _LON_RANGE / dlon)
    else:
        londata = np.linspace(_LON_MIN, _LON_MAX - dlon, _LON_RANGE / dlon)

    lats = iris.coords.DimCoord(
        latdata,
        standard_name='latitude',
        units='degrees_north',
        var_name='lat')
    lats.guess_bounds()

    lons = iris.coords.DimCoord(
        londata,
        standard_name='longitude',
        units='degrees_east',
        var_name='lon')
    lons.guess_bounds()

    # Construct the resultant stock cube, with dummy data.
    shape = (latdata.size, londata.size)
    dummy = np.empty(shape, dtype=np.dtype('int8'))
    coords_spec = [(lats, 0), (lons, 1)]
    cube = iris.cube.Cube(dummy, dim_coords_and_dims=coords_spec)

    return cube


def _attempt_irregular_regridding(cube, scheme):
    """Check if irregular regridding with ESMF should be used."""
    if scheme in ESMF_REGRID_METHODS:
        try:
            lat_dim = cube.coord('latitude').ndim
            lon_dim = cube.coord('longitude').ndim
            if lat_dim == lon_dim == 2:
                return True
        except iris.exceptions.CoordinateNotFoundError:
            pass
    return False


def regrid(cube, target_grid, scheme, lat_offset=True, lon_offset=True):
    """
    Perform horizontal regridding.

    Parameters
    ----------
    cube : cube
        The source cube to be regridded.
    target_grid : cube or str
        The cube that specifies the target or reference grid for the regridding
        operation. Alternatively, a string cell specification may be provided,
        of the form 'MxN', which specifies the extent of the cell, longitude by
        latitude (degrees) for a global, regular target grid.
    scheme : str
        The regridding scheme to perform, choose from
        'linear',
        'linear_extrapolate',
        'nearest',
        'area_weighted',
        'unstructured_nearest'.
    lat_offset : bool
        Offset the grid centers of the latitude coordinate w.r.t. the
        pole by half a grid step. This argument is ignored if `target_grid`
        is a cube or file.
    lon_offset : bool
        Offset the grid centers of the longitude coordinate w.r.t. Greenwich
        meridian by half a grid step.
        This argument is ignored if `target_grid` is a cube or file.

    Returns
    -------
    cube

    See Also
    --------
    extract_levels : Perform vertical regridding.

    """
    if HORIZONTAL_SCHEMES.get(scheme.lower()) is None:
        emsg = 'Unknown regridding scheme, got {!r}.'
        raise ValueError(emsg.format(scheme))

    if isinstance(target_grid, six.string_types):
        if os.path.isfile(target_grid):
            target_grid = iris.load_cube(target_grid)
        else:
            # Generate a target grid from the provided cell-specification,
            # and cache the resulting stock cube for later use.
            target_grid = _CACHE.setdefault(
                target_grid,
                _stock_cube(target_grid, lat_offset, lon_offset),
            )
            # Align the target grid coordinate system to the source
            # coordinate system.
            src_cs = cube.coord_system()
            xcoord = target_grid.coord(axis='x', dim_coords=True)
            ycoord = target_grid.coord(axis='y', dim_coords=True)
            xcoord.coord_system = src_cs
            ycoord.coord_system = src_cs

    if not isinstance(target_grid, iris.cube.Cube):
        raise ValueError('Expecting a cube, got {}.'.format(target_grid))

    # Unstructured regridding requires x2 2d spatial coordinates,
    # so ensure to purge any 1d native spatial dimension coordinates
    # for the regridder.
    if scheme == 'unstructured_nearest':
        for axis in ['x', 'y']:
            coords = cube.coords(axis=axis, dim_coords=True)
            if coords:
                [coord] = coords
                cube.remove_coord(coord)

    # Perform the horizontal regridding.
    if _attempt_irregular_regridding(cube, scheme):
        cube = esmpy_regrid(cube, target_grid, scheme)
    else:
        cube = cube.regrid(target_grid, HORIZONTAL_SCHEMES[scheme])

    return cube


def _create_cube(src_cube, data, levels):
    """
    Generate a new cube with the interpolated data.

    The resultant cube is seeded with `src_cube` metadata and coordinates,
    excluding any source coordinates that span the associated vertical
    dimension. The `levels` of interpolation are used along with the
    associated source cube vertical coordinate metadata to add a new
    vertical coordinate to the resultant cube.

    Parameters
    ----------
    src_cube : cube
        The source cube that was vertically interpolated.
    data : array
        The payload resulting from interpolating the source cube
        over the specified levels.
    levels : array
        The vertical levels of interpolation.

    Returns
    -------
    cube

    .. note::

        If there is only one level of interpolation, the resultant cube
        will be collapsed over the associated vertical dimension, and a
        scalar vertical coordinate will be added.

    """
    # Get the source cube vertical coordinate and associated dimension.
    src_levels = src_cube.coord(axis='z', dim_coords=True)
    z_dim, = src_cube.coord_dims(src_levels)

    if data.shape[z_dim] != levels.size:
        emsg = ('Mismatch between data and levels for data dimension {!r}, '
                'got data shape {!r} with levels shape {!r}.')
        raise ValueError(emsg.format(z_dim, data.shape, levels.shape))

    # Construct the resultant cube with the interpolated data
    # and the source cube metadata.
    kwargs = deepcopy(src_cube.metadata)._asdict()
    result = iris.cube.Cube(data, **kwargs)

    # Add the appropriate coordinates to the cube, excluding
    # any coordinates that span the z-dimension of interpolation.
    for coord in src_cube.dim_coords:
        [dim] = src_cube.coord_dims(coord)
        if dim != z_dim:
            result.add_dim_coord(coord.copy(), dim)

    for coord in src_cube.aux_coords:
        dims = src_cube.coord_dims(coord)
        if z_dim not in dims:
            result.add_aux_coord(coord.copy(), dims)

    for coord in src_cube.derived_coords:
        dims = src_cube.coord_dims(coord)
        if z_dim not in dims:
            result.add_aux_coord(coord.copy(), dims)

    # Construct the new vertical coordinate for the interpolated
    # z-dimension, using the associated source coordinate metadata.
    kwargs = deepcopy(src_levels._as_defn())._asdict()

    try:
        coord = iris.coords.DimCoord(levels, **kwargs)
        result.add_dim_coord(coord, z_dim)
    except ValueError:
        coord = iris.coords.AuxCoord(levels, **kwargs)
        result.add_aux_coord(coord, z_dim)

    # Collapse the z-dimension for the scalar case.
    if levels.size == 1:
        slicer = [slice(None)] * result.ndim
        slicer[z_dim] = 0
        result = result[tuple(slicer)]

    return result


def _vertical_interpolate(cube, levels, interpolation, extrapolation):
    """Perform vertical interpolation."""
    # Determine the source levels and axis for vertical interpolation.
    src_levels = cube.coord(axis='z', dim_coords=True)
    z_axis, = cube.coord_dims(src_levels)

    # Broadcast the 1d source cube vertical coordinate to fully
    # describe the spatial extent that will be interpolated.
    broadcast_shape = cube.shape[z_axis:]
    reshape = [1] * len(broadcast_shape)
    reshape[0] = cube.shape[z_axis]
    src_levels_reshaped = src_levels.points.reshape(reshape)
    src_levels_broadcast = np.broadcast_to(src_levels_reshaped,
                                           broadcast_shape)

    # force mask onto data as nan's
    if np.ma.is_masked(cube.data):
        cube.data[cube.data.mask] = np.nan

    # Now perform the actual vertical interpolation.
    new_data = stratify.interpolate(
        levels,
        src_levels_broadcast,
        cube.data,
        axis=z_axis,
        interpolation=interpolation,
        extrapolation=extrapolation)

    # Calculate the mask based on the any NaN values in the interpolated data.
    mask = np.isnan(new_data)

    if np.any(mask):
        # Ensure that the data is masked appropriately.
        new_data = np.ma.array(new_data, mask=mask, fill_value=_MDI)

    # Construct the resulting cube with the interpolated data.
    return _create_cube(cube, new_data, levels.astype(float))


def extract_levels(cube, levels, scheme):
    """
    Perform vertical interpolation.

    Parameters
    ----------
    cube : cube
        The source cube to be vertically interpolated.
    levels : array
        One or more target levels for the vertical interpolation. Assumed
        to be in the same S.I. units of the source cube vertical dimension
        coordinate.
    scheme : str
        The vertical interpolation scheme to use. Choose from
        'linear',
        'nearest',
        'nearest_horizontal_extrapolate_vertical',
        'linear_horizontal_extrapolate_vertical'.

    Returns
    -------
    cube

    See Also
    --------
    regrid : Perform horizontal regridding.

    """
    if scheme not in VERTICAL_SCHEMES:
        emsg = 'Unknown vertical interpolation scheme, got {!r}. '
        emsg += 'Possible schemes: {!r}'
        raise ValueError(emsg.format(scheme, VERTICAL_SCHEMES))

    # This allows us to put level 0. to load the ocean surface.
    extrap_scheme = 'nan'
    if scheme == 'nearest_horizontal_extrapolate_vertical':
        scheme = 'nearest'
        extrap_scheme = 'nearest'

    if scheme == 'linear_horizontal_extrapolate_vertical':
        scheme = 'linear'
        extrap_scheme = 'nearest'

    # Ensure we have a non-scalar array of levels.
    levels = np.array(levels, ndmin=1)

    # Get the source cube vertical coordinate, if available.
    src_levels = cube.coord(axis='z', dim_coords=True)

    if (src_levels.shape == levels.shape
            and np.allclose(src_levels.points, levels)):
        # Only perform vertical extraction/interploation if the source
        # and target levels are not "similar" enough.
        result = cube
    elif set(levels).issubset(set(src_levels.points)):
        # If all target levels exist in the source cube, simply extract them.
        name = src_levels.name()
        coord_values = {name: lambda cell: cell.point in set(levels)}
        constraint = iris.Constraint(coord_values=coord_values)
        result = cube.extract(constraint)
        # Ensure the constraint did not fail.
        if not result:
            emsg = 'Failed to extract levels {!r} from cube {!r}.'
            raise ValueError(emsg.format(list(levels), name))
    else:
        # As a last resort, perform vertical interpolation.
        result = _vertical_interpolate(cube, levels, scheme, extrap_scheme)

    return result


def get_cmor_levels(cmor_table, coordinate):
    """Get level definition from a CMOR coordinate.

    Parameters
    ----------
    cmor_table: str
        CMOR table name
    coordinate: str
        CMOR coordinate name

    Returns
    -------
    list[int]

    Raises
    ------
    ValueError:
        If the CMOR table is not defined, the coordinate does not specify any
        levels or the string is badly formatted.

    """
    if cmor_table not in CMOR_TABLES:
        raise ValueError(
            "Level definition cmor_table '{}' not available".format(
                cmor_table))

    if coordinate not in CMOR_TABLES[cmor_table].coords:
        raise ValueError('Coordinate {} not available for {}'.format(
            coordinate, cmor_table))

    cmor = CMOR_TABLES[cmor_table].coords[coordinate]

    if cmor.requested:
        return [float(level) for level in cmor.requested]
    if cmor.value:
        return [float(cmor.value)]

    raise ValueError(
        'Coordinate {} in {} does not have requested values'.format(
            coordinate, cmor_table))


def get_reference_levels(filename,
                         project,
                         dataset,
                         short_name,
                         fix_dir):
    """Get level definition from a CMOR coordinate.

    Parameters
    ----------
    filename: str
        Path to the reference file

    Returns
    -------
    list[float]

    Raises
    ------
    ValueError:
        If the dataset is not defined, the coordinate does not specify any
        levels or the string is badly formatted.

    """
    filename = fix_file(filename, short_name, project, dataset, fix_dir)
    cubes = load(filename, callback=concatenate_callback)
    cubes = fix_metadata(cubes, short_name, project, dataset)
    cube = cubes[0]
    try:
        coord = cube.coord(axis='Z')
    except iris.exceptions.CoordinateNotFoundError:
        raise ValueError('z-coord not available in {}'.format(filename))
    return coord.points.tolist()
