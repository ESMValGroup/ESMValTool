"""
A package for performing horizontal regridding, and vertical level extraction
or vertical level interpolation.

"""

from __future__ import (absolute_import, division, print_function)
from six.moves import (filter, input, map, range, zip)  # noqa
import six

from copy import deepcopy
from glob import glob
import os
import re

import iris
import numpy as np
import stratify


# Regular expression to parse a "MxN" cell-specification.
_CELL_SPEC = re.compile(r'''\A
                            \s*(?P<dx>\d+(\.\d+)?)\s*
                            x
                            \s*(?P<dy>\d+(\.\d+)?)\s*
                            \Z
                         ''', re.IGNORECASE | re.VERBOSE)


# Default fill-value.
_MDI = 1e+20


# A cached stock of standard horizontal target grids.
_cache = dict()


# Supported horizontal regridding schemes.
schemes = dict(linear=iris.analysis.Linear(),
               nearest=iris.analysis.Nearest(),
               areaweighted=iris.analysis.AreaWeighted(),
               unstructurednearest=iris.analysis.UnstructuredNearest())


def _stock_cube(spec):
    """
    Create a global cube with M degree-east by N degree-north regular grid
    cells.

    The longitude range is from 0 to 360 degrees. The latitude range is from
    -90 to 90 degrees. Each cell grid point is calculated as the mid-point of
    the associated MxN cell.

    Paramaters
    ----------
    spec : str
        Specifies the 'MxN' degree cell-specification for the global grid.

    Returns
    -------
        A :class:`~iris.cube.Cube`.

    """
    # Parse the MxN cell specification string.
    cell_match = _CELL_SPEC.match(spec)
    if cell_match is None:
        emsg = 'Invalid MxN cell specification for stock cube, got {!r}.'
        raise ValueError(emsg.format(spec))

    cell_group = cell_match.groupdict()
    dx = float(cell_group['dx'])
    dy = float(cell_group['dy'])
    mid_dx, mid_dy = dx / 2, dy / 2

    cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)

    ydata = np.linspace(-90.0 + mid_dy, 90.0 - mid_dy, 180.0 / dy)
    lats = iris.coords.DimCoord(ydata,
                                standard_name='latitude',
                                units='degrees_north',
                                coord_system=cs)
    lats.guess_bounds()

    xdata = np.linspace(mid_dx, 360.0 - mid_dx, 360.0 / dx)
    lons = iris.coords.DimCoord(xdata,
                                standard_name='longitude',
                                units='degrees_east',
                                coord_system=cs)
    lons.guess_bounds()

    shape = (ydata.size, xdata.size)
    dummy = np.empty(shape, dtype=np.dtype('int8'))
    coords = [(lats, 0), (lons, 1)]
    cube = iris.cube.Cube(dummy, dim_coords_and_dims=coords)

    return cube


def regrid(src_cube, tgt_grid, scheme):
    """
    Perform horizontal regridding.

    Parameters
    ----------
    src_cube : cube
        The source cube to be regridded.
    tgt_cube : cube or str
        The cube that specifies the target or reference grid for the regridding
        operation. Alternatively, a string cell specification may be provided,
        of the form 'MxN', which specifies the extent of the cell, longitude by
        latitude (degrees) for a global, regular target grid.
    scheme : str
        The regridding scheme to perform, see `regrid.schemes`.

    Returns
    -------
    cube

    See Also
    --------
    vinterp : Perform vertical regridding.

    """
    if tgt_grid is None and scheme is None:
        # nop
        return src_cube

    if tgt_grid is None:
        emsg = 'A target grid must be specified for horizontal regridding.'
        raise ValueError(emsg)

    if scheme is None:
        emsg = 'A scheme must be specified for horizontal regridding.'
        raise ValueError(emsg)

    if schemes.get(scheme.lower()) is None:
        emsg = 'Unknown regridding scheme, got {!r}.'
        raise ValueError(emsg.format(scheme))

    if isinstance(tgt_grid, six.string_types):
        # Generate a target grid from the provided cell-specification,
        # and cache the resulting stock cube for later use.
        tgt_grid = _cache.setdefault(tgt_grid, _stock_cube(tgt_grid))
    elif not isinstance(tgt_grid, iris.cube.Cube):
        emsg = 'Expecting a cube or cell-specification, got {!r}.'
        raise ValueError(emsg.format(tgt_grid))

    # Perform the horizontal regridding.
    result = src_cube.regrid(tgt_grid, schemes[scheme])

    return result


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
    if coord.shape == (1,):
        slicer = [slice(None)] * result.ndim
        slicer[z_dim] = 0
        result = result[tuple(slicer)]

    return result


def vinterp(src_cube, levels):
    """
    Perform vertical interpolation.

    Paramaters
    ----------
    src_cube : cube
        The source cube to be vertically interpolated.
    levels : array
        One or more target levels for the vertical interpolation. Assumed
        to be in the same S.I. units of the source cube vertical dimension
        coordinate.

    Returns
    -------
    cube

    See Also
    --------
    regrid : Perform horizontal regridding.

    """
    # Default to no interpolation and pass-thru the original source cube.
    result = src_cube

    if levels is not None:
        # Ensure we have a non-scalar array of levels.
        levels = np.array(levels, ndmin=1)

        # Get the source cube vertical coordinate, if available.
        src_levels = src_cube.coord(axis='z', dim_coords=True)

        # Don't perform vertial interploation if the source and
        # target values are "similar" enough.
        if src_levels.shape == levels.shape and \
           np.allclose(src_levels.points, levels):
            return result

        # Determine whether we can simply extract the target levels,
        # if they *all* exist in the source cube, otherwise
        # perform vertical interpolation.
        if set(levels).issubset(set(src_levels.points)):
            name = src_levels.name()
            coord_values = {name: lambda cell: cell.point in set(levels)}
            constraint = iris.Constraint(coord_values=coord_values)
            result = src_cube.extract(constraint)
        else:
            # Determine the source axis for vertical interpolation.
            z_axis, = src_cube.coord_dims(src_levels)

            # Broadcast the 1d source cube vertical coordinate to fully
            # describe the spatial extent that will be interpolated.
            broadcast_shape = src_cube.shape[z_axis:]
            reshape = [1] * len(broadcast_shape)
            reshape[0] = src_cube.shape[z_axis]
            src_levels_reshaped = src_levels.points.reshape(reshape)
            src_levels_broadcast = np.broadcast_to(src_levels_reshaped,
                                                   broadcast_shape)

            # Now perform the actual vertical interpolation.
            new_data = stratify.interpolate(levels,
                                            src_levels_broadcast,
                                            src_cube.data,
                                            axis=z_axis,
                                            interpolation='linear',
                                            extrapolation='nan')

            # Determine whether we need to fill any extrapolated NaN values.
            mask = np.isnan(new_data)

            if np.any(mask):
                # Check the fill-value is appropriate for the result dtype.
                try:
                    [fill_value] = np.array([_MDI], dtype=new_data.dtype)
                except OverflowError:
                    emsg = 'Fill value of {!r} invalid for result {!r}.'
                    raise ValueError(emsg.format(_MDI, new_data.dtype))

                # Replace the NaN values with the fill-value.
                new_data[mask] = fill_value

            # Construct the resulting cube with the interpolated data.
            result = _create_cube(src_cube, new_data, levels)

    return result
