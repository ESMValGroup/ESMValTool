from __future__ import (absolute_import, division, print_function)
from six.moves import (filter, input, map, range, zip)  # noqa
import six

from glob import glob
import os

import iris
import iris.coords
import iris.coord_systems
import iris.cube
import iris.fileformats.pp
import numpy as np
import stratify


# Supported horizontal regridding schemes.
_schemes = dict(linear=iris.analysis.Linear(),
                nearest=iris.analysis.Nearest(),
                areaweighted=iris.analysis.AreaWeighted(),
                unstructurednearest=iris.analysis.UnstructuredNearest())


def _stock_1x1():
    """Generate a 1 degree by 1 degree global stock cube."""
    cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)

    ydata = np.linspace(-89.5, 89.5, 180)
    lats = iris.coords.DimCoord(ydata,
                                standard_name='latitude',
                                units='degrees_north',
                                coord_system=cs)

    xdata = np.linspace(0.5, 359.5, 360)
    lons = iris.coords.DimCoord(xdata,
                                standard_name='longitude',
                                units='degrees_east',
                                coord_system=cs)

    shape = (ydata.size, xdata.size)
    coords = [(lats, 0), (lons, 1)]
    cube = iris.cube.Cube(np.zeros(shape),
                          dim_coords_and_dims=coords)
    return cube


# Cache a stock of standard horizontal target grids.
# XXX: This is only for horizontal grids at the moment
#      we need to extend this concept to include standard
#      vertical pressure levels (perhaps).
_cache = {'1x1': _stock_1x1()}


def _vinterp(src_cube, tgt_grid):
    # Default to no interpolation and pass-thru the original source cube.
    result = src_cube

    # Get the target grid pressure coordinate, if available.
    tgt_plev = tgt_grid.coords('air_pressure')

    if tgt_plev:
        # Unpack the target pressure coordinate.
        tgt_plev, = tgt_plev
        tgt_plev_shape = tgt_plev.shape
        tgt_plev_units = tgt_plev.units

        # Get the source cube pressure coordinate, if available.
        src_plev = src_cube.coords('air_pressure')

        if src_plev:
            # Unpack the source pressure coordinate.
            src_plev, = src_plev
            # Get the source cube pressure coordinate points.
            src_plev_points = src_plev.points

            # Perform unit conversion if necessary, but CMOR checking
            # should have caught this hopefully!
            if src_plev.units != tgt_plev_units:
                src_plev_points = src_plev.units.convert(src_plev_points,
                                                         tgt_plev_units)

            if src_plev.shape == tgt_plev_shape and \
               np.allclose(src_plev_points, tgt_plev.points):
                # Don't perform vertical interpolation if the target and
                # source pressure coordinates are "similar" enough.
                return result

            # Determine the source axis for vertical interpolation.
            z_axis, = src_cube.coord_dims(src_plev)
            x_axis, = src_cube.coord_dims(src_cube.coord(axis='x'))
            y_axis, = src_cube.coord_dims(src_cube.coord(axis='y'))

            # Broadcast the 1d source pressure coordinate to fully describe
            # the spatial extent that will be interpolated.
            axes = (z_axis, y_axis, x_axis)
            min_axis, max_axis = min(axes), max(axes)
            reshape = [1] * (max_axis - min_axis + 1)
            broadcast_shape = [None] * len(reshape)
            for i in range(min_axis, max_axis+1):
                offset, N = i - min_axis, src_cube.shape[i]
                if i == z_axis:
                    reshape[offset] = N
                broadcast_shape[offset] = N

            src_plev_reshape = src_plev_points.reshape(reshape)
            src_plev_broadcast = np.broadcast_to(src_plev_reshape,
                                                 broadcast_shape)

            # Now perform the actual vertical interpolation.
            new_data = stratify.interpolate(tgt_plev.points,
                                            src_plev_broadcast,
                                            src_cube.data,
                                            axis=z_axis)

            result = tgt_grid.copy(data=new_data)

    return result


def regrid(src_cube, tgt_grid, scheme):
    if scheme is None:
        if tgt_grid is not None:
            emsg = 'Target grid must be None if no scheme is given.'
            raise ValueError(emsg)
        return src_cube

    if isinstance(scheme, six.string_types):
        scheme = scheme.lower()

    if _scheme.get(scheme) is None:
        emsg = 'Unknown regridding scheme, got {!r}.'
        raise ValueError(emsg.format(scheme))

    if tgt_grid is None:
        raise ValueError('Target grid must not be None.')
    elif tgt_grid == '1x1':
        tgt_grid = _cache[tgt_grid]
    elif not isinstance(tgt_grid, iris.cube.Cube):
        emsg = 'Expecting a {!r} instance or "1x1", got {!r}.'
        raise ValueError(emsg.format(iris.cube.Cube.__name__,
                                     type(tgt_grid)))

    # Perform horizontal regridding.
    cube = src_cube.regrid(tgt_grid, _schemes[scheme])

    # Perform vertical regridding.
    result_cube = _vinterp(cube, tgt_cube)

    return result_cube
