"""
Unit tests for the :mod:`esmvalcore.preprocessor.regrid` module.

"""

import iris
import numpy as np
from iris.coords import AuxCoord, CellMethod, DimCoord


def _make_vcoord(data, dtype=None):
    """
    Create a synthetic test vertical coordinate.

    """
    if dtype is None:
        dtype = np.dtype('int8')

    if isinstance(data, int):
        data = np.arange(data, dtype=dtype)
    elif not isinstance(data, np.ndarray):
        data = np.asarray(data, dtype=dtype)

    # Create a pressure vertical coordinate.
    kwargs = dict(
        standard_name='air_pressure',
        long_name='Pressure',
        var_name='plev',
        units='hPa',
        attributes=dict(positive='down'),
        coord_system=None)

    try:
        zcoord = DimCoord(data, **kwargs)
    except ValueError:
        zcoord = AuxCoord(data, **kwargs)

    return zcoord


def _make_cube(data, aux_coord=True, dim_coord=True, dtype=None):
    """
    Create a 3d synthetic test cube.

    """
    if dtype is None:
        dtype = np.dtype('int8')

    if not isinstance(data, np.ndarray):
        data = np.empty(data, dtype=dtype)

    z, y, x = data.shape

    # Create the cube.
    cm = CellMethod(
        method='mean', coords='time', intervals='20 minutes', comments=None)
    kwargs = dict(
        standard_name='air_temperature',
        long_name='Air Temperature',
        var_name='ta',
        units='K',
        attributes=dict(cube='attribute'),
        cell_methods=(cm, ))
    cube = iris.cube.Cube(data, **kwargs)

    # Create a synthetic test vertical coordinate.
    if dim_coord:
        cube.add_dim_coord(_make_vcoord(z, dtype=dtype), 0)

    # Create a synthetic test latitude coordinate.
    data = np.arange(y, dtype=dtype) + 1
    cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
    kwargs = dict(
        standard_name='latitude',
        long_name='Latitude',
        var_name='lat',
        units='degrees_north',
        attributes=dict(latitude='attribute'),
        coord_system=cs)
    ycoord = DimCoord(data, **kwargs)
    if data.size > 1:
        ycoord.guess_bounds()
    cube.add_dim_coord(ycoord, 1)

    # Create a synthetic test longitude coordinate.
    data = np.arange(x, dtype=dtype) + 1
    kwargs = dict(
        standard_name='longitude',
        long_name='Longitude',
        var_name='lon',
        units='degrees_east',
        attributes=dict(longitude='attribute'),
        coord_system=cs)
    xcoord = DimCoord(data, **kwargs)
    if data.size > 1:
        xcoord.guess_bounds()
    cube.add_dim_coord(xcoord, 2)

    # Create a synthetic test 2d auxiliary coordinate
    # that spans the vertical dimension.
    if aux_coord:
        data = np.arange(np.prod((z, y)), dtype=dtype).reshape(z, y)
        kwargs = dict(
            standard_name=None,
            long_name='Pressure Slice',
            var_name='aplev',
            units='hPa',
            attributes=dict(positive='down'),
            coord_system=None)
        zycoord = AuxCoord(data, **kwargs)
        cube.add_aux_coord(zycoord, (0, 1))

    return cube
