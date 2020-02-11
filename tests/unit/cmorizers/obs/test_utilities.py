"""Tests for the module :mod:`esmvaltool.cmorizers.obs.utilities`."""

import dask.array as da
import iris
import numpy as np
import pytest

from cf_units import Unit
import esmvaltool.cmorizers.obs.utilities as utils


def np_to_da(array, lazy):
    if not lazy:
        return array
    if array is None:
        return array
    return da.from_array(array)


def is_lazy(cube):
    if not cube.has_lazy_data():
        return False
    for coord in cube.coords(dim_coords=False):
        if not coord.has_lazy_points():
            return False
        if coord.has_bounds():
            if not coord.has_lazy_bounds():
                return False
    return True


def cubes_generator(lazy=True):
    cube_datas = [
        np.array([[0, 1], [-1, 0]], dtype=np.int),
        np.array([[0.0, 1.0], [-1.0, 0.0]], dtype=np.float32),
        np.array([[0.0, 1.0], [-1.0, 0.0]], dtype=np.float64),
        np.ma.masked_equal([[0, 1], [2, 3]], 3).astype(np.int),
        np.ma.masked_values([[0.0, 1.0], [2.0, 3.0]], 3.0).astype(np.float32),
        np.ma.masked_values([[0.0, 1.0], [2.0, 3.0]], 3.0).astype(np.float64),
    ]
    x_coords = [
        (np.array([1, 3], dtype=np.int), None),
        (np.array([1, 3],
                  dtype=np.int), np.array([[0, 2], [2, 4]], dtype=np.int)),
        (np.array([1.0, 3.0], dtype=np.float32),
         np.array([[0.0, 2.0], [2.0, 4.0]], dtype=np.float32)),
        (np.array([1.0, 3.0], dtype=np.float64), None),
        (np.array([1.0, 3.0], dtype=np.float64),
         np.array([[0.0, 2.0], [2.0, 4.0]], dtype=np.float64)),
    ]
    y_coords = [
        (np.array([1, 3], dtype=np.int),
         np.array([[0.0, 2.0], [2.0, 4.0]], dtype=np.float32)),
        (np.array([1.0, 3.0], dtype=np.float32),
         np.array([[0.0, 2.0], [2.0, 4.0]], dtype=np.float64)),
        (np.array([1.0, 3.0],
                  dtype=np.float64), np.array([[0, 2], [2, 4]], dtype=np.int)),
    ]
    for cube_data in cube_datas:
        cube_data = np_to_da(cube_data, lazy)
        for x_val in x_coords:
            x_val = (np_to_da(x_val[0], lazy), np_to_da(x_val[1], lazy))
            x_coord = iris.coords.DimCoord(x_val[0],
                                           bounds=x_val[1],
                                           var_name='x')
            for y_val in y_coords:
                y_val = (np_to_da(y_val[0], lazy), np_to_da(y_val[1], lazy))
                y_coord = iris.coords.DimCoord(y_val[0],
                                               bounds=y_val[1],
                                               var_name='y')
                aux_coord = iris.coords.AuxCoord(y_val[0],
                                                 bounds=y_val[1],
                                                 var_name='aux')
                cube = iris.cube.Cube(
                    cube_data,
                    var_name='test_var',
                    dim_coords_and_dims=[(x_coord, 0), (y_coord, 1)],
                    aux_coords_and_dims=[(aux_coord, 0)],
                )
                yield cube


@pytest.mark.parametrize('cube', cubes_generator(lazy=True))
def test_fix_dtype_lazy(cube):
    assert is_lazy(cube)
    utils._fix_dtype(cube)
    assert cube.dtype == np.float32
    for coord in cube.coords():
        assert coord.dtype == np.float64
        if coord.has_bounds():
            assert coord.bounds_dtype == np.float64
    assert is_lazy(cube)


@pytest.mark.parametrize('cube', cubes_generator(lazy=False))
def test_fix_dtype_not_lazy(cube):
    assert not is_lazy(cube)
    utils._fix_dtype(cube)
    assert cube.dtype == np.float32
    for coord in cube.coords():
        assert coord.dtype == np.float64
        if coord.has_bounds():
            assert coord.bounds_dtype == np.float64
    assert not is_lazy(cube)


def _create_sample_cube():
    """Create a quick CMOR-compliant sample cube."""
    coord_sys = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
    cube_data = np.ones((2, 3, 2, 2))
    cube_data[1, 1, 1, 1] = 22.
    time = iris.coords.DimCoord([15, 45],
                                standard_name='time',
                                bounds=[[1., 30.], [30., 60.]],
                                units=Unit(
                                    'days since 1950-01-01',
                                    calendar='gregorian'))
    zcoord = iris.coords.DimCoord([0.5, 5., 50.],
                                  var_name='lev',
                                  standard_name='depth',
                                  bounds=[[0., 2.5], [2.5, 25.],
                                          [25., 250.]],
                                  units='m',
                                  attributes={'positive': 'down'})
    lons = iris.coords.DimCoord([1.5, 2.5],
                                standard_name='longitude',
                                bounds=[[1., 2.], [2., 3.]],
                                units='degrees_east',
                                coord_system=coord_sys)
    lats = iris.coords.DimCoord([1.5, 2.5],
                                standard_name='latitude',
                                bounds=[[1., 2.], [2., 3.]],
                                units='degrees_north',
                                coord_system=coord_sys)
    coords_spec = [(time, 0), (zcoord, 1), (lats, 2), (lons, 3)]
    cube = iris.cube.Cube(cube_data, dim_coords_and_dims=coords_spec)
    return cube


def test_add_scalar_height_coord():
    """Test add height aux coord."""
    cube = _create_sample_cube()
    utils.add_scalar_height_coord(cube, height=10.)
    assert cube.coord("height").points[0] == 10.
    assert "positive" in cube.coord("height").attributes
    assert cube.coord("height").attributes["positive"] == "up"


def test_convert_time_units():
    """Test convert time units functionlity."""
    cube = _create_sample_cube()
    cube.coord("time").units = 'months since 0000-01-01 00:00:00'
    utils.convert_timeunits(cube, "1950")
    converted_units = cube.coord("time").units
    assert converted_units == 'months since 1950-01-01 00:00:00'
    cube.coord("time").units = 'days since 0000-01-01 00:00:00'
    utils.convert_timeunits(cube, "1950")
    converted_units = cube.coord("time").units
    assert converted_units == 'days since 1950-01-01 00:00:00'
    cube.coord("time").units = 'days since 1950-1-1'
    utils.convert_timeunits(cube, "1950")
    converted_units = cube.coord("time").units
    assert converted_units == 'days since 1950-1-1 00:00:00'
    cube = _create_sample_cube()
    utils.convert_timeunits(cube, "1950")
    not_converted_units = cube.coord("time").units
    assert not_converted_units == 'days since 1950-01-01 00:00:00'


def test_fix_coords():
    """Test fix coordinates."""
    cube = _create_sample_cube()
    cube.coord("time").bounds = None
    cube.coord('time').convert_units(
        Unit('days since 1850-1-1 00:00:00', calendar='gregorian'))
    cube.coord("longitude").bounds = None
    cube.coord("latitude").bounds = None
    cube.coord("depth").bounds = None
    cube.coord("longitude").points = cube.coord("longitude").points - 3.
    utils.fix_coords(cube)
    assert cube.coord("time").has_bounds()
    assert cube.coord("time").bounds[0][1] == 30.
    assert cube.coord("time").units == 'days since 1950-1-1 00:00:00'
    assert cube.coord("time").units.calendar == "gregorian"
    assert cube.coord("longitude").points[0] == 178.5
    assert cube.coord("longitude").points[1] == 179.5
    assert cube.coord("longitude").has_bounds()
    assert cube.coord("longitude").bounds[1][1] == 180.
    assert cube.data[1, 1, 1, 0] == 22.
    assert cube.coord("latitude").has_bounds()
    assert cube.coord("depth").has_bounds()







