"""Tests for the module :mod:`esmvaltool.cmorizers.obs.utilities`."""

import dask.array as da
import iris
import numpy as np
import pytest

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
