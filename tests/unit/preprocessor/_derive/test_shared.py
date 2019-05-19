"""Tests for the shared functions of the derive preprocessor."""
import copy

import numpy as np
import pytest
import iris
from iris.cube import CubeList
from cf_units import Unit

import esmvaltool.preprocessor._derive._shared as shared

O_NAME = 'sea_surface_temperature'
L_NAME = 'air_temperature'
SFTOF_CUBE = iris.cube.Cube(
    [100.0, 0.0, 50.0, 70.0],
    var_name='sftof',
    standard_name='sea_area_fraction',
    units=Unit('%'),
)
SFTLF_CUBE = iris.cube.Cube(
    [10.0, 0.0, 100.0],
    var_name='sftlf',
    standard_name='land_area_fraction',
    units=Unit('%'),
)
O_CUBE_1 = iris.cube.Cube(
    [1.0, 2.0, -1.0, 2.0],
    standard_name=O_NAME,
)
O_CUBE_2 = iris.cube.Cube(
    [1.0, -1.0, 3.0],
    standard_name=O_NAME,
)
L_CUBE = iris.cube.Cube(
    [10.0, 20.0, 0.0],
    standard_name=L_NAME,
)
FRAC_O = np.array([0.0, 1.0, 0.5, 0.3])
FRAC_L = np.array([0.1, 0.0, 1.0])

GET_LAND_FRACTION = [
    (CubeList([L_CUBE]), L_NAME, False, None),
    (CubeList([L_CUBE]), L_NAME, True, None),
    (CubeList([SFTLF_CUBE, L_CUBE]), L_NAME, False, FRAC_L),
    (CubeList([SFTLF_CUBE, O_CUBE_1]), O_NAME, False, None),
    (CubeList([SFTLF_CUBE, O_CUBE_1]), O_NAME, True, None),
    (CubeList([SFTLF_CUBE, O_CUBE_2]), O_NAME, False, FRAC_L),
    (CubeList([SFTLF_CUBE, O_CUBE_2]), O_NAME, True, FRAC_L),
    (CubeList([SFTOF_CUBE, L_CUBE]), L_NAME, False, None),
    (CubeList([SFTOF_CUBE, L_CUBE]), L_NAME, True, None),
    (CubeList([SFTOF_CUBE, O_CUBE_1]), O_NAME, False, None),
    (CubeList([SFTOF_CUBE, O_CUBE_1]), O_NAME, True, FRAC_O),
    (CubeList([SFTOF_CUBE, O_CUBE_2]), O_NAME, True, None),
    (CubeList([SFTOF_CUBE, SFTLF_CUBE, O_CUBE_2]), O_NAME, True, FRAC_L),
    (CubeList([SFTOF_CUBE, SFTLF_CUBE, O_CUBE_2]), O_NAME, False, FRAC_L),
]


@pytest.mark.parametrize('cubes,std_name,ocean_var,out', GET_LAND_FRACTION)
def test_get_land_fraction(cubes, std_name, ocean_var, out):
    """Test retrieving of land fraction from list of cubes."""
    land_fraction = shared._get_land_fraction(
        cubes, std_name, derive_from_ocean_fraction=ocean_var)
    if land_fraction is None or out is None:
        assert land_fraction is out
        return
    land_fraction = np.array(land_fraction)
    assert np.allclose(land_fraction, out)


SHAPES_TO_BROADCAST = [
    ((), (1, ), True),
    ((), (10, 10), True),
    ((1, ), (10, ), True),
    ((1, ), (10, 10), True),
    ((2, ), (10, ), False),
    ((10, ), (), True),
    ((10, ), (1, ), True),
    ((10, ), (10, ), True),
    ((10, ), (10, 10), True),
    ((10, ), (7, 1), True),
    ((10, ), (10, 7), False),
    ((10, ), (7, 1, 10), True),
    ((10, ), (7, 1, 1), True),
    ((10, ), (7, 1, 7), False),
    ((10, ), (7, 10, 7), False),
    ((10, 1), (1, 1), True),
    ((10, 1), (1, 100), True),
    ((10, 1), (10, 7), True),
    ((10, 12), (10, 1), True),
    ((10, 12), (), True),
    ((10, 12), (1, ), True),
    ((10, 12), (12, ), True),
    ((10, 12), (1, 1), True),
    ((10, 12), (1, 12), True),
    ((10, 12), (10, 10, 1), True),
    ((10, 12), (10, 12, 1), False),
    ((10, 12), (10, 12, 12), False),
    ((10, 12), (10, 10, 12), True),
]


@pytest.mark.parametrize('shape_1,shape_2,out', SHAPES_TO_BROADCAST)
def test_shape_is_broadcastable(shape_1, shape_2, out):
    """Test check if two shapes are broadcastable."""
    is_broadcastable = shared._shape_is_broadcastable(shape_1, shape_2)
    assert is_broadcastable == out


O_CUBE_1_OUT = O_CUBE_1.copy([1.0, 0.0, -0.5, 1.4])
O_CUBE_2_OUT = O_CUBE_2.copy([0.9, -1.0, 0.0])
O_CUBE_2_OUT_WRONG = O_CUBE_2.copy([0.1, 0.0, 3.0])
L_CUBE_OUT = L_CUBE.copy([1.0, 0.0, 0.0])
L_CUBE_OUT_WRONG = L_CUBE.copy([9.0, 20.0, 0.0])

CUBES_GRID_AREA_CORRECTION = [
    (CubeList([L_CUBE]), L_NAME, False, L_CUBE),
    (CubeList([L_CUBE]), L_NAME, True, L_CUBE),
    (CubeList([SFTLF_CUBE, L_CUBE]), L_NAME, False, L_CUBE_OUT),
    (CubeList([SFTLF_CUBE, L_CUBE]), L_NAME, True, L_CUBE_OUT_WRONG),
    (CubeList([SFTLF_CUBE, O_CUBE_1]), O_NAME, False, O_CUBE_1),
    (CubeList([SFTLF_CUBE, O_CUBE_1]), O_NAME, True, O_CUBE_1),
    (CubeList([SFTLF_CUBE, O_CUBE_2]), O_NAME, False, O_CUBE_2_OUT_WRONG),
    (CubeList([SFTLF_CUBE, O_CUBE_2]), O_NAME, True, O_CUBE_2_OUT),
    (CubeList([SFTOF_CUBE, O_CUBE_1]), O_NAME, False, O_CUBE_1),
    (CubeList([SFTOF_CUBE, O_CUBE_1]), O_NAME, True, O_CUBE_1_OUT),
    (CubeList([SFTOF_CUBE, O_CUBE_2]), O_NAME, False, O_CUBE_2),
    (CubeList([SFTOF_CUBE, O_CUBE_2]), O_NAME, True, O_CUBE_2),
    (CubeList([SFTOF_CUBE, SFTLF_CUBE, O_CUBE_1]), O_NAME, False, O_CUBE_1),
    (CubeList([SFTOF_CUBE, SFTLF_CUBE, O_CUBE_1]), O_NAME, True, O_CUBE_1_OUT),
    (CubeList([SFTOF_CUBE, SFTLF_CUBE, O_CUBE_2]), O_NAME, False,
     O_CUBE_2_OUT_WRONG),
    (CubeList([SFTOF_CUBE, SFTLF_CUBE, O_CUBE_2]), O_NAME, True, O_CUBE_2_OUT),
]


@pytest.mark.parametrize('cubes,std_name,ocean_var,out',
                         CUBES_GRID_AREA_CORRECTION)
def test_grid_area_correction(cubes, std_name, ocean_var, out):
    """Test grid area correction."""
    cubes = copy.deepcopy(cubes)
    cube = shared.grid_area_correction(cubes, std_name, ocean_var=ocean_var)
    assert cube == out
