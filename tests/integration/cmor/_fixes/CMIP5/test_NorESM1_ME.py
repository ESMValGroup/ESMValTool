"""Tests for fixes of NorESM1-ME (CMIP5)."""
import pytest
import iris
from iris.cube import CubeList

from esmvalcore.cmor._fixes.CMIP5.NorESM1_ME import tas

DIM_COORD_SHORT = iris.coords.DimCoord(
    [1.0, 2.0, 3.0],
    bounds=[[0.5, 1.5], [1.5, 2.5], [2.5, 3.5]],
    var_name='dim_coord',
)
DIM_COORD_LONG = iris.coords.DimCoord(
    [1.1234567891011, 2.1234567891011, 3.1234567891011],
    bounds=[
        [0.51234567891011, 1.51234567891011],
        [1.51234567891011, 2.51234567891011],
        [2.51234567891011, 3.51234567891011],
    ],
    var_name='dim_coord',
)
DIM_COORD_ROUNDED = iris.coords.DimCoord(
    [1.123456789101, 2.123456789101, 3.123456789101],
    bounds=[
        [0.512345678910, 1.512345678910],
        [1.512345678910, 2.512345678910],
        [2.512345678910, 3.512345678910],
    ],
    var_name='dim_coord',
)
AUX_COORD = iris.coords.AuxCoord(
    [1.1284712947128749498712, 2.12421841274128947982, 3.12787129852141124214],
    var_name='aux_coord',
)

CUBE_IN_SHORT = iris.cube.Cube(
    [3.14, 6.28, 9.42],
    dim_coords_and_dims=[(DIM_COORD_SHORT, 0)],
    aux_coords_and_dims=[(AUX_COORD, 0)],
)
CUBE_IN_LONG = iris.cube.Cube(
    [3.14, 6.28, 9.42],
    dim_coords_and_dims=[(DIM_COORD_LONG, 0)],
    aux_coords_and_dims=[(AUX_COORD, 0)],
)
CUBE_OUT_LONG = iris.cube.Cube(
    [3.14, 6.28, 9.42],
    dim_coords_and_dims=[(DIM_COORD_ROUNDED, 0)],
    aux_coords_and_dims=[(AUX_COORD, 0)],
)

CUBES_TO_FIX = [
    (CubeList([CUBE_IN_SHORT]), CubeList([CUBE_IN_SHORT])),
    (CubeList([CUBE_IN_LONG]), CubeList([CUBE_OUT_LONG])),
    (CubeList([CUBE_IN_LONG, CUBE_IN_SHORT]),
     CubeList([CUBE_OUT_LONG, CUBE_IN_SHORT])),
]


@pytest.mark.parametrize('cubes_in,cubes_out', CUBES_TO_FIX)
def test_tas(cubes_in, cubes_out):
    fix = tas()
    new_cubes = fix.fix_metadata(cubes_in)
    assert new_cubes is cubes_in
    assert new_cubes == cubes_out
