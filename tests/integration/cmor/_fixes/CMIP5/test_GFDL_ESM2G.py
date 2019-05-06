"""Tests for fixes of GFDL-ESM2G (CMIP5)."""
import unittest

import iris
import mock
import pytest
from cf_units import Unit

from esmvaltool.cmor._fixes.CMIP5.GFDL_ESM2G import (_get_and_remove, allvars,
                                                     co2, fgco2)

CUBE_1 = iris.cube.Cube([1.0], long_name='to_be_rm')
CUBE_2 = iris.cube.Cube([1.0], long_name='not_to_be_rm')
CUBES_LISTS = [
    (iris.cube.CubeList([CUBE_1]), iris.cube.CubeList([])),
    (iris.cube.CubeList([CUBE_1, CUBE_2]), iris.cube.CubeList([CUBE_2])),
    (iris.cube.CubeList([CUBE_2]), iris.cube.CubeList([CUBE_2])),
]


@pytest.mark.parametrize('cubes_in,cubes_out', CUBES_LISTS)
def test_get_and_remove(cubes_in, cubes_out):
    _get_and_remove(cubes_in, 'to_be_rm')
    assert cubes_in is not cubes_out
    assert cubes_in == cubes_out


CUBES = iris.cube.CubeList([CUBE_1, CUBE_2])


@mock.patch(
    'esmvaltool.cmor._fixes.CMIP5.GFDL_ESM2G._get_and_remove', autospec=True)
def test_allvars(mock_get_and_remove):
    fix = allvars()
    fix.fix_metadata(CUBES)
    assert mock_get_and_remove.call_count == 3
    assert mock_get_and_remove.call_args_list == [
        mock.call(CUBES, 'Start time for average period'),
        mock.call(CUBES, 'End time for average period'),
        mock.call(CUBES, 'Length of average period'),
    ]


@mock.patch(
    'esmvaltool.cmor._fixes.CMIP5.GFDL_ESM2G._get_and_remove', autospec=True)
def test_fgco2(mock_get_and_remove):
    fix = fgco2()
    fix.fix_metadata(CUBES)
    assert mock_get_and_remove.call_count == 2
    assert mock_get_and_remove.call_args_list == [
        mock.call(CUBES, 'Latitude of tracer (h) points'),
        mock.call(CUBES, 'Longitude of tracer (h) points'),
    ]


class TestCo2(unittest.TestCase):
    def setUp(self):
        self.cube = iris.cube.Cube([1.0], var_name='co2', units='J')
        self.fix = co2()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 1e6)
        self.assertEqual(cube.units, Unit('J'))
