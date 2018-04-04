import unittest

from iris.coords import DimCoord
from iris.cube import Cube

from esmvaltool.cmor._fixes.OBS.BDBP import tro3prof


class TestTro3prof(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1, 2], var_name='tro3prof', units='J')
        self.cube.add_dim_coord(
            DimCoord([1, 2], standard_name='air_pressure', units='hPa'), 0)
        self.fix = tro3prof()

    def test_fix_metadata(self):
        cube = self.fix.fix_metadata(self.cube)
        self.assertEqual(cube.coord('air_pressure').units.origin, 'Pa')
        self.assertEqual(cube.coord('air_pressure').points[0], 100)
        self.assertEqual(cube.coord('air_pressure').points[1], 200)
        self.assertIsNone(cube.coord('air_pressure').bounds)

    def test_fix_metadata_with_bounds(self):
        self.cube = Cube([1, 2], var_name='tro3prof', units='J')
        self.cube.add_dim_coord(
            DimCoord(
                [1, 2],
                standard_name='air_pressure',
                units='hPa',
                bounds=[[0.5, 1.5], [1.5, 2.5]]), 0)
        cube = self.fix.fix_metadata(self.cube)

        plev = cube.coord('air_pressure')
        self.assertEqual(plev.units.origin, 'Pa')
        self.assertEqual(plev.points[0], 100)
        self.assertEqual(plev.points[1], 200)
        self.assertEqual(plev.bounds[0][0], 50)
        self.assertEqual(plev.bounds[0][1], 150)
        self.assertEqual(plev.bounds[1][0], 150)
        self.assertEqual(plev.bounds[1][1], 250)
