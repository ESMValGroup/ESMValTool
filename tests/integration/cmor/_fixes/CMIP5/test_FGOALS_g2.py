import unittest

from cf_units import Unit
from iris.coords import DimCoord
from iris.cube import Cube

from esmvalcore.cmor._fixes.CMIP5.FGOALS_g2 import allvars


class TestAll(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1.0, 2.0], var_name='co2', units='J')
        self.cube.add_dim_coord(
            DimCoord(
                [0.0, 1.0],
                standard_name='time',
                units=Unit('days since 0001-01', calendar='gregorian')),
            0)
        self.fix = allvars()

    def test_fix_metadata(self):
        cube = self.fix.fix_metadata([self.cube])[0]

        time = cube.coord('time')
        self.assertEqual(time.units.origin,
                         'day since 1-01-01 00:00:00.000000')
        self.assertEqual(time.units.calendar, 'gregorian')
