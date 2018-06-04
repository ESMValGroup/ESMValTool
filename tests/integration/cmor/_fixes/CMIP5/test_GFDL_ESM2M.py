import unittest

from cf_units import Unit
from iris.coords import DimCoord
from iris.cube import Cube

from esmvaltool.cmor._fixes.CMIP5.GFDL_ESM2M import allvars, co2, sftof


class TestAll(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1, 2], var_name='co2', units='J')
        self.cube.add_dim_coord(
            DimCoord(
                [0, 1],
                standard_name='time',
                units=Unit(
                    'days since 0001-01-01 00:00:00', calendar='gregorian')),
            0)
        self.fix = allvars()

    def test_fix_data(self):
        cube = self.fix.fix_metadata(self.cube)

        time = cube.coord('time')
        self.assertEqual(time.units.origin, 'days since 1850-01-01 00:00:00')
        self.assertEqual(time.units.calendar, 'gregorian')


class TestSftof(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='sftof', units='J')
        self.fix = sftof()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 100)
        self.assertEqual(cube.units, Unit('J'))


class TestCo2(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='co2', units='J')
        self.fix = co2()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 1e6)
        self.assertEqual(cube.units, Unit('J'))
