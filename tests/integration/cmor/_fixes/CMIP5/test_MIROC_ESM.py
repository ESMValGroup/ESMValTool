import unittest

from cf_units import Unit
from iris.coords import DimCoord
from iris.cube import Cube
from iris.exceptions import CoordinateNotFoundError

from esmvaltool.cmor._fixes.CMIP5.MIROC_ESM import allvars, co2, gpp, tro3


class TestCo2(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='co2', units='J')
        self.fix = co2()

    def test_fix_metadata(self):
        cube = self.fix.fix_metadata(self.cube)
        self.assertEqual(cube.data[0], 1)
        self.assertEqual(cube.units, Unit('1e-6'))


class TestTro3(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='tro3', units='J')
        self.fix = tro3()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 1000)
        self.assertEqual(cube.units, Unit('J'))


class TestGpp(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='gpp', units='J')
        self.fix = gpp()

    def test_fix_metadata(self):
        cube = self.fix.fix_metadata(self.cube)
        self.assertEqual(cube.data[0], 1)
        self.assertEqual(cube.units, Unit('g m-2 day-1'))


class TestAll(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([[1, 2], [3, 4]], var_name='co2', units='J')
        self.cube.add_dim_coord(
            DimCoord(
                [0, 1],
                standard_name='time',
                units=Unit(
                    'days since 0000-01-01 00:00:00', calendar='gregorian')),
            0)
        self.cube.add_dim_coord(DimCoord([0, 1], long_name='AR5PL35'), 1)

        self.fix = allvars()

    def test_fix_metadata(self):
        cube = self.fix.fix_metadata(self.cube)
        time = cube.coord('time')
        self.assertEqual(time.units.origin, 'days since 1849-01-01 00:00:00')
        self.assertEqual(time.units.calendar, 'gregorian')

    def test_fix_metadata_1_1(self):
        time = self.cube.coord('time')
        time.units = Unit("days since 1-1-1", time.units.calendar)
        cube = self.fix.fix_metadata(self.cube)

        time = cube.coord('time')
        self.assertEqual(time.units.origin, 'days since 1850-01-01 00:00:00')
        self.assertEqual(time.units.calendar, 'gregorian')

    def test_fix_metadata_plev(self):
        time = self.cube.coord('time')
        time.units = Unit("days since 1-1-1", time.units.calendar)
        cube = self.fix.fix_metadata(self.cube)
        cube.coord('air_pressure')

    def test_fix_metadata_no_plev(self):
        self.cube.remove_coord('AR5PL35')
        cube = self.fix.fix_metadata(self.cube)
        with self.assertRaises(CoordinateNotFoundError):
            cube.coord('air_pressure')


# if (iscoord(var, "time")) then
#     if (isatt(var&time,"units"))then
#         if (var&time@units.eq."days since 0000-01-01 00:00:00") then
#             var&time@units ="days since 1849-01-01 00:00:00"
#             ret = 0
#         end if
#         if (var&time@units.eq."days since 1-1-1")then
#             var&time@units ="days since 1850-01-01 00:00:00"
#             ret = 0
#         end if
#     end if
# end if
#
# if (iscoord(var, "AR5PL35")) then
#     var!1 = "plev"
#     ret = 0
# end if
