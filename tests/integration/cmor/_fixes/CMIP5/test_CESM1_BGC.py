import os
import shutil
import tempfile
import unittest

import netCDF4
from cf_units import Unit
from iris.coords import DimCoord
from iris.cube import Cube

from esmvaltool.cmor._fixes.CMIP5.CESM1_BGC import allvars, co2, nbp


class TestAll(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1, 2], var_name='co2', units='J')
        self.cube.add_dim_coord(
            DimCoord(
                [0, 1],
                standard_name='time',
                units=Unit(
                    'days since 0001-01-01 00:00:00.0000000 UTC',
                    calendar='gregorian')),
            0)
        self.fix = allvars()

    def test_fix_metadata(self):
        cube = self.fix.fix_metadata(self.cube)

        time = cube.coord('time')
        self.assertEqual(time.units.origin, 'days since 1850-01-01 00:00:00')
        self.assertEqual(time.units.calendar, 'gregorian')

    def test_fix_metadata_good_units(self):
        self.cube.coord('time').units = Unit('days since 1950-01-01 00:00:00',
                                             calendar='gregorian')
        cube = self.fix.fix_metadata(self.cube)

        time = cube.coord('time')
        self.assertEqual(time.units.origin, 'days since 1950-01-01 00:00:00')
        self.assertEqual(time.units.calendar, 'gregorian')


class TestCo2(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='co2', units='J')
        self.fix = co2()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 28.966 / 44.0)
        self.assertEqual(cube.units, Unit('J'))
