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
                    'days since 0001-01-01 00:00:00', calendar='gregorian')),
            0)
        self.fix = allvars()

    def test_fix_data(self):
        cube = self.fix.fix_metadata(self.cube)

        time = cube.coord('time')
        self.assertEqual(time.units.origin, 'days since 1850-01-01 00:00:00')
        self.assertEqual(time.units.calendar, 'gregorian')


class TestCo2(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='co2', units='J')
        self.fix = co2()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 28.966 / 44.0)
        self.assertEqual(cube.units, Unit('J'))


class TestNbp(unittest.TestCase):
    def setUp(self):
        """Prepare temp folder for test"""
        self.fix = nbp()
        self.temp_folder = tempfile.mkdtemp()

    def tearDown(self):
        """Prepare temp folder for test"""
        shutil.rmtree(self.temp_folder)

    def test_fix_data(self):
        """Test fix on nbp files to set correct missing and fill values"""
        temp_handler, temp_path = tempfile.mkstemp('.nc', dir=self.temp_folder)
        os.close(temp_handler)
        output_dir = os.path.join(self.temp_folder, 'fixed')

        dataset = netCDF4.Dataset(temp_path, "w")
        var = dataset.createVariable('nbp', float, fill_value=1.0e20)
        var.missing_value = 1.0e20
        dataset.close()

        new_file = self.fix.fix_file(temp_path, output_dir)

        self.assertNotEqual(
            os.path.realpath(temp_path), os.path.realpath(new_file))

        dataset = netCDF4.Dataset(new_file)
        var = dataset.variables['nbp']
        self.assertEqual(var.missing_value, 1.0e33)
        self.assertEqual(var._FillValue, 1.0e33)
