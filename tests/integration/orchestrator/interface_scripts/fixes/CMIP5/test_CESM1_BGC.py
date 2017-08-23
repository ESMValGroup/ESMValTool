import unittest
from iris.cube import Cube
from cf_units import Unit
from orchestrator.interface_scripts.fixes.CMIP5.CESM1_BGC import co2, nbp
import tempfile
import netCDF4


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
        self.fix = nbp()

    def test_fix_data(self):
        temp_handler, temp_path = tempfile.mkstemp('.nc')
        dataset = netCDF4.Dataset(temp_path, "w")
        var = dataset.createVariable('nbp', float, fill_value=1.0e20)
        var.missing_value = 1.0e20
        dataset.close()

        new_file = self.fix.fix_file(temp_path)

        dataset = netCDF4.Dataset(new_file)
        var = dataset.variables['nbp']
        self.assertEqual(new_file, temp_path)
        self.assertEqual(var.missing_value, 1.0e33)
        self.assertEqual(var._FillValue, 1.0e33)
