import unittest
from iris.cube import Cube
from cf_units import Unit
from orchestrator.interface_scripts.fixes.CMIP5.GFDL_ESM2M import sftof, co2


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
