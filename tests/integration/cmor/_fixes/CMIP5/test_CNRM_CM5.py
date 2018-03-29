import unittest

from cf_units import Unit
from iris.cube import Cube

from esmvaltool.cmor._fixes.CMIP5.CNRM_CM5 import msftmyz, msftmyzba


class TestMsftmyz(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='msftmyz', units='J')
        self.fix = msftmyz()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 1e6)
        self.assertEqual(cube.units, Unit('J'))


class TestMsftmyzba(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='msftmyzba', units='J')
        self.fix = msftmyzba()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 1e6)
        self.assertEqual(cube.units, Unit('J'))
