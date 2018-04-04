import unittest

from cf_units import Unit
from iris.cube import Cube

from esmvaltool.cmor._fixes.CMIP5.inmcm4 import gpp, lai


class TestGpp(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='gpp', units='J')
        self.fix = gpp()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], -1)
        self.assertEqual(cube.units, Unit('J'))


class TestLai(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='lai', units='J')
        self.fix = lai()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 1.0 / 100.0)
        self.assertEqual(cube.units, Unit('J'))
