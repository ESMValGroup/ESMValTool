import unittest

from cf_units import Unit
from iris.cube import Cube

from esmvaltool.cmor._fixes.CMIP5.EC_EARTH import sftlf, sic


class TestSic(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='sic', units='J')
        self.fix = sic()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 100)
        self.assertEqual(cube.units, Unit('J'))


class TestSftlf(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='sftlf', units='J')
        self.fix = sftlf()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 100)
        self.assertEqual(cube.units, Unit('J'))
