import unittest

from cf_units import Unit
from iris.cube import Cube

from esmvaltool.cmor._fixes.CMIP5.GFDL_ESM2G import co2


class TestCo2(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='co2', units='J')
        self.fix = co2()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 1e6)
        self.assertEqual(cube.units, Unit('J'))
