import unittest

from cf_units import Unit
from iris.cube import Cube

from esmvaltool.cmor._fixes.CMIP5.GFDL_CM2p1 import sftof


class TestSftof(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='sftof', units='J')
        self.fix = sftof()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 100)
        self.assertEqual(cube.units, Unit('J'))
