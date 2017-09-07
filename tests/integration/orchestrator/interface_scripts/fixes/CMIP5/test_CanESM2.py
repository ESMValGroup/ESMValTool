import unittest
from iris.cube import Cube
from cf_units import Unit
from orchestrator.interface_scripts.fixes.CMIP5.CanESM2 import fgco2


class TestCanESM2Fgco2(unittest.TestCase):

    def setUp(self):
        self.cube = Cube([1], var_name='fgco2', units='J')
        self.fix = fgco2()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 12.0 / 44.0)
        self.assertEqual(cube.units, Unit('J'))

