import unittest
from iris.cube import Cube
from cf_units import Unit
from orchestrator.interface_scripts.fixes.CMIP5.CESM1_BGC import co2


class TestCo2(unittest.TestCase):

    def setUp(self):
        self.cube = Cube([1], var_name='co2', units='J')
        self.fix = co2()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 28.966 / 44.0)
        self.assertEqual(cube.units, Unit('J'))

