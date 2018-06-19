import unittest

from cf_units import Unit
from iris.cube import Cube

from esmvaltool.cmor._fixes.CMIP5.MPI_ESM_LR import pctisccp


class TestPctisccp2(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='pctisccp', units='J')
        self.fix = pctisccp()

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 100)
        self.assertEqual(cube.units, Unit('J'))
