import unittest

from cf_units import Unit
from iris.cube import Cube

from esmvaltool.cmor._fixes.CMIP5.BNU_ESM import ch4, co2, fgco2, spco2


class TestCo2(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='co2', units='J')
        self.fix = co2()

    def test_fix_metadata(self):
        cube = self.fix.fix_metadata(self.cube)
        self.assertEqual(cube.units, Unit('1e-6'))
        self.assertEqual(cube.data[0], 1)

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 29.0 / 44.0 * 1.e6)
        self.assertEqual(cube.units, Unit('J'))


class Testfgco2(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='fgco2', units='J')
        self.fix = fgco2()

    def test_fix_metadata(self):
        cube = self.fix.fix_metadata(self.cube)
        self.assertEqual(cube.units, Unit('kg m-2 s-1'))
        self.assertEqual(cube.data[0], 1)

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 12.0 / 44.0)
        self.assertEqual(cube.units, Unit('J'))


class TestCh4(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='ch4', units='J')
        self.fix = ch4()

    def test_fix_metadata(self):
        cube = self.fix.fix_metadata(self.cube)
        self.assertEqual(cube.units, Unit('1e-9'))
        self.assertEqual(cube.data[0], 1)

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 29.0 / 16.0 * 1.e9)
        self.assertEqual(cube.units, Unit('J'))


class Testspco2(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='spco2', units='J')
        self.fix = spco2()

    def test_fix_metadata(self):
        cube = self.fix.fix_metadata(self.cube)
        self.assertEqual(cube.units, Unit('J'))
        self.assertEqual(cube.data[0], 1)

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 1.e6)
        self.assertEqual(cube.units, Unit('J'))
