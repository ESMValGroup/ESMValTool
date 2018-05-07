import unittest

from iris.cube import Cube

from esmvaltool.cmor._fixes.OBS.ESACCI_GHG import (
    xch4Stddev, xch4Stderr, xco2Stddev, xco2Stderr)


class Testxco2Stderr(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='xco2Stderr', units='1')
        self.fix = xco2Stderr()

    def test_fix_metadata(self):
        cube = self.fix.fix_metadata(self.cube)
        self.assertEqual(cube.units.origin, '1.0e-6')
        self.assertEqual(cube.data[0], 1)

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 1e6)


class Testxco2Stddev(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='xco2Stddev', units='1')
        self.fix = xco2Stddev()

    def test_fix_metadata(self):
        cube = self.fix.fix_metadata(self.cube)
        self.assertEqual(cube.units.origin, '1.0e-6')
        self.assertEqual(cube.data[0], 1)

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 1e6)


class Testxch4Stddev(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='xch4Stddev', units='1')
        self.fix = xch4Stddev()

    def test_fix_metadata(self):
        cube = self.fix.fix_metadata(self.cube)
        self.assertEqual(cube.units.origin, '1.0e-9')
        self.assertEqual(cube.data[0], 1)

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 1e9)


class Testxch4Stderr(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='xch4Stderr', units='1')
        self.fix = xch4Stderr()

    def test_fix_metadata(self):
        cube = self.fix.fix_metadata(self.cube)
        self.assertEqual(cube.units.origin, '1.0e-9')
        self.assertEqual(cube.data[0], 1)

    def test_fix_data(self):
        cube = self.fix.fix_data(self.cube)
        self.assertEqual(cube.data[0], 1e9)
