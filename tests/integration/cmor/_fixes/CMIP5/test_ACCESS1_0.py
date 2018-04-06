import unittest

from iris.cube import Cube

from esmvaltool.cmor._fixes.CMIP5.ACCESS1_0 import rlut


class TestRlut(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='rlut')
        self.fix = rlut()

    def test_fix_metadata(self):
        cube = self.fix.fix_metadata(self.cube)
        self.assertEqual(cube.attributes['positive'], 'up')

class TestRlutcs(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([1], var_name='rlutcs')
        self.fix = rlut()

    def test_fix_metadata(self):
        cube = self.fix.fix_metadata(self.cube)
        self.assertEqual(cube.attributes['positive'], 'up')
