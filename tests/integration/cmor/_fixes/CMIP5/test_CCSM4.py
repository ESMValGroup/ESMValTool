import unittest

import numpy as np
from iris.coords import DimCoord
from iris.cube import Cube

from esmvaltool.cmor._fixes.CMIP5.CCSM4 import rlut, rlutcs, so


class TestsRlut(unittest.TestCase):
    """Test for rlut fixes"""

    def setUp(self):
        """Prepare tests"""
        self.cube = Cube([1.0, 2.0], var_name='rlut')
        self.cube.add_dim_coord(
            DimCoord([0.50001, 1.499999],
                     standard_name='latitude',
                     bounds=[
                         [0.00001, 0.999999],
                         [1.00001, 1.999999],
                     ]), 0)
        self.fix = rlut()

    def test_fix_metadata(self):
        """Check that latitudes values are rounded"""
        cube = self.fix.fix_metadata([self.cube])[0]

        latitude = cube.coord('latitude')
        self.assertTrue(np.all(latitude.points == np.array([0.5000, 1.5000])))
        self.assertTrue(
            np.all(latitude.bounds == np.array([[0.0000, 1.0000],
                                                [1.0000, 2.0000]])))


class TestsRlutcs(unittest.TestCase):
    """Test for rlutcs fixes"""

    def setUp(self):
        """Prepare tests"""
        self.cube = Cube([1.0, 2.0], var_name='rlutcs')
        self.cube.add_dim_coord(
            DimCoord([0.50001, 1.499999],
                     standard_name='latitude',
                     bounds=[
                         [0.00001, 0.999999],
                         [1.00001, 1.999999],
                     ]), 0)
        self.fix = rlutcs()

    def test_fix_metadata(self):
        """Check that latitudes values are rounded"""
        cube = self.fix.fix_metadata([self.cube])[0]

        latitude = cube.coord('latitude')
        self.assertTrue(np.all(latitude.points == np.array([0.5000, 1.5000])))
        self.assertTrue(
            np.all(latitude.bounds == np.array([[0.0000, 1.0000],
                                                [1.0000, 2.0000]])))


class TestSo(unittest.TestCase):
    """Tests for so fixes"""

    def setUp(self):
        """Prepare tests"""
        self.cube = Cube([1.0, 2.0], var_name='so', units='1.0')
        self.fix = so()

    def test_fix_metadata(self):
        """Checks that units are changed to the correct value"""
        cube = self.fix.fix_metadata([self.cube])[0]
        self.assertEqual('1e3', cube.units.origin)
