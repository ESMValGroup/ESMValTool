"""Unit test for the :func:`esmvaltool.preprocessor._mask` function"""

import unittest

import iris
import numpy as np

import tests
from esmvaltool.preprocessor._mask import (
    mask_above_threshold, mask_below_threshold, mask_inside_range,
    mask_outside_range)


class Test(tests.Test):
    """Test class for _mask"""

    def setUp(self):
        """Prepare tests"""
        coord_sys = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        self.data2 = np.array([[0., 1.], [2., 3.]])
        lons2 = iris.coords.DimCoord([1.5, 2.5],
                                     standard_name='longitude',
                                     bounds=[[1., 2.], [2., 3.]],
                                     units='degrees_east',
                                     coord_system=coord_sys)
        lats2 = iris.coords.DimCoord([1.5, 2.5],
                                     standard_name='latitude',
                                     bounds=[[1., 2.], [2., 3.]],
                                     units='degrees_north',
                                     coord_system=coord_sys)
        coords_spec3 = [(lats2, 0), (lons2, 1)]
        self.arr = iris.cube.Cube(self.data2, dim_coords_and_dims=coords_spec3)

    def test_mask_above_threshold(self):
        """Test to mask above a threshold."""
        result = mask_above_threshold(self.arr, 1.5)
        expected = np.ma.array(self.data2, mask=[[False, False], [True, True]])
        self.assertArrayEqual(result.data, expected)

    def test_mask_below_threshold(self):
        """Test to mask below a threshold."""
        result = mask_below_threshold(self.arr, 1.5)
        expected = np.ma.array(self.data2, mask=[[True, True], [False, False]])
        self.assertArrayEqual(result.data, expected)

    def test_mask_inside_range(self):
        """Test to mask inside a range."""
        result = mask_inside_range(self.arr, 0.5, 2.5)
        expected = np.ma.array(self.data2, mask=[[False, True], [True, False]])
        self.assertArrayEqual(result.data, expected)

    def test_mask_outside_range(self):
        """Test to mask outside a range."""
        result = mask_outside_range(self.arr, 0.5, 2.5)
        expected = np.ma.array(self.data2, mask=[[True, False], [False, True]])
        self.assertArrayEqual(result.data, expected)


if __name__ == '__main__':
    unittest.main()
