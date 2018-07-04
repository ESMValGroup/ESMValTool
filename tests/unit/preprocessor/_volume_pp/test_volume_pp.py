"""
Unit tests for the :func:`esmvaltool.preprocessor.regrid.regrid` function.

"""

from __future__ import absolute_import, division, print_function

import unittest

import iris
import numpy as np

import tests
from esmvaltool.preprocessor._volume_pp import volume_slice
from esmvaltool.preprocessor._volume_pp import volume_average
from esmvaltool.preprocessor._volume_pp import depth_integration
from esmvaltool.preprocessor._volume_pp import extract_transect
from esmvaltool.preprocessor._volume_pp import extract_trajectory


class Test(tests.Test):
    def setUp(self):
        """ Prepare tests """
        cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        data2 = np.ones((3, 2, 2))

        grid_3d = iris.cube.Cube(data2)
        lons2 = iris.coords.DimCoord(
            [1.5, 2.5],
            standard_name='longitude',
            bounds=[[1., 2.], [2., 3.]],
            units='degrees_east',
            coord_system=cs)
        lats2 = iris.coords.DimCoord(
            [1.5, 2.5],
            standard_name='latitude',
            bounds=[[1., 2.], [2., 3.]],
            units='degrees_north',
            coord_system=cs)
        depth = iris.coords.DimCoord(
            [0.5, 5., 50.],
            standard_name='depth',
            bounds=[[0., 2.5], [2.5, 25.], [25., 250.]],
            units='m')
        coords_spec3 = [(depth, 0), (lats2, 1), (lons2, 2)]
        self.grid_3d = iris.cube.Cube(data2, dim_coords_and_dims=coords_spec3)

    def test_area_volume_slice(self):
        """
        Test to extract the top two layers of a 3 layer depth column.
        """
        result = volume_slice(self.grid_3d, 0., 10.)
        expected = np.ones((2, 2, 2))
        print(result.data, expected.data)
        self.assertArrayEqual(result.data, expected)

    def test_area_volume_average(self):
        """
        Test to take the volume weighted average of a (3,2,2) cube.
        """
        cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        result = volume_average(self.grid_3d, 'depth', 'latitude', 'longitude')
        expected = np.array([1.])
        self.assertArrayEqual(result.data, expected)

    def test_area_depth_integration_1D(self):
        """
        Test to take the depth integration of a 3 layer cube.
        """
        result = depth_integration(self.grid_3d[:, 0, 0], 'depth')
        expected = np.ones((1, 1)) * 250.
        print(result.data, expected.data)
        self.assertArrayEqual(result.data, expected)

    def test_area_depth_integration_3D(self):
        """
        Test to take the depth integration of a 3 layer cube.
        """
        result = depth_integration(self.grid_3d, 'depth')
        expected = np.ones((2, 2)) * 250.
        print(result.data, expected.data)
        self.assertArrayEqual(result.data, expected)

    def test_area_extract_transect_latitude(self):
        """
        Test to extract a transect from a (3, 2, 2) cube.
        """
        result = extract_transect(self.grid_3d, latitude=1.5)
        expected = np.ones((3, 2))
        self.assertArrayEqual(result.data, expected)

    def test_area_extract_transect_longitude(self):
        """
        Test to extract a transect from a (3, 2, 2) cube.
        """
        result = extract_transect(self.grid_3d, longitude=1.5)
        expected = np.ones((3, 2))
        self.assertArrayEqual(result.data, expected)

    def test_area_extract_trajectory(self):
        """
        Test to extract a trajectory from a (3, 2, 2) cube.
        """
        result = extract_trajectory(self.grid_3d, [1.5, 2.5], [2., 2.], 2)
        expected = np.ones((3, 2))
        self.assertArrayEqual(result.data, expected)


if __name__ == '__main__':
    unittest.main()
