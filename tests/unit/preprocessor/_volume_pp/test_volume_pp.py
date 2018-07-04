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
    def test_area_volume_slice(self):
        """
        Test to extract the top two layers of a 3 layer depth column.
        """    
        cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        data = np.arange(3)
        grid = iris.cube.Cube(data)

        depth = iris.coords.DimCoord(
            [0.5, 5., 50.],
            standard_name='depth',
            bounds=[[0., 2.5], [2.5, 25.], [25., 250.]],
            units='m',
            coord_system=cs)

        coords_spec = [(depth, 0), ]
        grid = iris.cube.Cube(data, dim_coords_and_dims=coords_spec)

        result = volume_slice(grid, 0., 10.)

        expected = np.array([0., 1.])
        self.assertArrayEqual(result.data, expected)

    def test_area_volume_average(self):
        """
        Test to take the volume weighted average of a (3,2,2) cube.
        """
        cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        data2 = np.ones((3, 2, 2))

        grid2 = iris.cube.Cube(data2)
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
            units='m',
            coord_system=cs)
        coords_spec3 = [(depth, 0), (lats2, 1), (lons2, 2)]
        grid2 = iris.cube.Cube(data2, dim_coords_and_dims=coords_spec3)

        result = volume_average(grid2, 'depth', 'latitude', 'longitude')

        expected = np.array([1.])
        self.assertArrayEqual(result.data, expected)

    def test_area_depth_integration(self):
        """
        Test to take the depth integration of a 3 layer cube.
        """    
        cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        data = np.ones((3, ))
        grid = iris.cube.Cube(data)

        depth = iris.coords.DimCoord(
            [0.5, 5., 50.],
            standard_name='depth',
            bounds=[[0., 2.5], [2.5, 25.], [25., 250.]],
            units='m',
            coord_system=cs)

        coords_spec = [(depth, 0), ]
        grid = iris.cube.Cube(data, dim_coords_and_dims=coords_spec)

        result = depth_integration(grid, 'depth')

        expected = np.array([250., ])
        self.assertArrayEqual(result.data, expected)
        
    def test_area_extract_transect(self):
        """
        Test to extract a transect from a (3,2,2) cube.
        """
        cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        data2 = np.ones((3, 2, 2))

        grid2 = iris.cube.Cube(data2)
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
            units='m',
            coord_system=cs)
        coords_spec3 = [(depth, 0), (lats2, 1), (lons2, 2)]
        grid2 = iris.cube.Cube(data2, dim_coords_and_dims=coords_spec3)

        result = extract_transect(grid2, latitude= 1.5)

        expected = np.ones((3, 2))
        self.assertArrayEqual(result.data, expected)
        
    def test_area_extract_trajectory(self):
        """
        Test to extract a trajectory from a (3,2,2) cube.
        """
        cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        data2 = np.ones((3, 2, 2))

        grid2 = iris.cube.Cube(data2)
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
            units='m',
            coord_system=cs)
        coords_spec3 = [(depth, 0), (lats2, 1), (lons2, 2)]
        grid2 = iris.cube.Cube(data2, dim_coords_and_dims=coords_spec3)

        result = extract_trajectory(grid2, [1.5, 2.5], [2., 2.], 2)

        expected = np.ones((3, 2))
        self.assertArrayEqual(result.data, expected)
        
        

if __name__ == '__main__':
    unittest.main()
