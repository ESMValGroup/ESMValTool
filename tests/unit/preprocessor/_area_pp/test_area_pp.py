"""
Unit tests for the :func:`esmvaltool.preprocessor.regrid.regrid` function.

"""

from __future__ import absolute_import, division, print_function

import unittest

import iris
import mock
import numpy as np

import tests

from esmvaltool.preprocessor._area_pp import area_slice as extract_region
from esmvaltool.preprocessor._area_pp import area_average as average_region


class Test(tests.Test):

    def test_area_area_average(self):
        cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        data = np.ones((1, 1))
        grid = iris.cube.Cube(data)
        lons = iris.coords.DimCoord(
            [1.5],
            standard_name='longitude',
            bounds=[[1., 2.]],
            units='degrees_east',
            coord_system=cs)
        lats = iris.coords.DimCoord(
            [1.5],
            standard_name='latitude',
            bounds=[[1., 2.]],
            units='degrees_north',
            coord_system=cs)
        coords_spec = [(lats, 0), (lons, 1)]
        grid0d = iris.cube.Cube(data, dim_coords_and_dims=coords_spec)

        result = average_region(grid0d, 'latitude', 'longitude')
        expected = np.array([[[1.]]])
        self.assertArrayEqual(result.data, expected)

    def test_area_area_average_2d(self):
        cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        data2 = np.ones((2, 2))
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
        coords_spec2 = [(lats2, 0), (lons2, 1)]
        grid = iris.cube.Cube(data2, dim_coords_and_dims=coords_spec2)

        result = average_region(grid, 'latitude', 'longitude')
        expected = np.array([1.])
        self.assertArrayEqual(result.data, expected)

    def test_area_extract_region(self):
        cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        # data to sample
        data = np.ones((5, 5))
        grid = iris.cube.Cube(data)
        lons = iris.coords.DimCoord(
            [i + .5 for i in range(5)],
            standard_name='longitude',
            bounds=[[i, i + 1.] for i in range(5)],  # [0,1] to [4,5]
            units='degrees_east',
            coord_system=cs)
        lats = iris.coords.DimCoord(
            [i + .5 for i in range(5)],
            standard_name='latitude',
            bounds=[[i, i + 1.] for i in range(5)],
            units='degrees_north',
            coord_system=cs)
        coords_spec = [(lats, 0), (lons, 1)]
        grid = iris.cube.Cube(data, dim_coords_and_dims=coords_spec)
        result = extract_region(grid, 1.5, 2.5, 1.5, 2.5)

        # expected outcome
        data2 = np.ones((2, 2))
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
        coords_spec2 = [(lats2, 0), (lons2, 1)]
        expected = iris.cube.Cube(data2, dim_coords_and_dims=coords_spec2)

        self.assertArrayEqual(result, expected)

    def test_area_extract_region2(self):
        cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        data = np.ones((6, 6))
        grid = iris.cube.Cube(data)

        lons = iris.coords.DimCoord(
            [i - 2.5 for i in range(6)],
            standard_name='longitude',
            bounds=[[i - 3., i - 2.] for i in range(6)],  # [3,2] to [4,5]
            units='degrees_east',
            coord_system=cs)
        lats = iris.coords.DimCoord(
            [i - 2.5 for i in range(6)],
            standard_name='latitude',
            bounds=[[i - 3., i - 2.] for i in range(6)],
            units='degrees_north',
            coord_system=cs)
        coords_spec = [(lats, 0), (lons, 1)]
        grid = iris.cube.Cube(data, dim_coords_and_dims=coords_spec)
        result = extract_region(grid, -0.5, 0.5, -0.5, 0.5)

        # expected outcome
        data2 = np.ones((2, 2))
        grid2 = iris.cube.Cube(data2)
        lons2 = iris.coords.DimCoord(
            [-0.5, 0.5],
            standard_name='longitude',
            bounds=[[-1., 0.], [0., 1.]],
            units='degrees_east',
            coord_system=cs)
        lats2 = iris.coords.DimCoord(
            [-0.5, 0.5],
            standard_name='latitude',
            bounds=[[-1., 0.], [0., 1.]],
            units='degrees_north',
            coord_system=cs)
        coords_spec2 = [(lats2, 0), (lons2, 1)]
        expected = iris.cube.Cube(data2, dim_coords_and_dims=coords_spec2)

        self.assertArrayEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
