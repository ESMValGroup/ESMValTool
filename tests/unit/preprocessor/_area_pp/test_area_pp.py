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
    def setUp(self):
        shape = (3, 2, 2)
        data = np.arange(np.prod(shape)).reshape(shape)

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
        self.grid0d = iris.cube.Cube(data, dim_coords_and_dims=coords_spec)

        data2 = np.ones((2, 2))
        grid2 = iris.cube.Cube(data)
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
        self.grid2d = iris.cube.Cube(data2, dim_coords_and_dims=coords_spec2)

    def test_area_area_average(self):
        result = average_region(self.grid0d, 'latitude', 'longitude')
        expected = np.array([[[1.]]])
        self.assertArrayEqual(result.data, expected)

    def test_area_area_average_2d(self):
        result = average_region(self.grid2d, 'latitude', 'longitude')
        expected = np.array([1.])
        self.assertArrayEqual(result.data, expected)

    def test_area_extract_region(self):
        cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
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
        self.assertArrayEqual(result, self.grid2d)

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
        self.assertArrayEqual(result.data, self.grid2d.data)


if __name__ == '__main__':
    unittest.main()
