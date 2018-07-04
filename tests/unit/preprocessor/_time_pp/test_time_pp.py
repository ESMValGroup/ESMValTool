"""
Unit tests for the :func:`esmvaltool.preprocessor.regrid.regrid` function.

"""

from __future__ import absolute_import, division, print_function

import unittest

import iris
import mock
import numpy as np
from cf_units import Unit

import tests

from esmvaltool.preprocessor._time_pp import time_average


class Test(tests.Test):

    def test_time_time_average(self):
        data = np.ones((3))
        cube = iris.cube.Cube(data)

        times = np.array([15., 45., 75.])
        bounds = np.array([[0., 30.], [30., 60.], [60., 90.]])
        time = iris.coords.DimCoord(times,
                                    bounds=bounds,
                                    standard_name='time',
                                    units=Unit('days since 1950-01-01',
                                               calendar='gregorian'))
        cube = iris.cube.Cube(data, dim_coords_and_dims=[(time, 0)])

        result = time_average(cube)
        expected = np.array([1.])
        self.assertArrayEqual(result.data, expected)

    def test_time_time_average_uneven(self):
        data = np.array([1., 5.])
        cube = iris.cube.Cube(data)

        times = np.array([5., 25.])
        bounds = np.array([[0., 1.], [1., 4.]])
        time = iris.coords.DimCoord(times,
                                    bounds=bounds,
                                    standard_name='time',
                                    units=Unit('days since 1950-01-01',
                                               calendar='gregorian'))
        cube = iris.cube.Cube(data, dim_coords_and_dims=[(time, 0)])

        result = time_average(cube)
        expected = np.array([4.])
        self.assertArrayEqual(result.data, expected)

    def test_time_time_average_365_day(self):
        data = np.ones((6,))
        cube = iris.cube.Cube(data)

        times = np.array([15, 45, 74, 105, 135, 166])
        bounds = np.array([[0, 31], [31, 59], [59, 90], [90, 120],
                           [120, 151], [151, 181]])
        time = iris.coords.DimCoord(times,
                                    bounds=bounds,
                                    standard_name='time',
                                    var_name='time',
                                    units=Unit('days since 1950-01-01',
                                               calendar='365_day'))
        cube = iris.cube.Cube(data, dim_coords_and_dims=[(time, 0)])

        result = time_average(cube)
        expected = np.array([1.])
        self.assertArrayEqual(result.data, expected)


if __name__ == '__main__':
    unittest.main()
