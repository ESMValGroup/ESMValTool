"""Unit tests for the :func:`esmvaltool.preprocessor._time_pp` module"""

from __future__ import absolute_import, division, print_function

import unittest

import iris
import iris.coord_categorisation
import iris.coords
import numpy as np
from cf_units import Unit
from iris.cube import Cube

import tests
from esmvaltool.preprocessor._time_area import (extract_month, extract_season,
                                                time_average)


def _create_sample_cube():
    cube = Cube(np.arange(1, 25), var_name='co2', units='J')
    cube.add_dim_coord(
        iris.coords.DimCoord(
            np.arange(15., 720., 30.),
            standard_name='time',
            units=Unit('days since 1950-01-01', calendar='gregorian')), 0)
    iris.coord_categorisation.add_month_number(cube, 'time')
    return cube


class TestExtractMonth(tests.Test):
    """Tests for extract_month`."""

    def setUp(self):
        """Prepare tests"""
        self.cube = _create_sample_cube()

    def test_get_january(self):
        """Test january extraction"""
        sliced = extract_month(self.cube, 1)
        print(sliced)
        self.assertTrue(
            (np.array([1, 1]) == sliced.coord('month_number').points).all())


class TestExtractSeason(tests.Test):
    """Tests for extract_season"""

    def setUp(self):
        """Prepare tests"""
        self.cube = _create_sample_cube()

    def test_get_djf(self):
        """Test function for winter"""
        sliced = extract_season(self.cube, 'djf')
        print(sliced)
        self.assertTrue(
            (np.array([1, 2, 12, 1, 2,
                       12]) == sliced.coord('month_number').points).all())

    def test_get_djf_caps(self):
        """Test function works when season specified in caps"""
        sliced = extract_season(self.cube, 'DJF')
        print(sliced)
        self.assertTrue(
            (np.array([1, 2, 12, 1, 2,
                       12]) == sliced.coord('month_number').points).all())

    def test_get_mam(self):
        """Test function for spring"""
        sliced = extract_season(self.cube, 'mam')
        print(sliced)
        self.assertTrue((np.array(
            [3, 4, 5, 3, 4, 5]) == sliced.coord('month_number').points).all())

    def test_get_jja(self):
        """Test function for summer"""
        sliced = extract_season(self.cube, 'jja')
        print(sliced)
        self.assertTrue((np.array(
            [6, 7, 8, 6, 7, 8]) == sliced.coord('month_number').points).all())

    def test_get_son(self):
        """Test function for summer"""
        sliced = extract_season(self.cube, 'son')
        print(sliced)
        self.assertTrue(
            (np.array([9, 10, 11, 9, 10,
                       11]) == sliced.coord('month_number').points).all())


class TestTimeAverage(tests.Test):
    """Test class for the :func:`esmvaltool.preprocessor._time_pp` module"""

    def test_time_average(self):
        """Test for time average of a 1D field."""
        data = np.ones((3))
        cube = iris.cube.Cube(data)

        times = np.array([15., 45., 75.])
        bounds = np.array([[0., 30.], [30., 60.], [60., 90.]])
        time = iris.coords.DimCoord(
            times,
            bounds=bounds,
            standard_name='time',
            units=Unit('days since 1950-01-01', calendar='gregorian'))
        cube = iris.cube.Cube(data, dim_coords_and_dims=[(time, 0)])

        result = time_average(cube)
        expected = np.array([1.])
        self.assertArrayEqual(result.data, expected)

    def test_time_average_uneven(self):
        """Test for time average of a 1D field with uneven time boundaries."""
        data = np.array([1., 5.])
        cube = iris.cube.Cube(data)

        times = np.array([5., 25.])
        bounds = np.array([[0., 1.], [1., 4.]])
        time = iris.coords.DimCoord(
            times,
            bounds=bounds,
            standard_name='time',
            units=Unit('days since 1950-01-01', calendar='gregorian'))
        cube = iris.cube.Cube(data, dim_coords_and_dims=[(time, 0)])

        result = time_average(cube)
        expected = np.array([4.])
        self.assertArrayEqual(result.data, expected)

    def test_time_average_365_day(self):
        """Test for time avg of a realisitc time axis and 365 day calendar"""
        data = np.ones((6, ))
        cube = iris.cube.Cube(data)

        times = np.array([15, 45, 74, 105, 135, 166])
        bounds = np.array([[0, 31], [31, 59], [59, 90], [90, 120], [120, 151],
                           [151, 181]])
        time = iris.coords.DimCoord(
            times,
            bounds=bounds,
            standard_name='time',
            var_name='time',
            units=Unit('days since 1950-01-01', calendar='365_day'))
        cube = iris.cube.Cube(data, dim_coords_and_dims=[(time, 0)])

        result = time_average(cube)
        expected = np.array([1.])
        self.assertArrayEqual(result.data, expected)


if __name__ == '__main__':
    unittest.main()
