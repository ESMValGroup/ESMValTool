"""Unit tests for the :func:`esmvalcore.preprocessor._time` module."""

import unittest

import iris
import iris.coord_categorisation
import iris.coords
import numpy as np
import pytest
from cf_units import Unit
from iris.cube import Cube
from numpy.testing import assert_array_equal

import tests
from esmvalcore.preprocessor._time import (annual_mean, extract_month,
                                           extract_season, extract_time,
                                           regrid_time, time_average)


def _create_sample_cube():
    cube = Cube(np.arange(1, 25), var_name='co2', units='J')
    cube.add_dim_coord(
        iris.coords.DimCoord(
            np.arange(15., 720., 30.),
            standard_name='time',
            units=Unit('days since 1950-01-01 00:00:00', calendar='gregorian'),
        ),
        0,
    )
    iris.coord_categorisation.add_month_number(cube, 'time')
    return cube


def add_auxiliary_coordinate(cubeList):
    """Add AuxCoords to cubes in cubeList."""
    for cube in cubeList:
        iris.coord_categorisation.add_day_of_month(cube, cube.coord('time'))
        iris.coord_categorisation.add_day_of_year(cube, cube.coord('time'))


class TestExtractMonth(tests.Test):
    """Tests for extract_month."""

    def setUp(self):
        """Prepare tests"""
        self.cube = _create_sample_cube()

    def test_get_january(self):
        """Test january extraction"""
        sliced = extract_month(self.cube, 1)
        print(sliced)
        assert_array_equal(
            np.array([1, 1]),
            sliced.coord('month_number').points)


class TestTimeSlice(tests.Test):
    """Tests for extract_time."""

    def setUp(self):
        """Prepare tests"""
        self.cube = _create_sample_cube()

    def test_extract_time(self):
        """Test extract_time."""
        sliced = extract_time(self.cube, 1950, 1, 1, 1950, 12, 31)
        print(sliced)
        assert_array_equal(
            np.arange(1, 13, 1),
            sliced.coord('month_number').points)

    def test_extract_time_no_slice(self):
        """Test fail of extract_time."""
        with self.assertRaises(ValueError):
            extract_time(self.cube, 2200, 1, 1, 2200, 12, 31)

    def test_extract_time_one_time(self):
        """Test extract_time with one time step."""
        cube = _create_sample_cube()
        cube = cube.collapsed('time', iris.analysis.MEAN)
        sliced = extract_time(cube, 1950, 1, 1, 1952, 12, 31)
        print(sliced)
        assert_array_equal(np.array([360.]), sliced.coord('time').points)

    def test_extract_time_no_time(self):
        """Test extract_time with no time step."""
        cube = _create_sample_cube()[0]
        sliced = extract_time(cube, 1950, 1, 1, 1950, 12, 31)
        print('sliced', sliced, sliced.shape)
        print('cube', cube, cube.shape)
        assert cube == sliced


class TestExtractSeason(tests.Test):
    """Tests for extract_season."""

    def setUp(self):
        """Prepare tests"""
        self.cube = _create_sample_cube()

    def test_get_djf(self):
        """Test function for winter"""
        sliced = extract_season(self.cube, 'djf')
        print(sliced)
        assert_array_equal(
            np.array([1, 2, 12, 1, 2, 12]),
            sliced.coord('month_number').points)

    def test_get_djf_caps(self):
        """Test function works when season specified in caps"""
        sliced = extract_season(self.cube, 'DJF')
        print(sliced)
        assert_array_equal(
            np.array([1, 2, 12, 1, 2, 12]),
            sliced.coord('month_number').points)

    def test_get_mam(self):
        """Test function for spring"""
        sliced = extract_season(self.cube, 'mam')
        print(sliced)
        assert_array_equal(
            np.array([3, 4, 5, 3, 4, 5]),
            sliced.coord('month_number').points)

    def test_get_jja(self):
        """Test function for summer"""
        sliced = extract_season(self.cube, 'jja')
        print(sliced)
        assert_array_equal(
            np.array([6, 7, 8, 6, 7, 8]),
            sliced.coord('month_number').points)

    def test_get_son(self):
        """Test function for summer"""
        sliced = extract_season(self.cube, 'son')
        print(sliced)
        assert_array_equal(
            np.array([9, 10, 11, 9, 10, 11]),
            sliced.coord('month_number').points)


class TestTimeAverage(tests.Test):
    """Test class for the :func:`esmvalcore.preprocessor._time_pp` module"""

    def test_time_average(self):
        """Test for time average of a 1D field."""
        data = np.ones((3))
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
        assert_array_equal(result.data, expected)

    def test_time_average_uneven(self):
        """Test for time average of a 1D field with uneven time boundaries."""
        data = np.array([1., 5.])
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
        assert_array_equal(result.data, expected)

    def test_time_average_365_day(self):
        """Test for time avg of a realisitc time axis and 365 day calendar"""
        data = np.ones((6, ))
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
        assert_array_equal(result.data, expected)


class TestRegridTimeMonthly(tests.Test):
    """Tests for regrid_time with monthly frequency."""

    def setUp(self):
        """Prepare tests"""
        self.cube_1 = _create_sample_cube()
        self.cube_2 = _create_sample_cube()
        self.cube_2.data = self.cube_2.data * 2.
        self.cube_2.remove_coord('time')
        self.cube_2.add_dim_coord(
            iris.coords.DimCoord(
                np.arange(14., 719., 30.),
                standard_name='time',
                units=Unit(
                    'days since 1950-01-01 00:00:00', calendar='360_day'),
            ),
            0,
        )
        add_auxiliary_coordinate([self.cube_1, self.cube_2])

    def test_regrid_time_mon(self):
        """Test changes to cubes."""
        # test monthly
        newcube_1 = regrid_time(self.cube_1, frequency='mon')
        newcube_2 = regrid_time(self.cube_2, frequency='mon')
        # no changes to core data
        assert_array_equal(newcube_1.data, self.cube_1.data)
        assert_array_equal(newcube_2.data, self.cube_2.data)
        # no changes to number of coords and aux_coords
        assert len(newcube_1.coords()) == len(self.cube_1.coords())
        assert len(newcube_1.aux_coords) == len(self.cube_1.aux_coords)
        # test difference; also diff is zero
        expected = self.cube_1.data
        diff_cube = newcube_2 - newcube_1
        assert_array_equal(diff_cube.data, expected)


class TestRegridTimeDaily(tests.Test):
    """Tests for regrid_time with daily frequency."""

    def setUp(self):
        """Prepare tests"""
        self.cube_1 = _create_sample_cube()
        self.cube_2 = _create_sample_cube()
        self.cube_2.data = self.cube_2.data * 2.
        self.cube_1.remove_coord('time')
        self.cube_2.remove_coord('time')
        self.cube_1.add_dim_coord(
            iris.coords.DimCoord(
                np.arange(14. * 24. + 6., 38. * 24. + 6., 24.),
                standard_name='time',
                units=Unit(
                    'hours since 1950-01-01 00:00:00', calendar='360_day'),
            ),
            0,
        )
        self.cube_2.add_dim_coord(
            iris.coords.DimCoord(
                np.arange(14. * 24. + 3., 38. * 24. + 3., 24.),
                standard_name='time',
                units=Unit(
                    'hours since 1950-01-01 00:00:00', calendar='360_day'),
            ),
            0,
        )
        add_auxiliary_coordinate([self.cube_1, self.cube_2])

    def test_regrid_time_day(self):
        """Test changes to cubes."""
        # test daily
        newcube_1 = regrid_time(self.cube_1, frequency='day')
        newcube_2 = regrid_time(self.cube_2, frequency='day')
        # no changes to core data
        self.assertArrayEqual(newcube_1.data, self.cube_1.data)
        self.assertArrayEqual(newcube_2.data, self.cube_2.data)
        # no changes to number of coords and aux_coords
        assert len(newcube_1.coords()) == len(self.cube_1.coords())
        assert len(newcube_1.aux_coords) == len(self.cube_1.aux_coords)
        # test difference; also diff is zero
        expected = self.cube_1.data
        diff_cube = newcube_2 - newcube_1
        self.assertArrayEqual(diff_cube.data, expected)


def make_time_series(number_years=2):
    """Make a cube with time only dimension."""
    times = np.array([i * 30 + 15 for i in range(0, 12 * number_years, 1)])
    bounds = np.array([i * 30 for i in range(0, 12 * number_years + 1, 1)])
    bounds = np.array(
        [[bnd, bounds[index + 1]] for index, bnd in enumerate(bounds[:-1])])
    data = np.ones_like(times)
    time = iris.coords.DimCoord(
        times,
        bounds=bounds,
        standard_name='time',
        units=Unit('days since 1950-01-01', calendar='360_day'))
    cube = iris.cube.Cube(data, dim_coords_and_dims=[(time, 0)])
    return cube


@pytest.mark.parametrize('existing_coord', [True, False])
def test_annual_average(existing_coord):
    """Test for annual average."""
    cube = make_time_series(number_years=2)
    if existing_coord:
        iris.coord_categorisation.add_year(cube, 'time')

    result = annual_mean(cube, decadal=False)
    expected = np.array([1., 1.])
    assert_array_equal(result.data, expected)
    expected_time = np.array([180., 540.])
    assert_array_equal(result.coord('time').points, expected_time)


@pytest.mark.parametrize('existing_coord', [True, False])
def test_decadal_average(existing_coord):
    """Test for decadal average."""
    cube = make_time_series(number_years=20)
    if existing_coord:

        def get_decade(coord, value):
            """Callback function to get decades from cube."""
            date = coord.units.num2date(value)
            return date.year - date.year % 10

        iris.coord_categorisation.add_categorised_coord(
            cube, 'decade', 'time', get_decade)

    result = annual_mean(cube, decadal=True)
    expected = np.array([1., 1.])
    assert_array_equal(result.data, expected)
    expected_time = np.array([1800., 5400.])
    assert_array_equal(result.coord('time').points, expected_time)


if __name__ == '__main__':
    unittest.main()
