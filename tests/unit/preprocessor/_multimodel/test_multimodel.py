"""Unit test for :func:`esmvalcore.preprocessor._volume`."""

import unittest

import cftime
import iris
import numpy as np
from cf_units import Unit
from numpy.testing import assert_array_equal, assert_equal

import tests
from esmvalcore.preprocessor._multimodel import (_assemble_overlap_data,
                                                 _compute_statistic,
                                                 _datetime_to_int_days,
                                                 _get_time_offset,
                                                 _put_in_cube)


class Test(tests.Test):
    """Test class for preprocessor/_multimodel.py."""

    def setUp(self):
        """Prepare tests."""
        coord_sys = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        data2 = np.ma.ones((2, 3, 2, 2))
        data3 = np.ma.ones((4, 3, 2, 2))
        mask3 = np.full((4, 3, 2, 2), False)
        mask3[0, 0, 0, 0] = True
        data3 = np.ma.array(data3, mask=mask3)

        time = iris.coords.DimCoord([15, 45],
                                    standard_name='time',
                                    bounds=[[1., 30.], [30., 60.]],
                                    units=Unit(
                                        'days since 1950-01-01',
                                        calendar='gregorian'))
        time2 = iris.coords.DimCoord([1., 2., 3., 4.],
                                     standard_name='time',
                                     bounds=[
                                         [0.5, 1.5],
                                         [1.5, 2.5],
                                         [2.5, 3.5],
                                         [3.5, 4.5], ],
                                     units=Unit(
                                         'days since 1950-01-01',
                                         calendar='gregorian'))

        zcoord = iris.coords.DimCoord([0.5, 5., 50.],
                                      standard_name='air_pressure',
                                      long_name='air_pressure',
                                      bounds=[[0., 2.5], [2.5, 25.],
                                              [25., 250.]],
                                      units='m',
                                      attributes={'positive': 'down'})
        lons = iris.coords.DimCoord([1.5, 2.5],
                                    standard_name='longitude',
                                    bounds=[[1., 2.], [2., 3.]],
                                    units='degrees_east',
                                    coord_system=coord_sys)
        lats = iris.coords.DimCoord([1.5, 2.5],
                                    standard_name='latitude',
                                    bounds=[[1., 2.], [2., 3.]],
                                    units='degrees_north',
                                    coord_system=coord_sys)

        coords_spec4 = [(time, 0), (zcoord, 1), (lats, 2), (lons, 3)]
        self.cube1 = iris.cube.Cube(data2, dim_coords_and_dims=coords_spec4)

        coords_spec5 = [(time2, 0), (zcoord, 1), (lats, 2), (lons, 3)]
        self.cube2 = iris.cube.Cube(data3, dim_coords_and_dims=coords_spec5)

    def test_get_time_offset(self):
        """Test time unit."""
        result = _get_time_offset("days since 1950-01-01")
        expected = cftime.real_datetime(1950, 1, 1, 0, 0)
        assert_equal(result, expected)

    def test_compute_statistic(self):
        """Test statistic."""
        datas = [self.cube1.data[0], self.cube2.data[0]]
        stat_mean = _compute_statistic(datas, "mean")
        stat_median = _compute_statistic(datas, "median")
        expected_mean = np.ma.ones((3, 2, 2))
        expected_median = np.ma.ones((3, 2, 2))
        assert_array_equal(stat_mean, expected_mean)
        assert_array_equal(stat_median, expected_median)

    def test_put_in_cube(self):
        """Test put in cube."""
        cube_data = np.ma.ones((2, 3, 2, 2))
        stat_cube = _put_in_cube(self.cube1, cube_data, "mean", t_axis=None)
        assert_array_equal(stat_cube.data, self.cube1.data)

    def test_datetime_to_int_days(self):
        """Test _datetime_to_int_days."""
        computed_dats = _datetime_to_int_days(self.cube1)
        expected_dats = [0, 31]
        assert_array_equal(computed_dats, expected_dats)

    def test_assemble_overlap_data(self):
        """Test overlap data."""
        comp_ovlap_mean = _assemble_overlap_data([self.cube1, self.cube1],
                                                 [0, 31], "mean")
        expected_ovlap_mean = np.ma.ones((2, 3, 2, 2))
        assert_array_equal(comp_ovlap_mean.data, expected_ovlap_mean)


if __name__ == '__main__':
    unittest.main()
