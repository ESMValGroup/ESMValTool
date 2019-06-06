"""Unit test for :func:`esmvalcore.preprocessor._volume`."""

import unittest

import iris
import numpy as np
from numpy.testing import assert_array_equal, assert_equal
import cftime
from cf_units import Unit

import tests
from esmvalcore.preprocessor._multimodel import (_get_time_offset,
                                                 _compute_statistic,
                                                 _put_in_cube,
                                                 _datetime_to_int_days,
                                                 _assemble_overlap_data)


class Test(tests.Test):
    """Test class for _volume_pp"""

    def setUp(self):
        """Prepare tests"""
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
                                         [3.5, 4.5],
                                     ],
                                     units=Unit(
                                         'days since 1950-01-01',
                                         calendar='gregorian'))

        zcoord = iris.coords.DimCoord([0.5, 5., 50.],
                                      long_name='zcoord',
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

        # allow iris to figure out the axis='z' coordinate
        iris.util.guess_coord_axis(self.cube1.coord('zcoord'))
        iris.util.guess_coord_axis(self.cube2.coord('zcoord'))

    def test_get_time_offset(self):
        """Test time unit."""
        result = _get_time_offset("days since 1950-01-01")
        expected = cftime.real_datetime(1950, 1, 1, 0, 0)
        assert_equal(result, expected)

    def test_compute_statistic(self):
        """Test time unit."""
        datas = [self.cube1.data[0], self.cube2.data[0]]
        stat_mean = _compute_statistic(datas, "mean")
        stat_median = _compute_statistic(datas, "median")
        expected_mean = np.ma.ones((3, 2, 2))
        expected_median = np.ma.ones((3, 2, 2))
        assert_array_equal(stat_mean, expected_mean)
        assert_array_equal(stat_median, expected_median)



if __name__ == '__main__':
    unittest.main()
