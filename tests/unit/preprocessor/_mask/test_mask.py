"""Unit test for the :func:`esmvalcore.preprocessor._mask` function."""

import unittest

import numpy as np
from numpy.testing import assert_array_equal, assert_equal

import iris
import tests
from cf_units import Unit
from esmvalcore.preprocessor._mask import (
    _apply_fx_mask, _check_dims,
    count_spells, _get_fx_mask,
    mask_above_threshold, mask_below_threshold, mask_inside_range,
    mask_outside_range)


class Test(tests.Test):
    """Test class for _mask."""

    def setUp(self):
        """Prepare tests."""
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
        self.time_cube = iris.cube.Cube(np.arange(1, 25),
                                        var_name='co2', units='J')
        self.time_cube.add_dim_coord(
            iris.coords.DimCoord(
                np.arange(15., 720., 30.),
                standard_name='time',
                units=Unit('days since 1950-01-01 00:00:00',
                           calendar='gregorian')), 0)
        self.fx_data = np.array([20., 60., 50.])

    def test_apply_fx_mask(self):
        """Test _apply_fx_mask func."""
        dummy_fx_mask = np.ma.array([33., 22., 33.], mask=[True, False, True],
                                    fill_value=1e+20)
        app_mask = _apply_fx_mask(dummy_fx_mask,
                                  self.time_cube.data[0:3].astype('float64'))
        fixed_mask = np.ma.array(self.time_cube.data[0:3].astype('float64'),
                                 mask=dummy_fx_mask.mask)
        assert_array_equal(fixed_mask, app_mask)

    def test_check_dims(self):
        """Test _check_dims func."""
        malformed_cube = self.arr[0]
        assert_equal(True, _check_dims(self.arr, self.arr))
        assert_equal(False, _check_dims(self.arr, malformed_cube))

    def test_count_spells(self):
        """Test count_spells func."""
        ref_spells = count_spells(self.time_cube.data, -1000., 0, 1)
        assert_equal(24, ref_spells)
        ref_spells = count_spells(self.time_cube.data, -1000., 0, 2)
        assert_equal(12, ref_spells)

    def test_get_fx_mask(self):
        """Test _get_fx_mask func."""
        # sftlf: land. sea
        computed = _get_fx_mask(self.fx_data, 'land', 'sftlf')
        expected = np.array([False, True, False])
        assert_array_equal(expected, computed)
        computed = _get_fx_mask(self.fx_data, 'sea', 'sftlf')
        expected = np.array([True, False, True])
        assert_array_equal(expected, computed)
        # sftof: land, sea
        computed = _get_fx_mask(self.fx_data, 'land', 'sftof')
        expected = np.array([True, False, False])
        assert_array_equal(expected, computed)
        computed = _get_fx_mask(self.fx_data, 'sea', 'sftof')
        expected = np.array([False, True, True])
        assert_array_equal(expected, computed)
        # sftgif: ice, landsea
        computed = _get_fx_mask(self.fx_data, 'ice', 'sftgif')
        expected = np.array([False, True, False])
        assert_array_equal(expected, computed)
        computed = _get_fx_mask(self.fx_data, 'landsea', 'sftgif')
        expected = np.array([True, False, True])
        assert_array_equal(expected, computed)

    def test_mask_above_threshold(self):
        """Test to mask above a threshold."""
        result = mask_above_threshold(self.arr, 1.5)
        expected = np.ma.array(self.data2, mask=[[False, False], [True, True]])
        assert_array_equal(result.data, expected)

    def test_mask_below_threshold(self):
        """Test to mask below a threshold."""
        result = mask_below_threshold(self.arr, 1.5)
        expected = np.ma.array(self.data2, mask=[[True, True], [False, False]])
        assert_array_equal(result.data, expected)

    def test_mask_inside_range(self):
        """Test to mask inside a range."""
        result = mask_inside_range(self.arr, 0.5, 2.5)
        expected = np.ma.array(self.data2, mask=[[False, True], [True, False]])
        assert_array_equal(result.data, expected)

    def test_mask_outside_range(self):
        """Test to mask outside a range."""
        result = mask_outside_range(self.arr, 0.5, 2.5)
        expected = np.ma.array(self.data2, mask=[[True, False], [False, True]])
        assert_array_equal(result.data, expected)


if __name__ == '__main__':
    unittest.main()
