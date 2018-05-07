"""
Unit tests for the :func:`esmvaltool.preprocessor.regrid._stock_cube`
function.

"""

from __future__ import absolute_import, division, print_function

import unittest

import iris
import mock
import numpy as np

import tests
from esmvaltool.preprocessor._regrid import _stock_cube as stock_cube
from esmvaltool.preprocessor._regrid import (
    _LAT_MAX, _LAT_MIN, _LAT_RANGE, _LON_MAX, _LON_MIN, _LON_RANGE)


class Test(tests.Test):
    def _check(self, dx, dy):
        # Generate the expected stock cube coordinate points.
        dx, dy = float(dx), float(dy)
        mid_dx, mid_dy = dx / 2, dy / 2
        expected_lat_points = np.linspace(_LAT_MIN + mid_dy, _LAT_MAX - mid_dy,
                                          _LAT_RANGE / dy)
        expected_lon_points = np.linspace(_LON_MIN + mid_dx, _LON_MAX - mid_dx,
                                          _LON_RANGE / dx)

        # Check the stock cube coordinates.
        self.assertEqual(self.mock_DimCoord.call_count, 2)
        call_lats, call_lons = self.mock_DimCoord.call_args_list

        # Check the latitude coordinate creation.
        [args], kwargs = call_lats
        self.assertArrayEqual(args, expected_lat_points)
        expected_lat_kwargs = dict(
            standard_name='latitude', units='degrees_north', var_name='lat')
        self.assertEqual(kwargs, expected_lat_kwargs)

        # Check the longitude coordinate creation.
        [args], kwargs = call_lons
        self.assertArrayEqual(args, expected_lon_points)
        expected_lon_kwargs = dict(
            standard_name='longitude', units='degrees_east', var_name='lon')
        self.assertEqual(kwargs, expected_lon_kwargs)

        # Check that the coordinate guess_bounds method has been called.
        expected_calls = [mock.call.guess_bounds()] * 2
        self.assertEqual(self.mock_coord.mock_calls, expected_calls)

        # Check the stock cube creation.
        self.mock_Cube.assert_called_once()
        _, kwargs = self.mock_Cube.call_args
        spec = [(self.mock_coord, 0), (self.mock_coord, 1)]
        expected_cube_kwargs = dict(dim_coords_and_dims=spec)
        self.assertEqual(kwargs, expected_cube_kwargs)

        # Reset the mocks to enable multiple calls per test-case.
        for mocker in self.mocks:
            mocker.reset_mock()

    def setUp(self):
        self.Cube = mock.sentinel.Cube
        self.mock_Cube = self.patch('iris.cube.Cube', return_value=self.Cube)
        self.mock_coord = mock.Mock(spec=iris.coords.DimCoord)
        self.mock_DimCoord = self.patch(
            'iris.coords.DimCoord', return_value=self.mock_coord)
        self.mocks = [self.mock_Cube, self.mock_coord, self.mock_DimCoord]

    def test_invalid_cell_spec__alpha(self):
        emsg = 'Invalid MxN cell specification'
        with self.assertRaisesRegex(ValueError, emsg):
            stock_cube('Ax1')

    def test_invalid_cell_spec__separator(self):
        emsg = 'Invalid MxN cell specification'
        with self.assertRaisesRegex(ValueError, emsg):
            stock_cube('1y1')

    def test_invalid_cell_spec__longitude(self):
        emsg = 'Invalid longitude delta in MxN cell specification'
        with self.assertRaisesRegex(ValueError, emsg):
            stock_cube('1.3x1')

    def test_invalid_cell_spec__latitude(self):
        emsg = 'Invalid latitude delta in MxN cell specification'
        with self.assertRaisesRegex(ValueError, emsg):
            stock_cube('1x2.3')

    def test_specs(self):
        specs = ['0.5x0.5', '1x1', '2.5x2.5', '5x5', '10x10']
        for spec in specs:
            result = stock_cube(spec)
            self.assertEqual(result, self.Cube)
            self._check(*list(map(float, spec.split('x'))))


if __name__ == '__main__':
    unittest.main()
