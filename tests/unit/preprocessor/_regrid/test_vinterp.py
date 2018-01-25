"""
Unit tests for the :func:`esmvaltool.preprocessor.regrid.vinterp`
function.

"""

from __future__ import absolute_import, division, print_function

import unittest

import mock
import numpy as np
from numpy import ma

import tests
from esmvaltool.preprocessor._regrid import _MDI, vertical_schemes, vinterp
from tests.unit.preprocessor._regrid import _make_cube, _make_vcoord


class Test(tests.Test):
    def setUp(self):
        self.shape = (3, 2, 1)
        self.z = self.shape[0]
        self.dtype = np.dtype('int8')
        data = np.arange(
            np.prod(self.shape), dtype=self.dtype).reshape(self.shape)
        self.cube = _make_cube(data, dtype=self.dtype)
        self.created_cube = mock.sentinel.created_cube
        self.mock_create_cube = self.patch(
            'esmvaltool.preprocessor._regrid._create_cube',
            return_value=self.created_cube)
        self.vinterp_schemes = ['linear', 'nearest']

    def test_nop(self):
        cube = mock.sentinel.cube
        result = vinterp(cube, None, None)
        self.assertEqual(result, cube)

    def test_invalid_levels__None(self):
        emsg = 'Target levels must be specified'
        with self.assertRaisesRegex(ValueError, emsg):
            vinterp(self.cube, None, 'linear')

    def test_invalid_scheme__None(self):
        levels = mock.sentinel.levels
        emsg = 'A scheme must be specified'
        with self.assertRaisesRegex(ValueError, emsg):
            vinterp(self.cube, levels, None)

    def test_invalid_scheme__unknown(self):
        levels = mock.sentinel.levels
        scheme = mock.sentinel.scheme
        emsg = 'Unknown vertical interpolation scheme'
        with self.assertRaisesRegex(ValueError, emsg):
            vinterp(self.cube, levels, scheme)

    def test_vertical_schemes(self):
        self.assertEqual(set(vertical_schemes), set(self.vinterp_schemes))

    def test_nop__levels_match(self):
        vcoord = _make_vcoord(self.z, dtype=self.dtype)
        self.assertEqual(self.cube.coord(axis='z', dim_coords=True), vcoord)
        levels = vcoord.points
        result = vinterp(self.cube, levels, 'linear')
        self.assertEqual(id(result), id(self.cube))
        self.assertEqual(result, self.cube)

    def test_extraction(self):
        levels = [0, 2]
        result = vinterp(self.cube, levels, 'linear')
        data = np.array([0, 1, 4, 5], dtype=self.dtype).reshape(2, 2, 1)
        expected = _make_cube(
            data, aux_coord=False, dim_coord=False, dtype=self.dtype)
        coord = self.cube.coord('Pressure Slice').copy()
        expected.add_aux_coord(coord[levels], (0, 1))
        coord = self.cube.coord('air_pressure').copy()
        expected.add_dim_coord(coord[levels], 0)
        self.assertEqual(result, expected)

    def test_extraction__failure(self):
        levels = [0, 2]
        with mock.patch('iris.cube.Cube.extract', return_value=None):
            emsg = 'Failed to extract levels'
            with self.assertRaisesRegex(ValueError, emsg):
                vinterp(self.cube, levels, 'linear')

    def test_interpolation(self):
        new_data = np.array(True)
        levels = np.array([0.5, 1.5])
        scheme = 'linear'
        with mock.patch(
                'stratify.interpolate', return_value=new_data) as mocker:
            result = vinterp(self.cube, levels, scheme)
            self.assertEqual(result, self.created_cube)
            args, kwargs = mocker.call_args
            # Check the stratify.interpolate args ...
            self.assertEqual(len(args), 3)
            self.assertArrayEqual(args[0], levels)
            pts = self.cube.coord(axis='z', dim_coords=True).points
            src_levels_broadcast = np.broadcast_to(
                pts.reshape(self.z, 1, 1), self.cube.shape)
            self.assertArrayEqual(args[1], src_levels_broadcast)
            self.assertArrayEqual(args[2], self.cube.data)
            # Check the stratify.interpolate kwargs ...
            self.assertEqual(kwargs,
                             dict(
                                 axis=0,
                                 interpolation=scheme,
                                 extrapolation='nan'))
        args, kwargs = self.mock_create_cube.call_args
        # Check the _create_cube args ...
        self.assertEqual(len(args), 3)
        self.assertArrayEqual(args[0], self.cube)
        self.assertArrayEqual(args[1], new_data)
        self.assertArrayEqual(args[2], levels)
        # Check the _create_cube kwargs ...
        self.assertEqual(kwargs, dict())

    def test_interpolation__extrapolated_NaN_filling(self):
        new_data = np.array([0, np.nan])
        levels = [0.5, 1.5]
        scheme = 'nearest'
        with mock.patch(
                'stratify.interpolate', return_value=new_data) as mocker:
            result = vinterp(self.cube, levels, scheme)
            self.assertEqual(result, self.created_cube)
            args, kwargs = mocker.call_args
            # Check the stratify.interpolate args ...
            self.assertEqual(len(args), 3)
            self.assertArrayEqual(args[0], levels)
            pts = self.cube.coord(axis='z', dim_coords=True).points
            src_levels_broadcast = np.broadcast_to(
                pts.reshape(self.z, 1, 1), self.cube.shape)
            self.assertArrayEqual(args[1], src_levels_broadcast)
            self.assertArrayEqual(args[2], self.cube.data)
            # Check the stratify.interpolate kwargs ...
            self.assertEqual(kwargs,
                             dict(
                                 axis=0,
                                 interpolation=scheme,
                                 extrapolation='nan'))
        args, kwargs = self.mock_create_cube.call_args
        # Check the _create_cube args ...
        self.assertEqual(len(args), 3)
        self.assertArrayEqual(args[0], self.cube)
        new_data[np.isnan(new_data)] = _MDI
        self.assertArrayEqual(args[1], new_data)
        self.assertArrayEqual(args[2], levels)
        # Check the _create_cube kwargs ...
        self.assertEqual(kwargs, dict())

    def test_interpolation__masked(self):
        new_data = np.empty(self.shape)
        levels = np.array([0.5, 1.5])
        scheme = 'linear'
        mask = [[[1], [0]], [[1], [0]], [[1], [0]]]
        masked = ma.empty(self.shape)
        masked.mask = mask
        cube = _make_cube(masked, dtype=self.dtype)
        with mock.patch(
                'stratify.interpolate', return_value=new_data) as mocker:
            result = vinterp(cube, levels, scheme)
            self.assertEqual(result, self.created_cube)
            args, kwargs = mocker.call_args
            # Check the stratify.interpolate args ...
            self.assertEqual(len(args), 3)
            self.assertArrayEqual(args[0], levels)
            pts = cube.coord(axis='z', dim_coords=True).points
            src_levels_broadcast = np.broadcast_to(
                pts.reshape(self.z, 1, 1), cube.shape)
            self.assertArrayEqual(args[1], src_levels_broadcast)
            self.assertArrayEqual(args[2], cube.data)
            # Check the stratify.interpolate kwargs ...
            self.assertEqual(kwargs,
                             dict(
                                 axis=0,
                                 interpolation=scheme,
                                 extrapolation='nan'))
        args, kwargs = self.mock_create_cube.call_args
        # Check the _create_cube args ...
        self.assertEqual(len(args), 3)
        self.assertArrayEqual(args[0], cube)
        self.assertArrayEqual(args[1], new_data)
        self.assertTrue(ma.isMaskedArray(args[1]))
        self.assertArrayEqual(args[1].mask, mask)
        self.assertArrayEqual(args[2], levels)
        # Check the _create_cube kwargs ...
        self.assertEqual(kwargs, dict())


if __name__ == '__main__':
    unittest.main()
