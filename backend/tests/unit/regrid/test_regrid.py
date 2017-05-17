"""
Unit tests for the :func:`esmvaltool.backend.regrid.regrid` function.

"""

from __future__ import (absolute_import, division, print_function)
from six.moves import (filter, input, map, range, zip)  # noqa
import six

import iris
import mock
import unittest

import backend.tests as tests
from backend.regrid import _cache, horizontal_schemes, regrid


class Test(tests.Test):
    def _check(self, tgt_grid, scheme, spec=False):
        if spec:
            self.assertIn(tgt_grid, _cache)
            self.assertEqual(_cache[tgt_grid], tgt_grid)

        self.mock_regrid.assert_called_once_with(tgt_grid,
                                                 horizontal_schemes[scheme])

        # Reset the mocks to enable multiple calls per test-case.
        for mocker in self.mocks:
            mocker.reset_mock()

    def setUp(self):
        self.cube = iris.cube.Cube(0)
        self.regrid_schemes = ['linear', 'nearest', 'area_weighted',
                               'unstructured_nearest']
        self.regridded_cube = mock.sentinel.regridded_cube
        self.mock_regrid = self.patch('iris.cube.Cube.regrid',
                                      return_value=self.regridded_cube)
        self.mock_stock = self.patch('backend.regrid._stock_cube',
                                     side_effect=lambda arg: arg)
        self.mocks = [self.mock_regrid, self.mock_stock]

    def test_nop(self):
        cube = mock.sentinel.cube
        result = regrid(cube, None, None)
        self.assertEqual(result, cube)

    def test_invalid_tgt_grid__None(self):
        dummy = mock.sentinel.dummy
        emsg = 'A target grid must be specified'
        with self.assertRaisesRegexp(ValueError, emsg):
            regrid(dummy, None, dummy)

    def test_invalid_tgt_grid__unknown(self):
        dummy = mock.sentinel.dummy
        scheme = 'linear'
        emsg = 'Expecting a cube or cell-specification'
        with self.assertRaisesRegexp(ValueError, emsg):
            regrid(self.cube, dummy, scheme)

    def test_invalid_scheme__None(self):
        dummy = mock.sentinel.dummy
        emsg = 'A scheme must be specified'
        with self.assertRaisesRegexp(ValueError, emsg):
            regrid(dummy, dummy, None)

    def test_invalid_scheme__unknown(self):
        dummy = mock.sentinel.dummy
        emsg = 'Unknown regridding scheme'
        with self.assertRaisesRegexp(ValueError, emsg):
            regrid(dummy, dummy, 'wibble')

    def test_horizontal_schemes(self):
        self.assertEqual(horizontal_schemes.viewkeys(),
                         set(self.regrid_schemes))

    def test_regrid__horizontal_schemes(self):
        tgt_grid = mock.Mock(spec=iris.cube.Cube)
        for scheme in self.regrid_schemes:
            result = regrid(self.cube, tgt_grid, scheme)
            self.assertEqual(result, self.regridded_cube)
            self._check(tgt_grid, scheme)

    def test_regrid__cell_specification(self):
        specs = ['one', 'two', 'three', 'four', 'five']
        scheme = 'linear'
        for spec in specs:
            result = regrid(self.cube, spec, scheme)
            self.assertEqual(result, self.regridded_cube)
            self._check(spec, scheme, spec=True)
        self.assertEqual(_cache.viewkeys(), set(specs))


if __name__ == '__main__':
    unittest.main()
