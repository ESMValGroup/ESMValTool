"""
Integration tests for the :func:`esmvaltool.preprocessor.regrid.extract_levels`
function.

"""

import unittest

import iris
import numpy as np

import tests
from esmvaltool.preprocessor._regrid import _MDI, extract_levels
from tests.unit.preprocessor._regrid import _make_cube, _make_vcoord


class Test(tests.Test):
    def setUp(self):
        shape = (3, 2, 2)
        self.z = shape[0]
        data = np.arange(np.prod(shape)).reshape(shape)
        cubes = iris.cube.CubeList()
        # Create first realization cube.
        cube = _make_cube(data)
        coord = iris.coords.DimCoord(0, standard_name='realization')
        cube.add_aux_coord(coord)
        cubes.append(cube)
        # Create second realization cube.
        cube = _make_cube(data + np.prod(shape))
        coord = iris.coords.DimCoord(1, standard_name='realization')
        cube.add_aux_coord(coord)
        cubes.append(cube)
        # Create a 4d synthetic test cube.
        self.cube = cubes.merge_cube()
        coord = self.cube.coord(axis='z', dim_coords=True)
        self.shape = list(self.cube.shape)
        [self.z_dim] = self.cube.coord_dims(coord)

    def test_nop__levels_match(self):
        vcoord = _make_vcoord(self.z)
        self.assertEqual(self.cube.coord(axis='z', dim_coords=True), vcoord)
        levels = vcoord.points
        result = extract_levels(self.cube, levels, 'linear')
        self.assertEqual(result, self.cube)
        self.assertEqual(id(result), id(self.cube))

    def test_interpolation__linear(self):
        levels = [0.5, 1.5]
        scheme = 'linear'
        result = extract_levels(self.cube, levels, scheme)
        expected = np.array([[[[2., 3.], [4., 5.]], [[6., 7.], [8., 9.]]],
                             [[[14., 15.], [16., 17.]], [[18., 19.],
                                                         [20., 21.]]]])
        self.assertArrayEqual(result.data, expected)
        self.shape[self.z_dim] = len(levels)
        self.assertEqual(result.shape, tuple(self.shape))

    def test_interpolation__nearest(self):
        levels = [0.49, 1.51]
        scheme = 'nearest'
        result = extract_levels(self.cube, levels, scheme)
        expected = np.array([[[[0., 1.], [2., 3.]], [[8., 9.], [10., 11.]]],
                             [[[12., 13.], [14., 15.]], [[20., 21.],
                                                         [22., 23.]]]])
        self.assertArrayEqual(result.data, expected)
        self.shape[self.z_dim] = len(levels)
        self.assertEqual(result.shape, tuple(self.shape))

    def test_interpolation__extrapolated_NaN_filling(self):
        levels = [-10, 1, 2, 10]
        scheme = 'nearest'
        result = extract_levels(self.cube, levels, scheme)
        expected = np.array(
            [[[[_MDI, _MDI], [_MDI, _MDI]], [[4., 5.], [6., 7.]],
              [[8., 9.], [10., 11.]], [[_MDI, _MDI], [_MDI, _MDI]]],
             [[[_MDI, _MDI], [_MDI, _MDI]], [[16., 17.], [18., 19.]],
              [[20., 21.], [22., 23.]], [[_MDI, _MDI], [_MDI, _MDI]]]])
        self.assertArrayEqual(result.data, expected)
        self.shape[self.z_dim] = len(levels)
        self.assertEqual(result.shape, tuple(self.shape))

    def test_interpolation__scalar_collapse(self):
        level = 1
        scheme = 'nearest'
        result = extract_levels(self.cube, level, scheme)
        expected = np.array([[[4., 5.], [6., 7.]], [[16., 17.], [18., 19.]]])
        self.assertArrayEqual(result.data, expected)
        del self.shape[self.z_dim]
        self.assertEqual(result.shape, tuple(self.shape))


if __name__ == '__main__':
    unittest.main()
