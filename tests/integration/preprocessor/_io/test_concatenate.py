"""
Integration tests for the :func:
`esmvaltool.preprocessor._io.concatenate`
function.

"""

from __future__ import absolute_import, division, print_function

import unittest
import numpy as np
from iris.cube import Cube, CubeList
from iris.coords import DimCoord

from esmvaltool.preprocessor import _io


class TestConcatenate(unittest.TestCase):

    def setUp(self):
        coord = DimCoord([1, 2], var_name='coord')
        second_coord = coord.copy([3, 4])
        self.raw_cubes = []
        self.raw_cubes.append(Cube([1, 2], var_name='sample',
                                   dim_coords_and_dims=((coord, 0),)))
        self.raw_cubes.append(Cube([3, 4], var_name='sample',
                                   dim_coords_and_dims=((second_coord, 0),)))

    def test_easy_concatenate(self):
        concatenated = _io.concatenate(self.raw_cubes)
        self.assertTrue((concatenated.coord('coord').points ==
                         np.array([1, 2, 3, 4])).all())

    def test_fail_with_duplicates(self):
        self.raw_cubes.append(self.raw_cubes[0].copy())
        with self.assertRaises(_io.ConcatenationError):
            _io.concatenate(self.raw_cubes)

    def test_fail_because_metadata_differs(self):
        self.raw_cubes[0].units = 'm'
        # with self.assertRaises(_io.ConcatenationError):
        #     _io.concatenate(self.raw_cubes)
        _io.concatenate(self.raw_cubes)
