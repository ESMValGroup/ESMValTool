"""
Integration tests for the :func:
`esmvaltool.preprocessor.regrid.get_cmor_levels`
function.

"""

from __future__ import absolute_import, division, print_function

import unittest
import iris
import iris.cube
import iris.coords
import numpy as np
import tempfile
import os

from esmvaltool.preprocessor import _regrid


class TestGetFileLevels(unittest.TestCase):

    def setUp(self):
        """Prepare the sample file for the test"""
        self.cube = iris.cube.Cube(np.ones([2, 2, 2]), var_name='var')
        self.cube.add_dim_coord(iris.coords.DimCoord(np.arange(0, 2),
                                                     var_name='coord'), 0)

        descriptor, self.path = tempfile.mkstemp('.nc')
        os.close(descriptor)
        iris.save(self.cube, self.path)

    def tearDown(self):
        """Remove the sample file for the test"""
        os.remove(self.path)

    def test_get_coord(self):
        self.assertListEqual(_regrid.get_reference_levels(self.path, 'coord'),
                             [0., 1])

    def test_bad_coord(self):
        with self.assertRaises(ValueError):
            _regrid.get_reference_levels(self.path, 'bad_coord')
