"""
Integration tests for the :func:
`esmvalcore.preprocessor.regrid.get_cmor_levels`
function.

"""

import os
import tempfile
import unittest

import iris
import iris.coords
import iris.cube
import numpy as np

from esmvalcore.preprocessor import _regrid


class TestGetFileLevels(unittest.TestCase):
    def setUp(self):
        """Prepare the sample file for the test"""
        self.cube = iris.cube.Cube(np.ones([2, 2, 2]), var_name='var')
        self.cube.add_dim_coord(
            iris.coords.DimCoord(np.arange(0, 2), var_name='coord'), 0)

        self.cube.coord('coord').attributes['positive'] = 'up'
        iris.util.guess_coord_axis(self.cube.coord('coord'))
        descriptor, self.path = tempfile.mkstemp('.nc')
        os.close(descriptor)
        print(self.cube)
        iris.save(self.cube, self.path)

    def tearDown(self):
        """Remove the sample file for the test"""
        os.remove(self.path)

    def test_get_coord(self):
        self.assertListEqual(
            _regrid.get_reference_levels(
                self.path, 'project', 'dataset', 'short_name', 'output_dir'),
            [0., 1]
        )
