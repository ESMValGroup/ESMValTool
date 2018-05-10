"""Integration tests for :func:`esmvaltool.preprocessor._io.save`"""

from __future__ import absolute_import, division, print_function

import unittest
import os
import tempfile
import numpy as np
import iris
from iris.cube import Cube
from iris.coords import DimCoord

from esmvaltool.preprocessor import _io


class TestSave(unittest.TestCase):
    """Tests for :func:`esmvaltool.preprocessor._io.save`"""

    def setUp(self):
        self.temp_files = []

    def tearDown(self):
        for temp_file in self.temp_files:
            if os.path.isfile(temp_file):
                os.remove(temp_file)

    def _create_sample_cube(self):
        coord = DimCoord(np.asarray([1, 2], np.single),
                         standard_name='latitude',
                         units='degrees_north')

        cube = Cube(np.asarray([1, 2], np.single),
                    var_name='sample',
                    units='1',
                    dim_coords_and_dims=((coord, 0),))

        descriptor, temp_file = tempfile.mkstemp('.nc')
        os.close(descriptor)
        cube.attributes['_filename'] = temp_file
        self.temp_files.append(temp_file)
        return cube

    def test_save(self):
        """Test save"""
        cube = self._create_sample_cube()
        paths = _io.save_cubes([cube])
        loaded_cube = iris.load_cube(paths[0])
        self.assertTrue((cube.data == loaded_cube.data).all())
        self.assertTrue((cube.coord('latitude').points ==
                         loaded_cube.coord('latitude').points).all())

    def test_save_debug(self):
        """Test save on debug mode"""
        cube = self._create_sample_cube()
        paths = _io.save_cubes([cube], debug=True)
        loaded_cube = iris.load_cube(paths[0])
        self.assertTrue((cube.data == loaded_cube.data).all())
        self.assertTrue((cube.coord('latitude').points ==
                         loaded_cube.coord('latitude').points).all())

    def test_fail_without_filename(self):
        """Test save fails if _filename is not added"""
        cube = self._create_sample_cube()
        del cube.attributes['_filename']
        with self.assertRaises(ValueError):
            _io.save_cubes([cube])
