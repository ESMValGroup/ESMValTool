"""Integration tests for :func:`esmvaltool.preprocessor._io.concatenate`"""

from __future__ import absolute_import, division, print_function

import unittest
import os
import tempfile
import numpy as np
import iris
from iris.cube import Cube
from iris.coords import DimCoord

from esmvaltool.preprocessor import _io


class TestLoad(unittest.TestCase):
    """Tests for :func:`esmvaltool.preprocessor._io.concatenate`"""

    def setUp(self):
        self.temp_files = []

    def tearDown(self):
        for temp_file in self.temp_files:
            os.remove(temp_file)

    def _create_sample_cube(self):
        coord = DimCoord([1, 2], standard_name='latitude',
                         units='degrees_north')
        cube = Cube([1, 2], var_name='sample',
                    dim_coords_and_dims=((coord, 0),))
        return cube

    def _save_cube(self, cube):
        descriptor, temp_file = tempfile.mkstemp('.nc')
        os.close(descriptor)
        iris.save(cube, temp_file)
        self.temp_files.append(temp_file)

    def test_load_multiple(self):
        """Test loading multiple files"""
        for num in range(2):
            cube = self._create_sample_cube()
            self._save_cube(cube)

        list = _io.load_cubes(self.temp_files, 'filename', None)
        cube = list[0]
        self.assertTrue((cube.data == np.array([1, 2])).all())
        self.assertTrue((cube.coord('latitude').points ==
                         np.array([1, 2])).all())
        self.assertEquals(cube.attributes['_filename'], 'filename')

    def test_callback_remove_attributtes(self):
        """Test callback remove unwanted attributes"""
        attributtes = ('history', 'creation_date', 'tracking_id')
        for x in range(2):
            cube = self._create_sample_cube()
            for attr in attributtes:
                cube.attributes[attr] = attr
            self._save_cube(cube)

        cubes = _io.load_cubes(self.temp_files, 'filename', None,
                               callback=_io.concatenate_callback)
        cube = cubes[0]
        self.assertTrue((cube.data == np.array([1, 2])).all())
        self.assertTrue((cube.coord('latitude').points ==
                         np.array([1, 2])).all())
        self.assertEqual(cube.attributes['_filename'], 'filename')
        for attr in attributtes:
            self.assertTrue(attr not in cube.attributes)

    def test_callback_fix_lat_units(self):
        """Test callback for fixing units"""
        cube = self._create_sample_cube()
        self._save_cube(cube)

        list = _io.load_cubes(self.temp_files, 'filename', None,
                              callback=_io.concatenate_callback)
        cube = list[0]
        self.assertTrue((cube.data == np.array([1, 2])).all())
        self.assertTrue((cube.coord('latitude').points ==
                         np.array([1, 2])).all())
        self.assertEquals(cube.attributes['_filename'], 'filename')
        self.assertEquals(cube.coord('latitude').units, 'degrees_north')
