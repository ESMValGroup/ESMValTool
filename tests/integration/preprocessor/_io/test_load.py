"""Integration tests for :func:`esmvaltool.preprocessor._io.load`."""

import os
import tempfile
import unittest

import iris
import numpy as np
from iris.coords import DimCoord
from iris.cube import Cube

from esmvaltool.preprocessor._io import concatenate_callback, load


def _create_sample_cube():
    coord = DimCoord([1, 2], standard_name='latitude', units='degrees_north')
    cube = Cube([1, 2], var_name='sample', dim_coords_and_dims=((coord, 0), ))
    return cube


class TestLoad(unittest.TestCase):
    """Tests for :func:`esmvaltool.preprocessor.load`."""

    def setUp(self):
        """Start tests."""
        self.temp_files = []

    def tearDown(self):
        """Finish tests."""
        for temp_file in self.temp_files:
            os.remove(temp_file)

    def _save_cube(self, cube):
        descriptor, temp_file = tempfile.mkstemp('.nc')
        os.close(descriptor)
        iris.save(cube, temp_file)
        self.temp_files.append(temp_file)
        return temp_file

    def test_load(self):
        """Test loading multiple files."""
        cube = _create_sample_cube()
        temp_file = self._save_cube(cube)

        cubes = load(temp_file)
        cube = cubes[0]
        self.assertEqual(1, len(cubes))
        self.assertEqual(temp_file, cube.attributes['source_file'])
        self.assertTrue((cube.data == np.array([1, 2])).all())
        self.assertTrue((cube.coord('latitude').points == np.array([1,
                                                                    2])).all())

    def test_callback_remove_attributes(self):
        """Test callback remove unwanted attributes."""
        attributes = ('history', 'creation_date', 'tracking_id')
        for _ in range(2):
            cube = _create_sample_cube()
            for attr in attributes:
                cube.attributes[attr] = attr
            self._save_cube(cube)
        for temp_file in self.temp_files:
            cubes = load(temp_file, callback=concatenate_callback)
            cube = cubes[0]
            self.assertEqual(1, len(cubes))
            self.assertTrue((cube.data == np.array([1, 2])).all())
            self.assertTrue(
                (cube.coord('latitude').points == np.array([1, 2])).all())
            for attr in attributes:
                self.assertTrue(attr not in cube.attributes)

    def test_callback_fix_lat_units(self):
        """Test callback for fixing units."""
        cube = _create_sample_cube()
        temp_file = self._save_cube(cube)

        cubes = load(temp_file, callback=concatenate_callback)
        cube = cubes[0]
        self.assertEqual(1, len(cubes))
        self.assertTrue((cube.data == np.array([1, 2])).all())
        self.assertTrue((cube.coord('latitude').points == np.array([1,
                                                                    2])).all())
        self.assertEquals(cube.coord('latitude').units, 'degrees_north')
