"""Integration tests for :func:`esmvaltool.preprocessor._io.save`"""

from __future__ import absolute_import, division, print_function

import unittest
import os
import tempfile
import numpy as np
import netCDF4
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
        lat = DimCoord(np.asarray([1, 2], np.single),
                       standard_name='latitude',
                       units='degrees_north')
        lon = DimCoord(np.asarray([1, 2], np.single),
                       standard_name='longitude',
                       units='degrees_east')
        time = DimCoord(np.asarray([1, 2], np.single),
                        standard_name='time',
                        units='days since 2000-1-1')

        cube = Cube(np.random.random_sample([2, 2, 2]),
                    var_name='sample',
                    units='1',
                    dim_coords_and_dims=((lat, 0), (lon, 1), (time, 2)))

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
        self._compare_cubes(cube, loaded_cube)

    def test_save_zlib(self):
        """Test save"""
        cube = self._create_sample_cube()
        paths = _io.save_cubes([cube], compress=True)
        loaded_cube = iris.load_cube(paths[0])
        self._compare_cubes(cube, loaded_cube)
        handler = netCDF4.Dataset(paths[0], 'r')
        sample_filters = handler.variables['sample'].filters()
        self.assertTrue(sample_filters['zlib'])
        self.assertTrue(sample_filters['shuffle'])
        self.assertEqual(sample_filters['complevel'], 4)
        handler.close()

    def test_save_debug(self):
        """Test save on debug mode"""
        cube = self._create_sample_cube()
        paths = _io.save_cubes([cube], debug=True)
        loaded_cube = iris.load_cube(paths[0])
        self._compare_cubes(cube, loaded_cube)

    def test_fail_without_filename(self):
        """Test save fails if _filename is not added"""
        cube = self._create_sample_cube()
        del cube.attributes['_filename']
        with self.assertRaises(ValueError):
            _io.save_cubes([cube])

    def test_save_optimized_map(self):
        """Test save"""
        cube = self._create_sample_cube()
        paths = _io.save_cubes([cube], optimize_access='map')
        loaded_cube = iris.load_cube(paths[0])
        self._compare_cubes(cube, loaded_cube)
        handler = netCDF4.Dataset(paths[0], 'r')
        chunking = handler.variables['sample'].chunking()
        handler.close()
        self.assertListEqual([2, 2, 1], chunking)

    def test_save_optimized_timeseries(self):
        """Test save"""
        cube = self._create_sample_cube()
        paths = _io.save_cubes([cube], optimize_access='timeseries')
        loaded_cube = iris.load_cube(paths[0])
        self._compare_cubes(cube, loaded_cube)
        handler = netCDF4.Dataset(paths[0], 'r')
        chunking = handler.variables['sample'].chunking()
        handler.close()
        self.assertListEqual([1, 1, 2], chunking)

    def test_save_optimized_lat(self):
        """Test save"""
        cube = self._create_sample_cube()
        paths = _io.save_cubes([cube], optimize_access='latitude')
        loaded_cube = iris.load_cube(paths[0])
        self._compare_cubes(cube, loaded_cube)
        handler = netCDF4.Dataset(paths[0], 'r')
        chunking = handler.variables['sample'].chunking()
        handler.close()
        self.assertListEqual([2, 1, 1], chunking)

    def test_save_optimized_lon_time(self):
        """Test save"""
        cube = self._create_sample_cube()
        paths = _io.save_cubes([cube], optimize_access='longitude time')
        loaded_cube = iris.load_cube(paths[0])
        self._compare_cubes(cube, loaded_cube)
        handler = netCDF4.Dataset(paths[0], 'r')
        chunking = handler.variables['sample'].chunking()
        handler.close()
        self.assertListEqual([1, 2, 2], chunking)

    def _compare_cubes(self, cube, loaded_cube):
        self.assertTrue((cube.data == loaded_cube.data).all())
        for coord in cube.coords():
            self.assertTrue((coord.points ==
                             loaded_cube.coord(coord.name()).points).all())
