"""Integration tests for :func:`esmvaltool.preprocessor.save`"""

import os
import tempfile
import unittest

import iris
import netCDF4
import numpy as np
from iris.coords import DimCoord
from iris.cube import Cube

from esmvaltool.preprocessor import save


class TestSave(unittest.TestCase):
    """Tests for :func:`esmvaltool.preprocessor.save`"""

    def setUp(self):
        self.temp_files = []

    def tearDown(self):
        for temp_file in self.temp_files:
            if os.path.isfile(temp_file):
                os.remove(temp_file)

    def _create_sample_cube(self):
        lat = DimCoord(
            np.asarray([1, 2], np.single),
            standard_name='latitude',
            units='degrees_north')
        lon = DimCoord(
            np.asarray([1, 2], np.single),
            standard_name='longitude',
            units='degrees_east')
        time = DimCoord(
            np.asarray([1, 2], np.single),
            standard_name='time',
            units='days since 2000-1-1')

        cube = Cube(
            np.random.random_sample([2, 2, 2]),
            var_name='sample',
            units='1',
            dim_coords_and_dims=((lat, 0), (lon, 1), (time, 2)))

        descriptor, filename = tempfile.mkstemp('.nc')
        os.close(descriptor)
        self.temp_files.append(filename)
        return cube, filename

    def test_save(self):
        """Test save"""
        cube, filename = self._create_sample_cube()
        path = save([cube], filename)
        loaded_cube = iris.load_cube(path)
        self._compare_cubes(cube, loaded_cube)

    def test_save_zlib(self):
        """Test save"""
        cube, filename = self._create_sample_cube()
        path = save([cube], filename, compress=True)
        loaded_cube = iris.load_cube(path)
        self._compare_cubes(cube, loaded_cube)
        handler = netCDF4.Dataset(path, 'r')
        sample_filters = handler.variables['sample'].filters()
        self.assertTrue(sample_filters['zlib'])
        self.assertTrue(sample_filters['shuffle'])
        self.assertEqual(sample_filters['complevel'], 4)
        handler.close()

    def test_fail_without_filename(self):
        """Test save fails if filename is not provided."""
        cube, _ = self._create_sample_cube()
        with self.assertRaises(TypeError):
            save([cube])

    def test_save_optimized_map(self):
        """Test save"""
        cube, filename = self._create_sample_cube()
        path = save([cube], filename, optimize_access='map')
        loaded_cube = iris.load_cube(path)
        self._compare_cubes(cube, loaded_cube)
        self._check_chunks(path, [2, 2, 1])

    def test_save_optimized_timeseries(self):
        """Test save"""
        cube, filename = self._create_sample_cube()
        path = save([cube], filename, optimize_access='timeseries')
        loaded_cube = iris.load_cube(path)
        self._compare_cubes(cube, loaded_cube)
        self._check_chunks(path, [1, 1, 2])

    def test_save_optimized_lat(self):
        """Test save"""
        cube, filename = self._create_sample_cube()
        path = save([cube], filename, optimize_access='latitude')
        loaded_cube = iris.load_cube(path)
        self._compare_cubes(cube, loaded_cube)
        expected_chunks = [2, 1, 1]
        self._check_chunks(path, expected_chunks)

    def _check_chunks(self, path, expected_chunks):
        handler = netCDF4.Dataset(path, 'r')
        chunking = handler.variables['sample'].chunking()
        handler.close()
        self.assertListEqual(expected_chunks, chunking)

    def test_save_optimized_lon_time(self):
        """Test save"""
        cube, filename = self._create_sample_cube()
        path = save([cube], filename, optimize_access='longitude time')
        loaded_cube = iris.load_cube(path)
        self._compare_cubes(cube, loaded_cube)
        self._check_chunks(path, [1, 2, 2])

    def _compare_cubes(self, cube, loaded_cube):
        self.assertTrue((cube.data == loaded_cube.data).all())
        for coord in cube.coords():
            self.assertTrue(
                (coord.points == loaded_cube.coord(coord.name()).points).all())
