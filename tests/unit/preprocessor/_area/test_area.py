"""Unit tests for the :func:`esmvaltool.preprocessor._area` module."""

import unittest

import iris
import numpy as np
from cf_units import Unit

import tests
from esmvaltool.preprocessor._area import (
    average_region, extract_named_regions, extract_region)


class Test(tests.Test):
    """Test class for the :func:`esmvaltool.preprocessor._area_pp` module."""

    def setUp(self):
        """Prepare tests."""
        self.coord_sys = iris.coord_systems.GeogCS(
            iris.fileformats.pp.EARTH_RADIUS)
        data = np.ones((5, 5))
        lons = iris.coords.DimCoord(
            [i + .5 for i in range(5)],
            standard_name='longitude',
            bounds=[[i, i + 1.] for i in range(5)],  # [0,1] to [4,5]
            units='degrees_east',
            coord_system=self.coord_sys)
        lats = iris.coords.DimCoord([i + .5 for i in range(5)],
                                    standard_name='latitude',
                                    bounds=[[i, i + 1.] for i in range(5)],
                                    units='degrees_north',
                                    coord_system=self.coord_sys)
        coords_spec = [(lats, 0), (lons, 1)]
        self.grid = iris.cube.Cube(data, dim_coords_and_dims=coords_spec)

        ndata = np.ones((6, 6))
        nlons = iris.coords.DimCoord(
            [i - 2.5 for i in range(6)],
            standard_name='longitude',
            bounds=[[i - 3., i - 2.] for i in range(6)],  # [3,2] to [4,5]
            units='degrees_east',
            coord_system=self.coord_sys)
        nlats = iris.coords.DimCoord(
            [i - 2.5 for i in range(6)],
            standard_name='latitude',
            bounds=[[i - 3., i - 2.] for i in range(6)],
            units='degrees_north',
            coord_system=self.coord_sys)
        coords_spec = [(nlats, 0), (nlons, 1)]
        self.negative_grid = iris.cube.Cube(
            ndata, dim_coords_and_dims=coords_spec)

    def test_average_region_mean(self):
        """Test for area average of a 2D field."""
        result = average_region(self.grid, 'latitude', 'longitude')
        expected = np.array([1.])
        self.assertArrayEqual(result.data, expected)

    def test_average_region_min(self):
        """Test for area average of a 2D field."""
        result = average_region(self.grid, 'latitude', 'longitude',
                                operator='min')
        expected = np.array([1.])
        self.assertArrayEqual(result.data, expected)

    def test_average_region_max(self):
        """Test for area average of a 2D field."""
        result = average_region(self.grid, 'latitude', 'longitude',
                                operator='max')
        expected = np.array([1.])
        self.assertArrayEqual(result.data, expected)

    def test_average_region_median(self):
        """Test for area average of a 2D field."""
        result = average_region(self.grid, 'latitude', 'longitude',
                                operator='median')
        expected = np.array([1.])
        self.assertArrayEqual(result.data, expected)

    def test_average_region_std_dev(self):
        """Test for area average of a 2D field."""
        result = average_region(self.grid, 'latitude', 'longitude',
                                operator='std_dev')
        expected = np.array([0.])
        self.assertArrayEqual(result.data, expected)

    def test_average_region_variance(self):
        """Test for area average of a 2D field."""
        result = average_region(self.grid, 'latitude', 'longitude',
                                operator='variance')
        expected = np.array([0.])
        self.assertArrayEqual(result.data, expected)

    def test_average_region_neg_lon(self):
        """Test for area average of a 2D field."""
        result = average_region(self.negative_grid, 'latitude', 'longitude')
        expected = np.array([1.])
        self.assertArrayEqual(result.data, expected)

    def test_extract_region(self):
        """Test for extracting a region from a 2D field."""
        result = extract_region(self.grid, 1.5, 2.5, 1.5, 2.5)
        # expected outcome
        expected = np.ones((2, 2))
        self.assertArrayEqual(result.data, expected)

    def test_extract_region_neg_lon(self):
        """Test for extracting a region with a negative longitude field."""
        result = extract_region(self.negative_grid, -0.5, 0.5, -0.5, 0.5)
        expected = np.ones((2, 2))
        self.assertArrayEqual(result.data, expected)

    def test_extract_named_region(self):
        """Test for extracting a named region."""
        # tests:
        # Create a cube with regions
        times = np.array([15., 45., 75.])
        bounds = np.array([[0., 30.], [30., 60.], [60., 90.]])
        time = iris.coords.DimCoord(
            times,
            bounds=bounds,
            standard_name='time',
            units=Unit('days since 1950-01-01', calendar='gregorian'))

        regions = ['region1', 'region2', 'region3']
        region = iris.coords.AuxCoord(
            regions,
            standard_name='region',
            units='1',
        )

        data = np.ones((3, 3))
        region_cube = iris.cube.Cube(
            data,
            dim_coords_and_dims=[(time, 0)],
            aux_coords_and_dims=[(region, 1)])

        # test string region
        result1 = extract_named_regions(region_cube, 'region1')
        expected = np.ones((3, ))
        self.assertArrayEqual(result1.data, expected)

        # test list of regions
        result2 = extract_named_regions(region_cube, ['region1', 'region2'])
        expected = np.ones((3, 2))
        self.assertArrayEqual(result2.data, expected)

        # test for expected failures:
        with self.assertRaises(ValueError):
            extract_named_regions(region_cube, 'reg_A')
            extract_named_regions(region_cube, ['region1', 'reg_A'])


if __name__ == '__main__':
    unittest.main()
