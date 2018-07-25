"""Unit tests for extract_month and season preprocessor function."""

from __future__ import absolute_import, division, print_function

import iris
from iris.cube import Cube
import iris.coords
import iris.coord_categorisation
from cf_units import Unit
import numpy as np

import tests
from esmvaltool.preprocessor._time_area import extract_month, extract_season


def _create_sample_cube():
    cube = Cube(np.arange(1, 25), var_name='co2', units='J')
    cube.add_dim_coord(
        iris.coords.DimCoord(np.arange(15., 720., 30.),
                             standard_name='time',
                             units=Unit('days since 1950-01-01',
                                        calendar='gregorian')),
        0)
    iris.coord_categorisation.add_month_number(cube, 'time')
    return cube


class TestExtractMonth(tests.Test):
    """Tests for extract_month`."""

    def setUp(self):
        """Prepare tests"""
        self.cube = _create_sample_cube()

    def test_get_january(self):
        """Test january extraction"""
        sliced = extract_month(self.cube, 1)
        print(sliced)
        self.assertTrue((np.array([1, 1]) ==
                         sliced.coord('month_number').points).all())


class TestExtractSeason(tests.Test):
    """Tests for extract_season"""

    def setUp(self):
        """Prepare tests"""
        self.cube = _create_sample_cube()

    def test_get_djf(self):
        """Test function for winter"""
        sliced = extract_season(self.cube, 'djf')
        print(sliced)
        self.assertTrue((np.array([1, 2, 12, 1, 2, 12]) ==
                         sliced.coord('month_number').points).all())

    def test_get_djf_caps(self):
        """Test function works when season specified in caps"""
        sliced = extract_season(self.cube, 'DJF')
        print(sliced)
        self.assertTrue((np.array([1, 2, 12, 1, 2, 12]) ==
                         sliced.coord('month_number').points).all())

    def test_get_mam(self):
        """Test function for spring"""
        sliced = extract_season(self.cube, 'mam')
        print(sliced)
        self.assertTrue((np.array([3, 4, 5, 3, 4, 5]) ==
                         sliced.coord('month_number').points).all())

    def test_get_jja(self):
        """Test function for summer"""
        sliced = extract_season(self.cube, 'jja')
        print(sliced)
        self.assertTrue((np.array([6, 7, 8, 6, 7, 8]) ==
                         sliced.coord('month_number').points).all())

    def test_get_son(self):
        """Test function for summer"""
        sliced = extract_season(self.cube, 'son')
        print(sliced)
        self.assertTrue((np.array([9, 10, 11, 9, 10, 11]) ==
                         sliced.coord('month_number').points).all())
