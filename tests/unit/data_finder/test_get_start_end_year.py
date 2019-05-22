"""Unit tests for :func:`esmvaltool._data_finder.regrid._stock_cube`"""

import unittest
import os
import tempfile

from esmvaltool._data_finder import get_start_end_year


class TestGetStartEndYear(unittest.TestCase):
    """Tests for get_start_end_year function"""

    def setUp(self):
        descriptor, self.temp_file = tempfile.mkstemp('.nc')
        os.close(descriptor)

    def tearDown(self):
        if os.path.isfile(self.temp_file):
            os.remove(self.temp_file)

    def test_years_at_the_end(self):
        """Test parse files with two years at the end"""
        start, end = get_start_end_year('var_whatever_1980-1981')
        self.assertEqual(1980, start)
        self.assertEqual(1981, end)

    def test_one_year_at_the_end(self):
        """Test parse files with one year at the end"""
        start, end = get_start_end_year('var_whatever_1980.nc')
        self.assertEqual(1980, start)
        self.assertEqual(1980, end)

    def test_full_dates_at_the_end(self):
        """Test parse files with two dates at the end"""
        start, end = get_start_end_year('var_whatever_19800101-19811231.nc')
        self.assertEqual(1980, start)
        self.assertEqual(1981, end)

    def test_one_fulldate_at_the_end(self):
        """Test parse files with one date at the end"""
        start, end = get_start_end_year('var_whatever_19800101.nc')
        self.assertEqual(1980, start)
        self.assertEqual(1980, end)

    def test_years_at_the_start(self):
        """Test parse files with two years at the start"""
        start, end = get_start_end_year('1980-1981_var_whatever.nc')
        self.assertEqual(1980, start)
        self.assertEqual(1981, end)

    def test_one_year_at_the_start(self):
        """Test parse files with one year at the start"""
        start, end = get_start_end_year('1980_var_whatever.nc')
        self.assertEqual(1980, start)
        self.assertEqual(1980, end)

    def test_full_dates_at_the_start(self):
        """Test parse files with two dates at the start"""
        start, end = get_start_end_year('19800101-19811231_var_whatever.nc')
        self.assertEqual(1980, start)
        self.assertEqual(1981, end)

    def test_one_fulldate_at_the_start(self):
        """Test parse files with one date at the start"""
        start, end = get_start_end_year('19800101_var_whatever.nc')
        self.assertEqual(1980, start)
        self.assertEqual(1980, end)

    def test_start_and_date_in_name(self):
        """Test parse one date at the start and one in experiment's name"""
        start, end = get_start_end_year(
            '19800101_var_control-1950_whatever.nc')
        self.assertEqual(1980, start)
        self.assertEqual(1980, end)

    def test_end_and_date_in_name(self):
        """Test parse one date at the end and one in experiment's name"""
        start, end = get_start_end_year(
            'var_control-1950_whatever_19800101.nc')
        self.assertEqual(1980, start)
        self.assertEqual(1980, end)

    def test_read_file_if_no_date_present(self):
        """Test raises if no date is present"""
        import iris
        from iris.cube import Cube
        from iris.coords import DimCoord
        cube = Cube([0,0], var_name='var')
        time = DimCoord([0, 366], 'time', units='days since 1-1-1990')
        cube.add_dim_coord(time, 0)
        iris.save(cube, self.temp_file)
        start, end = get_start_end_year(self.temp_file)
        self.assertEqual(1990, start)
        self.assertEqual(1991, end)

    def test_fails_if_no_date_present(self):
        """Test raises if no date is present"""
        with self.assertRaises(ValueError):
            get_start_end_year('var_whatever')
