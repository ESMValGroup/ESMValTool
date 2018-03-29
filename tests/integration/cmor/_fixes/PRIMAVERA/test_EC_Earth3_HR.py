"""Test for the fixes for EC-Earth3-HR model from PRIMAVERA project"""
import unittest

from iris.coords import DimCoord
from iris.cube import Cube

from esmvaltool.cmor._fixes.PRIMAVERA.EC_Earth3_HR import allvars


class TestAllVars(unittest.TestCase):
    """Test for the common fixes for all the variables"""
    def setUp(self):
        """Prepare to test"""
        self.cube = Cube([[1, 2], [3, 4]], var_name='var')
        self.cube.add_dim_coord(
            DimCoord([1, 2], standard_name='latitude',
                     var_name='latitude'), 0)
        self.cube.add_dim_coord(
            DimCoord([1, 2], standard_name='longitude',
                     var_name='longitude'), 1)
        self.fix = allvars()

    def test_fix_lat_lon_names(self):
        """
        Test latitude and longitude var names

        They should be lat and lon instead of the original latitude and
        longitude
        """
        cube = self.fix.fix_metadata(self.cube)
        self.assertEqual(cube.coord('latitude').var_name, 'lat')
        self.assertEqual(cube.coord('longitude').var_name, 'lon')
