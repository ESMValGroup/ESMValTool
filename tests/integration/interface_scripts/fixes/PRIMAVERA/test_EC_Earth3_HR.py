import unittest

from iris.coords import DimCoord
from iris.cube import Cube

from esmvaltool.interface_scripts.fixes.PRIMAVERA.EC_Earth3_HR import allvars


class TestAllVars(unittest.TestCase):
    def setUp(self):
        self.cube = Cube([[1, 2], [3, 4]], var_name='var')
        self.cube.add_dim_coord(
            DimCoord([1, 2], standard_name='latitude',
                     var_name='latitude'), 0)
        self.cube.add_dim_coord(
            DimCoord([1, 2], standard_name='longitude',
                     var_name='longitude'), 1)
        self.fix = allvars()

    def test_fix_lat_lon_names(self):
        cube = self.fix.fix_metadata(self.cube)
        self.assertEqual(cube.coord('latitude').var_name, 'lat')
        self.assertEqual(cube.coord('longitude').var_name, 'lon')
