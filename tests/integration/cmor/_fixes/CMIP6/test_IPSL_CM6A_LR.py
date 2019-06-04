import unittest

from iris.cube import Cube
from iris.exceptions import CoordinateNotFoundError

from esmvaltool.cmor._fixes.CMIP6.IPSL_CM6A_LR import allvars


class TestAllVars(unittest.TestCase):
    def setUp(self):
        
        self.fix = ch4()

    def test_fix_metadata_ocean_var(self):
        self.cube = Cube([1.0], var_name='ch4')
        self.cube.add_aux_coord(
            AuxCoord([0], var_name='nav_lat', standard_name='latitude')
        )
        self.cube.add_aux_coord(
            AuxCoord([0], var_name='nav_lon', standard_name='longitude')
        )
        self.cell_area = Cube([1.0], standard_name=='cell_area')
        cubes = self.fix.fix_metadata([self.cube, self.cell_area])

        self.assertEqual(len(cubes), 1)
        cube = cubes[0]
        self.assertEqual(cube.coord('latitude').var_name='lat')
        self.assertEqual(cube.coord('longitude').var_name='lon')
        self.cube.coord('cell_area')

    def test_fix_data_other_var(self):
        self.cube = Cube([1.0], var_name='ch4', units='J')
        self.cube.add_aux_coord(
            AuxCoord([0], var_name='nav_lat', standard_name='latitude')
        )
        self.cube.add_aux_coord(
            AuxCoord([0], var_name='nav_lon', standard_name='longitude')
        )
        cubes = self.fix.fix_metadata([self.cube])

        self.assertEqual(len(cubes), 1)
        cube = cubes[0]
        self.assertEqual(cube.coord('latitude').var_name='nav_lat')
        self.assertEqual(cube.coord('longitude').var_name='nav_lon')
        with self.assertRaises(Coordinate)
            self.cube.coord('cell_area')
