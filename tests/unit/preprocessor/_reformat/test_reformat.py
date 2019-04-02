"""Unit tests for the :func:`esmvaltool.preprocessor._reformat` module."""
import os
import unittest

import numpy as np
from esmvaltool.cmor.check import CMORCheckError
from esmvaltool.preprocessor._reformat import cmor_fix_fx

import iris
import tests


def set_paths(fx_var):
    """Set paths and make dirs."""
    in_path = os.path.join('test-reports', 'fx_files', fx_var)
    if not os.path.exists(in_path):
        os.makedirs(in_path)
    out_path = os.path.join(
        'test-reports', 'fx_files_fixed', fx_var)
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    return in_path, out_path


def make_dicts(cube, project, dataset, fx_var, file_name):
    """Make the variable and fx_files dicts."""
    var = {}
    var['project'] = project
    var['dataset'] = dataset
    var['cmor_table'] = project

    # make paths and dirs
    in_path, out_path = set_paths(fx_var)
    iris.save(cube, os.path.join(in_path, file_name))
    fx_files = {fx_var: os.path.join(in_path, file_name)}
    var['filename'] = os.path.join(out_path, file_name)
    return var, fx_files


class Test(tests.Test):
    """Test class for the :func:`esmvaltool.preprocessor._reformat` module."""

    def setUp(self):
        """Prepare tests."""
        self.coord_sys = iris.coord_systems.GeogCS(
            iris.fileformats.pp.EARTH_RADIUS)
        data = np.ones((5, 5))
        lons = iris.coords.DimCoord(
            [i + .5 for i in range(5)],
            standard_name='longitude', var_name='lon',
            bounds=[[i, i + 1.] for i in range(5)],
            units='degrees_east',
            coord_system=self.coord_sys)
        lats = iris.coords.DimCoord([i + .5 for i in range(5)],
                                    standard_name='latitude', var_name='lat',
                                    bounds=[[i, i + 1.] for i in range(5)],
                                    units='degrees_north',
                                    coord_system=self.coord_sys)
        coords_spec = [(lats, 0), (lons, 1)]
        self.good_cube = iris.cube.Cube(data, dim_coords_and_dims=coords_spec)
        self.good_cube.standard_name = 'land_area_fraction'

        bound_points = [86., 87., 88., 89., 90.]
        fixed_bounds = [[i - 0.3, i + 0.3] for i in bound_points]
        nlats = iris.coords.DimCoord(
            [86., 87., 88., 89., 92.000000073],
            standard_name='latitude', var_name='lat',
            bounds=fixed_bounds,
            units='degrees_north',
            coord_system=self.coord_sys)
        coords_spec = [(nlats, 0), (lons, 1)]
        self.to_fix_cube = iris.cube.Cube(
            data, dim_coords_and_dims=coords_spec)
        self.to_fix_cube.standard_name = \
            'soil_moisture_content_at_field_capacity'
        self.to_fix_cube.units = "kg m-2"

    def test_cmor_fix_fx_nofix(self):
        """Test the checks and fixes on the two fx test files."""
        var, fx_files = make_dicts(self.good_cube, 'CMIP5',
                                   'GFDL-CM3', 'sftlf',
                                   'sftlf_fx_GFDL-CM3_historical_r0i0p0.nc')
        fixed_cubes = cmor_fix_fx(fx_files, var)
        fixed_cube = iris.load_cube(fixed_cubes['sftlf'])
        self.assertArrayEqual(self.good_cube.data, fixed_cube.data)
        self.assertArrayEqual(self.good_cube.coord('longitude').points,
                              fixed_cube.coord('longitude').points)
        self.assertArrayEqual(self.good_cube.coord('latitude').points,
                              fixed_cube.coord('latitude').points)

    def test_cmor_fix_fx_withfix(self):
        """Test the checks and fixes on the two fx test files."""
        # break the units a bit and compute
        self.to_fix_cube.units = "kg"
        var, fx_files = make_dicts(
            self.to_fix_cube, 'CMIP5',
            'HadGEM2-ES', 'mrsofc',
            'mrsofc_fx_HadGEM2-ES_historical_r0i0p0.nc')
        with self.assertRaises(CMORCheckError):
            cmor_fix_fx(fx_files, var)

        # re-fix the units and re-test
        self.to_fix_cube.units = "kg m-2"
        var, fx_files = make_dicts(
            self.to_fix_cube, 'CMIP5',
            'HadGEM2-ES', 'mrsofc',
            'mrsofc_fx_HadGEM2-ES_historical_r0i0p0.nc')
        fixed_cubes = cmor_fix_fx(fx_files, var)
        fixed_cube = iris.load_cube(fixed_cubes['mrsofc'])
        self.assertArrayEqual(self.to_fix_cube.data, fixed_cube.data)
        self.assertArrayEqual(self.to_fix_cube.coord('longitude').points,
                              fixed_cube.coord('longitude').points)
        self.assertArrayEqual(np.array([86., 87., 88., 89., 90.]),
                              fixed_cube.coord('latitude').points)

    def test_cmor_fix_fx_fail(self):
        """Test the checks and fixes on the two fx test files."""
        var, fx_files = make_dicts(
            self.good_cube, 'CMIP6',
            'HadGEM2-ES', 'mrsofc',
            'mrsofc_fx_HadGEM2-ES_historical_r0i0p0.nc')
        with self.assertRaises(CMORCheckError):
            cmor_fix_fx(fx_files, var)

if __name__ == '__main__':
    unittest.main()
