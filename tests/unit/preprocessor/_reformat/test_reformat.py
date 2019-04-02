"""Unit tests for the :func:`esmvaltool.preprocessor._reformat` module."""
import os
import unittest

import numpy as np
from esmvaltool.cmor.check import CMORCheckError
from esmvaltool.preprocessor._reformat import cmor_fix_fx

import iris
import tests


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
        var = {}
        var['project'] = 'CMIP5'
        var['dataset'] = 'GFDL-CM3'
        var['cmor_table'] = 'CMIP5'
        var['filename'] = os.path.join(
            'test-reports', 'fx_files', 'sftlf',
            'sftlf_fx_GFDL-CM3_historical_r0i0p0.nc')
        if not os.path.exists(os.path.dirname(var['filename'])):
            os.makedirs(os.path.dirname(var['filename']))
        iris.save(self.good_cube,
                  os.path.join(os.path.dirname(var['filename']),
                               'sftlf_fx_GFDL-CM3_historical_r0i0p0.nc'))
        fx_files = {
            'sftlf': os.path.join(os.path.dirname(var['filename']),
                                  'sftlf_fx_GFDL-CM3_historical_r0i0p0.nc')}
        cmor_fix_fx(fx_files, var)
        saved_cube = iris.load_cube(var['filename'])
        self.assertArrayEqual(self.good_cube.data, saved_cube.data)
        self.assertArrayEqual(self.good_cube.coord('longitude').points,
                              saved_cube.coord('longitude').points)
        self.assertArrayEqual(self.good_cube.coord('latitude').points,
                              saved_cube.coord('latitude').points)

    def test_cmor_fix_fx_withfix(self):
        """Test the checks and fixes on the two fx test files."""
        var = {}
        var['project'] = 'CMIP5'
        var['dataset'] = 'HadGEM2-ES'
        var['cmor_table'] = 'CMIP5'
        var['filename'] = os.path.join(
            'test-reports', 'fx_files', 'mrsofc',
            'mrsofc_fx_HadGEM2-ES_historical_r0i0p0.nc')
        if not os.path.exists(os.path.dirname(var['filename'])):
            os.makedirs(os.path.dirname(var['filename']))

        # break the units a bit
        self.to_fix_cube.units = "kg"

        iris.save(self.to_fix_cube,
                  os.path.join(os.path.dirname(var['filename']),
                               'mrsofc_fx_HadGEM2-ES_historical_r0i0p0.nc'))
        fx_files = {
            'mrsofc': os.path.join(
                os.path.dirname(var['filename']),
                'mrsofc_fx_HadGEM2-ES_historical_r0i0p0.nc')}
        with self.assertRaises(CMORCheckError):
            cmor_fix_fx(fx_files, var)

        # re-fix the units and re-test
        self.to_fix_cube.units = "kg m-2"
        iris.save(self.to_fix_cube,
                  os.path.join(os.path.dirname(var['filename']),
                               'mrsofc_fx_HadGEM2-ES_historical_r0i0p0.nc'))

        saved_cube = iris.load_cube(var['filename'])
        self.assertArrayEqual(self.to_fix_cube.data, saved_cube.data)
        self.assertArrayEqual(self.to_fix_cube.coord('longitude').points,
                              saved_cube.coord('longitude').points)
        # TODO cmor checks in CMIP5 are loose: this case should return
        # a CMORCheckError as it does for CMIP6
        self.assertArrayEqual(np.array([86., 87., 88., 89., 92.000000073]),
                              saved_cube.coord('latitude').points)

    def test_cmor_fix_fx_fail(self):
        """Test the checks and fixes on the two fx test files."""
        var = {}
        var['project'] = 'CMIP6'
        var['dataset'] = 'HadGEM2-ES'
        var['cmor_table'] = 'CMIP6'
        var['filename'] = os.path.join(
            'test-reports', 'fx_files', 'mrsofc',
            'mrsofc_fx_HadGEM2-ES_historical_r0i0p0.nc')
        if not os.path.exists(os.path.dirname(var['filename'])):
            os.makedirs(os.path.dirname(var['filename']))

        iris.save(self.to_fix_cube,
                  os.path.join(os.path.dirname(var['filename']),
                               'mrsofc_fx_HadGEM2-ES_historical_r0i0p0.nc'))
        fx_files = {
            'mrsofc': os.path.join(
                os.path.dirname(var['filename']),
                'mrsofc_fx_HadGEM2-ES_historical_r0i0p0.nc')}

        with self.assertRaises(CMORCheckError):
            cmor_fix_fx(fx_files, var)

if __name__ == '__main__':
    unittest.main()
