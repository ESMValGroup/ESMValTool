"""
Test mask

Integration tests for the :func:`esmvaltool.preprocessor._mask`
module.

"""

from __future__ import absolute_import, division, print_function

import unittest

import os
import iris
import numpy as np

import tests
from esmvaltool.preprocessor import _mask as mask


class Test(tests.Test):
    """Test class"""

    def test_mask_landsea(self):
        """Test mask_landocean func"""
        fx_data = np.empty((3, 3))
        fx_data[:] = 60.
        new_cube_data = np.empty((3, 3))
        new_cube_data[:] = 200.
        crd_sys = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
        lons = iris.coords.DimCoord(
            [0, 1.5, 3],
            standard_name='longitude',
            bounds=[[0, 1], [1, 2], [2, 3]],
            units='degrees_east',
            coord_system=crd_sys)
        lats = iris.coords.DimCoord(
            [0, 1.5, 3],
            standard_name='latitude',
            bounds=[[0, 1], [1, 2], [2, 3]],
            units='degrees_north',
            coord_system=crd_sys)
        coords_spec = [(lats, 0), (lons, 1)]
        fx_mask = iris.cube.Cube(fx_data, dim_coords_and_dims=coords_spec)
        iris.save(fx_mask, 'sftlf_test.nc')
        new_cube_land = iris.cube.Cube(new_cube_data,
                                       dim_coords_and_dims=coords_spec)
        new_cube_sea = iris.cube.Cube(new_cube_data,
                                      dim_coords_and_dims=coords_spec)
        # mask with fx files
        result_land = mask.mask_landsea(new_cube_land,
                                        ['sftlf_test.nc'], 'land')
        result_sea = mask.mask_landsea(new_cube_sea, ['sftlf_test.nc'], 'sea')
        expected = np.ma.empty((3, 3))
        expected.data[:] = 200.
        expected.mask = np.ones((3, 3), bool)
        # set fillvalues so we are sure they are equal
        np.ma.set_fill_value(result_land.data, 1e+20)
        np.ma.set_fill_value(result_sea.data, 1e+20)
        np.ma.set_fill_value(expected, 1e+20)
        self.assertArrayEqual(result_land.data.mask, expected.mask)
        expected.mask = np.zeros((3, 3), bool)
        self.assertArrayEqual(result_sea.data, expected)
        # remove the fx.nc temporary file
        os.remove('sftlf_test.nc')

        # mask with shp files
        new_cube_land = iris.cube.Cube(new_cube_data,
                                       dim_coords_and_dims=coords_spec)
        new_cube_sea = iris.cube.Cube(new_cube_data,
                                      dim_coords_and_dims=coords_spec)
        # bear in mind all points are in the ocean
        result_land = mask.mask_landsea(new_cube_land, None, 'land')
        np.ma.set_fill_value(result_land.data, 1e+20)
        expected.mask = np.zeros((3, 3), bool)
        self.assertArrayEqual(result_land.data, expected)


if __name__ == '__main__':
    unittest.main()
