"""Fixes for EC-Earth3-HR PRIMAVERA project data"""
from netCDF4 import Dataset
import iris.coords
import iris.util
from iris.cube import CubeList
from ..fix import Fix


class zg(Fix):
    """Fixes common to all variables"""

    def fix_metadata(self, cube):
        """
        Fixes cube metadata

        Parameters
        ----------
        cube: Cube
            Cube to fix

        Returns
        -------
        Cube:
            Fixed cube. It is the same instance that was received
        """
        zg = cube.coord('air_pressure')
        zg.var_name = 'plev'
        slices = CubeList(reversed(
            [lat_slice for lat_slice in cube.slices_over('latitude')]
        ))
        cube = slices.merge_cube()
        return cube
