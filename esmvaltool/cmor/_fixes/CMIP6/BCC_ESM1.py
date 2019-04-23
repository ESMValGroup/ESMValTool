"""Fixes for BCC-ESM1."""
# import numpy as np
# import iris
from ..fix import Fix


class allvars(Fix):
    """Common fixes to all vars"""

    def fix_metadata(self, cubes):
        """
        Fix metadata.

        Fixes error in time coordinate, sometimes contains trailing zeros

        Parameters
        ----------
        cube: iris.cube.CubeList

        Returns
        -------
        iris.cube.CubeList

        """
        # Remove 1D lat and lon as they are incorrect, this is an irregular
        # 2D grid and lat/lon need to be 2d.
        for cube in cubes:
            coords_to_remove = []
            for coord in cube.coords():
                if coord.var_name in ['time', ]:
                    continue
                if coord.var_name in ['lat', 'lon']:
                    coords_to_remove.append(coord)

            for coord in coords_to_remove:
                cube.remove_coord(coord)

            latitudes = cube.coord('latitude')
            longitudes = cube.coord('longitude')
            latitudes.var_name = 'lat'
            longitudes.var_name = 'lon'

        return cubes
