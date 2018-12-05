"""Fixes for EC-Earth3-HR PRIMAVERA project data"""
from ..fix import Fix


class allvars(Fix):
    """Fixes common to all variables"""

    def fix_metadata(self, cubes):
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
        latitude = cubes[0].coord('latitude')
        latitude.var_name = 'lat'

        longitude = cubes[0].coord('longitude')
        longitude.var_name = 'lon'
        return cubes
