# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for EC-Earth3-HR PRIMAVERA project data"""
from ..fix import Fix


class allvars(Fix):
    """Fixes common to all variables."""

    def fix_metadata(self, cubes):
        """
        Fix cube metadata.

        Parameters
        ----------
        cube: Cube
            Cube to fix

        Returns
        -------
        Cube:
            Fixed cube. It is the same instance that was received
        """
        for cube in cubes:
            latitude = cube.coord('latitude')
            latitude.var_name = 'lat'

            longitude = cube.coord('longitude')
            longitude.var_name = 'lon'
        return cubes
