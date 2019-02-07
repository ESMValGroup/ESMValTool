# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for SSMI model."""
from ..fix import Fix


class prw(Fix):
    """Fixes for prw."""

    def fix_metadata(self, cubes):
        for cube in cubes:
            latitude = cube.coord('latitude')
            latitude.var_name = 'lat'

            longitude = cube.coord('longitude')
            longitude.var_name = 'lon'
        return cubes
