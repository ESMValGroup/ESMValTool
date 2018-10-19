"""Fixes for HadGEM3-GC31-MM CMIP6 project data."""
from ..fix import Fix


class allvars(Fix):
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
        latitude = cube.coord('latitude')
        latitude.var_name = 'lat'

        longitude = cube.coord('longitude')
        longitude.var_name = 'lon'
        return cube
