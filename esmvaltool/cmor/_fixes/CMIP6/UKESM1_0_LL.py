# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for UKESM1-0-LL."""
# import numpy as np
# import iris

from ..fix import Fix


class allvars(Fix):
    """Fixes for all vars."""

    def fix_metadata(self, cubelist):
        """
        Fix non-standard dimension names.

        Parameters
        ----------
        cubelist: iris CubeList
            List of cubes to fix

        Returns
        -------
        iris.cube.CubeList
        """
        for cube in cubelist:
            lats = cube.coord('latitude')
            lons = cube.coord('longitude')
            lats.var_name = 'lat'
            lons.var_name = 'lon'
        return cubelist
