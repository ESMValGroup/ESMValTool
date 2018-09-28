"""Fixes for HWSD."""
from ..fix import Fix


class cSoil(Fix):
    """Class to fix cSoil."""

    def fix_data(self, cube):
        """
        Fix data for cSoil.

        Data lives on lon = [0, 360] grid, but lon = 0 is East Asia, not
        Europe. To fix this, the data needs to be rotated.

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        lon_coord = cube.coord('longitude')
        lon_coord.points = lon_coord.points - 180.0
        cube = cube.intersection(longitude=(0.0, 360.0))
        return cube
