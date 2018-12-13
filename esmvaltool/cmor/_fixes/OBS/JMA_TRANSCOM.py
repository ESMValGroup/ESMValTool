"""Fixes for JMA-TRANSCOM."""
from ..fix import Fix


class fgco2(Fix):
    """Class to fix fgco2."""

    def fix_metadata(self, cube):
        """Fix metadata for fgco2.

        Longitude is [-180, 180], needs to be [0, 360].

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        cube = cube.intersection(longitude=(0.0, 360.0))
        return cube
