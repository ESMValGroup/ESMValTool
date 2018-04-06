"""Fixes for ACCESS1-0 model"""
from cf_units import Unit

from ..fix import Fix


class rlut(Fix):
    """Fixes for rlut"""

    def fix_metadata(self, cube):
        """
        Fix metadata

        Fixes positive attributte

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        # cube.attributes['positive'] = 'up'
        return cube


class rlutcs(rlut):
    """Fixes for rlutcs"""
    pass
