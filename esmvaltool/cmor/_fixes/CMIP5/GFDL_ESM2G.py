"""Fixes for GFDL ESM2G"""
from ..fix import Fix


class co2(Fix):
    """Fixes for co2"""

    def fix_data(self, cube):
        """
        Fix data

        Fixes discrepancy between declared units and real units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        return cube * 1e6
