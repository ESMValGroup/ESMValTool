"""Fixes for GFDL CM2p1 model"""
from ..fix import Fix


class sftof(Fix):
    """Fixes for sftof"""

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
        return cube * 100
