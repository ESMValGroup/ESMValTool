"""Fixes for MPI ESM LR model"""
from ..fix import Fix


class pctisccp(Fix):
    """Fixes for pctisccp"""

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
