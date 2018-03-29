"""Fixes for CanESM2 model"""
from ..fix import Fix


# noinspection PyPep8Naming
class fgco2(Fix):
    """Fixes for fgco2"""

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
        return cube * 12.0 / 44.0
