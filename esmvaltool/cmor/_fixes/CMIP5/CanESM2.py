# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for CanESM2 model."""
from ..fix import Fix


# noinspection PyPep8Naming
class fgco2(Fix):
    """Fixes for fgco2."""

    def fix_data(self, cube):
        """
        Fix data.

        Fixes discrepancy between declared units and real units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        metadata = cube.metadata
        cube *= 12.0 / 44.0
        cube.metadata = metadata
        return cube
