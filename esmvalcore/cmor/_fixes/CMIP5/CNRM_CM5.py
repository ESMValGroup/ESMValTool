# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for CNRM-CM5 model."""
from ..fix import Fix


class msftmyz(Fix):
    """Fixes for msftmyz."""

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
        cube *= 1e6
        cube.metadata = metadata
        return cube


class msftmyzba(msftmyz):
    """Fixes for msftmyzba."""
