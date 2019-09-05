# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for FIO ESM model."""
from ..fix import Fix


class co2(Fix):
    """Fixes for co2."""

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
        cube *= 29. / 44. * 1.e6
        cube.metadata = metadata
        return cube


class ch4(Fix):
    """Fixes for ch4."""

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
        cube *= 29. / 16. * 1.e9
        cube.metadata = metadata
        return cube
