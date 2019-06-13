# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for CESM1-BGC model."""

from cf_units import Unit

from ..fix import Fix


class co2(Fix):
    """Fixes for co2 variable."""

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
        cube *= 28.966 / 44.0
        cube.metadata = metadata
        return cube
