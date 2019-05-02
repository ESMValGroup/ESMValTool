# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for FGOALS-g2 model"""
from cf_units import Unit

from ..fix import Fix


class allvars(Fix):
    """Fixes common to all vars"""

    def fix_metadata(self, cubes):
        """
        Fix metadata.

        Fixes time units

        Parameters
        ----------
        cube: iris.cube.CubeList

        Returns
        -------
        iris.cube.Cube

        """
        for cube in cubes:
            time = cube.coord('time')
            time.units = Unit(time.units.name, time.units.calendar)
        return cubes
