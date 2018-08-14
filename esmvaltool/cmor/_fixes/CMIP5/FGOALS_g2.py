"""Fixes for FGOALS-g2 model"""
from cf_units import Unit

from ..fix import Fix


class allvars(Fix):
    """Fixes common to all vars"""

    def fix_metadata(self, cube):
        """
        Fix metadata

        Fixes time units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        time = cube.coord('time')
        time.units = Unit(time.units.name, time.units.calendar)
        return cube
