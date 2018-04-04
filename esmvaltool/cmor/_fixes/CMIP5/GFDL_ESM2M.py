"""Fixes for GFDL ESM2M"""
from cf_units import Unit

from ..fix import Fix


class allvars(Fix):
    """Fixes common to all vars"""

    def fix_metadata(self, cube):
        """
        Fix metadata

        Fixes errors in time units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        time = cube.coord('time')
        if time.units.origin in ("days since 0001-01-01 00:00:00",
                                 "days since 1-01-01 00:00:00"):
            time.units = Unit('days since 1850-01-01 00:00:00',
                              time.units.calendar)
        return cube


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
