"""Fixes for CESM1-BGC model"""

from cf_units import Unit

from ..fix import Fix


class co2(Fix):
    """Fixes for co2 variable"""

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
        metadata = cube.metadata
        cube *= 28.966 / 44.0
        cube.metadata = metadata
        return cube


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
        if time.units.name == 'day since 1-01-01 00:00:00.000000 UTC':
            time.units = Unit('days since 1850-01-01 00:00:00',
                              time.units.calendar)
        return cube
