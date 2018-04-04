"""Fixes for MIROC ESM model"""
import cf_units
from iris.coords import DimCoord
from iris.exceptions import CoordinateNotFoundError

from ..fix import Fix


class tro3(Fix):
    """Fixes for tro3"""

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
        return cube * 1000


class co2(Fix):
    """Fixes for co2"""

    def fix_metadata(self, cube):
        """
        Fix metadata

        Fixes error in cube units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        cube.units = cf_units.Unit('1.0e-6')
        return cube


class gpp(Fix):
    """Fixes for gpp"""

    def fix_metadata(self, cube):
        """
        Fix metadata

        Fixes error in cube units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        # Fixing the metadata, automatic unit conversion should do the trick
        cube.units = cf_units.Unit('g m-2 day-1')
        return cube


class allvars(Fix):
    """Common fixes to all vars"""

    def fix_metadata(self, cube):
        """
        Fix metadata

        Fixes errors in time units and correct air_pressure coordinate,
        sometimes called AR5PL35

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        try:
            time = cube.coord('time')
            if time.units.calendar:
                calendar = time.units.calendar
            else:
                calendar = 'standard'

            if time.units.origin == 'days since 0000-01-01 00:00:00':
                time.units = cf_units.Unit(
                    'days since 1849-01-01 00:00:00', calendar=calendar)
            elif time.units.origin == 'days since 1-1-1':
                time.units = cf_units.Unit(
                    'days since 1850-01-01 00:00:00', calendar=calendar)
        except CoordinateNotFoundError:
            pass

        try:
            old = cube.coord('AR5PL35')
            dims = cube.coord_dims(old)
            cube.remove_coord(old)

            plev = DimCoord.from_coord(old)
            plev.var_name = plev
            plev.standard_name = 'air_pressure'
            plev.long_name = 'Pressure '
            cube.add_dim_coord(plev, dims)
        except CoordinateNotFoundError:
            pass

        return cube
