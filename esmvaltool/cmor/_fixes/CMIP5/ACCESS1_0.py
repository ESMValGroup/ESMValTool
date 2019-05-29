"""Fixes for the ACCESS1-0 model."""
# pylint: disable=invalid-name
from cf_units import Unit
import iris
from ..fix import Fix


# noinspection PyPep8
class allvars(Fix):
    """Common fixes to all vars."""

    def fix_metadata(self, cubes):
        """
        Fix metadata.

        Fixes wrong calendar 'gregorian' instead of 'proleptic_gregorian'

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        for cube in cubes:
            try:
                time = cube.coord('time')
            except iris.exceptions.CoordinateNotFoundError:
                continue
            else:
                time.units = Unit(time.units.name, 'gregorian')
        return cubes
