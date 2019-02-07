# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for GFDL ESM2G"""
import iris
from iris.coords import AuxCoord
from ..fix import Fix


class allvars(Fix):
    """Common fixes."""

    def fix_metadata(self, cubes):
        """
        Fix metadata.

        Fixes bad standard names

        Parameters
        ----------
        cubes: iris.cube.CubeList

        Returns
        -------
        iris.cube.Cube

        """
        self._get_and_remove(cubes, 'Start time for average period')
        self._get_and_remove(cubes, 'End time for average period')
        self._get_and_remove(cubes, 'Length of average period')
        return cubes

    def _get_and_remove(self, cubes, long_name):
        try:
            cube = cubes.extract_strict(long_name)
            cubes.remove(cube)
        except iris.exceptions.ConstraintMismatchError:
            pass


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
        cube *= 1e6
        cube.metadata = metadata
        return cube
