# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for GFDL ESM2G."""

import iris

from ..fix import Fix


def _get_and_remove(cubes, long_name):
    try:
        cube = cubes.extract_strict(long_name)
        cubes.remove(cube)
    except iris.exceptions.ConstraintMismatchError:
        pass


class allvars(Fix):
    """Common fixes."""

    def fix_metadata(self, cubes):
        """Fix metadata.

        Fixes bad standard names.

        Parameters
        ----------
        cubes: iris.cube.CubeList

        Returns
        -------
        iris.cube.CubeList

        """
        _get_and_remove(cubes, 'Start time for average period')
        _get_and_remove(cubes, 'End time for average period')
        _get_and_remove(cubes, 'Length of average period')
        return cubes


class co2(Fix):
    """Fixes for co2."""

    def fix_data(self, cube):
        """Fix data.

        Fixes discrepancy between declared units and real units.

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


class fgco2(Fix):
    """Fixes for fgco2."""

    def fix_metadata(self, cubes):
        """Fix metadata.

        Remove unnecessary variables prohibiting cube concatenation.

        Parameters
        ----------
        cubes: iris.cube.CubeList

        Returns
        -------
        iris.cube.CubeList

        """
        _get_and_remove(cubes, 'Latitude of tracer (h) points')
        _get_and_remove(cubes, 'Longitude of tracer (h) points')
        return cubes
