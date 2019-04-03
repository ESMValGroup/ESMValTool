# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for CCSM4 model."""
import numpy as np

from ..fix import Fix


# noinspection PyPep8Naming
class rlut(Fix):
    """Fixes for rlut."""

    def fix_metadata(self, cubes):
        """
        Fix data.

        Fixes discrepancy between declared units and real units

        Parameters
        ----------
        cube: iris.cube.CubeList

        Returns
        -------
        iris.cube.Cube

        """
        cube = self.get_cube_from_list(cubes)
        lat = cube.coord('latitude')
        lat.points = np.round(lat.points, 3)
        lat.bounds = np.round(lat.bounds, 3)
        return cubes


class rlutcs(rlut):
    """Fixes for rlutcs."""


class rsut(rlut):
    """Fixes for rsut."""


class rsutcs(rlut):
    """Fixes for rsutcs."""


class so(Fix):
    """Fixes for so."""

    def fix_metadata(self, cubes):
        """
        Fix data.

        Fixes discrepancy between declared units and real units

        Parameters
        ----------
        cube: iris.cube.CubeList

        Returns
        -------
        iris.cube.Cube

        """
        self.get_cube_from_list(cubes).units = '1e3'

        return cubes
