"""Fixes for CCSM4 model."""
# pylint: disable=invalid-name, no-self-use, too-few-public-methods
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


# noinspection PyPep8Naming
class rlutcs(rlut):
    """Fixes for rlutcs."""


# noinspection PyPep8Naming
class rsut(rlut):
    """Fixes for rsut."""


# noinspection PyPep8Naming
class rsutcs(rlut):
    """Fixes for rsutcs."""


# noinspection PyPep8Naming
class rlus(rlut):
    """Fixes for rlus."""


# noinspection PyPep8Naming
class rsus(rlut):
    """Fixes for rsus."""


# noinspection PyPep8Naming
class rsuscs(rlut):
    """Fixes for rsuscs."""


# noinspection PyPep8Naming
class rlds(rlut):
    """Fixes for rlds."""


# noinspection PyPep8Naming
class rldscs(rlut):
    """Fixes for rldscs."""


# noinspection PyPep8Naming
class rsds(rlut):
    """Fixes for rsds."""


# noinspection PyPep8Naming
class rsdscs(rlut):
    """Fixes for rsdscs."""


# noinspection PyPep8Naming
class rsdt(rlut):
    """Fixes for rsdt."""


# noinspection PyPep8Naming
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
