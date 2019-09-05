# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for GFDL ESM2M."""
from ..CMIP5.GFDL_ESM2G import allvars as base_allvars
from ..fix import Fix


class allvars(base_allvars):
    """Fixes for all variables"""


class sftof(Fix):
    """Fixes for sftof."""

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
        cube *= 100
        cube.metadata = metadata
        return cube


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


class tos(Fix):
    """Fixes for tos."""

    def fix_metadata(self, cubes):
        """Fix metadata.

        Fixes cube `standard_name`.

        Parameters
        ----------
        cube: iris.cube.CubeList

        Returns
        -------
        iris.cube.Cube

        """
        cube = self.get_cube_from_list(cubes)
        cube.standard_name = 'sea_surface_temperature'
        return [cube]
