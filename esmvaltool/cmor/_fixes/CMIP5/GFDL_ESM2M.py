# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for GFDL ESM2M"""
from cf_units import Unit

from ..fix import Fix
from ..CMIP5.GFDL_ESM2G import allvars as base_allvars


class allvars(base_allvars):
    """Fixes for all variables"""


class sftof(Fix):
    """Fixes for sftof"""

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
        cube *= 100
        cube.metadata = metadata
        return cube


class co2(Fix):
    """Fixes for co2"""

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
