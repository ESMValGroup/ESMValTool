"""Fixes for CCSM4 model"""
import numpy as np

from ..fix import Fix


# noinspection PyPep8Naming
class rlut(Fix):
    """Fixes for rlut"""

    def fix_metadata(self, cube):
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
        lat = cube.coord('latitude')
        lat.points = np.round(lat.points, 3)
        lat.bounds = np.round(lat.bounds, 3)
        return cube


class rlutcs(rlut):
    """Fixes for rlut"""

    pass


class so(Fix):
    """Fixes for so"""

    def fix_metadata(self, cube):
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
        cube.units = '1e3'

        return cube
