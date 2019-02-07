"""Fixes for CCSM4 model"""
import numpy as np

from ..fix import Fix


# noinspection PyPep8Naming
class prw(Fix):
    """Fixes for prw"""

    def fix_metadata(self, cubes):
        """
        Fix metadata

        Remove error and number of observations cubes

        Parameters
        ----------
        cube: iris.cube.CubeList

        Returns
        -------
        iris.cube.Cube

        """
        cube = self.get_cube_from_list(cubes)
        return CubeList([cube])