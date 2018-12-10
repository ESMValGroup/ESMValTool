"""Fixes for MRI-CGCM3 model"""
import numpy as np

from ..fix import Fix


class msftmyz(Fix):
    """Fixes for msftmyz."""
    def fix_data(self, cube):
        """
        Fix data

        Fixes mask

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        cube.data = np.ma.array(cube.data)
        cube.data = np.ma.masked_where(cube.data.mask + (cube.data == 0.),
                                       cube.data)

        return cube
