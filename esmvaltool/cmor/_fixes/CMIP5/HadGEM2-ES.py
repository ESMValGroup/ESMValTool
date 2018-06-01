"""Fixes for HadGEM2-ES"""
from ..fix import Fix
import numpy as np


class lat(Fix):
    """Fixes for lat"""

    def fix_data(self, cube):
        """
        Fix data
        Fixes problem where HadGEM2-ES data has latitude >  90. N
        Parameters
        ----------
        cube: iris.cube.Cube
        Returns
        -------
        iris.cube.Cube
        """
        lat = cube.coord('latitude')
        lat.points = np.clip(lat.points,-90.,90.)
	lat.bounds = np.clip(lat.bounds,-90.,90.)
	
        return cube * 1e6
