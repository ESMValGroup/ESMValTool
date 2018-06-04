"""Fixes for HadGEM2-ES"""
from ..fix import Fix
import numpy as np

       
class allvars(Fix):
    """Fixes common to all vars"""

    def fix_metadata(self, cube):
        """
        Fixes latitude

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        #print("Fixing HadGEM2-CC Latitude metadata")
        lat = cube.coord('latitude')
        lat.points = np.clip(lat.points,-90.,90.)
        lat.bounds = np.clip(lat.bounds,-90.,90.)
        
        return cube
        
       
class tos(Fix):
    """Fixes common to all vars"""

    def fix_metadata(self, cube):
        """
        Fixes latitude

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        #print("Fixing HadGEM2-CC Latitude metadata")
        lat = cube.coord('latitude')
        lat.points = np.clip(lat.points,-90.,90.)
        lat.bounds = np.clip(lat.bounds,-90.,90.)
        
        return cube
               
