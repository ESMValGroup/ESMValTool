git"""Fixes for EC-Earth model"""
from ..fix import Fix
import iris

class sic(Fix):
    """Fixes for sic"""

    def fix_data(self, cube):
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
        metadata = cube.metadata
        cube *= 100
        cube.metadata = metadata
        return cube

class tas(Fix):
    """ Fixes for tas"""

    def fix_metadata(self, cube):
        """
        Fix data

        Includes height = 2m attribute

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        # This doesn't seem to work yet!
        height = iris.coords.AuxCoord(2.0, "height", "height", units = "m")
        cube.add_aux_coord(height, data_dims = None)
        return cube

class sftlf(Fix):
    """Fixes for sftlf"""

    def fix_data(self, cube):
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
        metadata = cube.metadata
        cube *= 100
        cube.metadata = metadata
        return cube
