"""Fixes for ESA-CCI GHG"""
import numpy as np
import cf_units

from iris.cube import CubeList

from ..fix import Fix

class xco2(Fix):
    """Fixes for xco2"""

    def fix_metadata(self, cubes):
        """
        Fix metadata

        Fix cube units

        Parameters
        ----------
        cube: iris.cube.CubeList
            Cube to fix

        Returns
        -------
        iris.cube.CubeList

        """
        print("*** FIX XCO2 METADATA ***")
        cube = cubes[0]
        cube.units = cf_units.Unit('1.0e-6')
        cube = cube.intersection(longitude=(0, 360))
        return CubeList(cube)

    def fix_data(self, cube):
        """
        Fix data

        Fix cube units
        Fix longitudes (change from -180...180 to 0...360)

        Parameters
        ----------
        cube: iris.cube.Cube
            Cube to fix

        Returns
        -------
        iris.cube.Cube

        """
        print("*** FIX XCO2 DATA ***")
        metadata = cube.metadata
        cube *= 1.0e6
        cube.metadata = metadata

        return cube
