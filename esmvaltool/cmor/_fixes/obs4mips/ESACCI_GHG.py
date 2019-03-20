"""Fixes for ESA-CCI GHG"""
import numpy as np
import cf_units

from ..fix import Fix


class xco2(Fix):
    """Fixes for xco2"""

    def fix_metadata(self, cube):
        """
        Fix metadata

        Fix cube units

        Parameters
        ----------
        cube: iris.cube.Cube
            Cube to fix

        Returns
        -------
        iris.cube.Cube

        """
        cube.units = cf_units.Unit('1.0e-6')
        rlon = cube.coord('longitude').points
        if (np.amin(rlon) < 0.):
            n = int(len(rlon) / 2)
            rlon = np.where(rlon < 0., rlon + 360., rlon)
            cube.coord('longitude').points = np.roll(rlon, n)

        print("*** FIX XCO2 METADATA ***")
        return cube

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

        rlon = cube.coord('longitude').points
        n = int(len(rlon) / 2)
        # longitude is assumed to be the last dimension in the numpy array
        ax = len(cube.shape) - 1
        cube.data = np.roll(cube.data, n, axis=ax)

        cube.metadata = metadata

        return cube
