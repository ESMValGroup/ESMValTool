"""Fixes for ghgcci (ESA-CCI GHG)"""
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
        # the obs4MIPs file contains different variables
        # (xco2, xco2_nobs, xco2_stddev, xco2_stderr), so
        # we have to manually pick the "xco2" cube.
        cube = cubes.extract('dry_atmosphere_mole_fraction_of_carbon_dioxide')[0]
        cube.units = cf_units.Unit('1.0e-6')
        cube = cube.intersection(longitude=(0, 360))
        return CubeList([cube])

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
        metadata = cube.metadata
        cube *= 1.0e6
        cube.metadata = metadata

        return cube
