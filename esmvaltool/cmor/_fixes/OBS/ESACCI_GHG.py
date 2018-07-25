"""Fixes for ESA-CCI GHG"""
import cf_units

from ..fix import Fix


class xco2Stderr(Fix):
    """Fixes for xco2Stderr"""

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
        return cube

    def fix_data(self, cube):
        """
        Fix data

        Fix cube units

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


class xco2Stddev(xco2Stderr):
    """Fixes for xco2Stddev"""

    pass


class xch4Stderr(Fix):
    """Fixes for xch4Stderr"""

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
        cube.units = cf_units.Unit('1.0e-9')
        return cube

    def fix_data(self, cube):
        """
        Fix data

        Fix cube units

        Parameters
        ----------
        cube: iris.cube.Cube
            Cube to fix

        Returns
        -------
        iris.cube.Cube

        """
        metadata = cube.metadata
        cube *= 1.0e9
        cube.metadata = metadata
        return cube


class xch4Stddev(xch4Stderr):
    """Fixes for xch4Stddev"""

    pass
