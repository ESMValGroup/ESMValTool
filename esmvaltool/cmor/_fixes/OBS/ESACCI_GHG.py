"""Fixes for ESA-CCI GHG"""
import cf_units

from ..fix import Fix


class xco2Stderr(Fix):
    """Fixes for xco2Stderr"""

    def fix_metadata(self, cubes):
        """
        Fix metadata

        Fix cube units

        Parameters
        ----------
        cubes: iris.cube.CubeList
            Cube to fix

        Returns
        -------
        iris.cube.CubeList

        """
        cubes[0].units = cf_units.Unit('1.0e-6')
        return cubes

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

    def fix_metadata(self, cubes):
        """
        Fix metadata

        Fix cube units

        Parameters
        ----------
        cubes: iris.cube.CubeList
            Cubes to fix

        Returns
        -------
        iris.cube.CubeList

        """
        cubes[0].units = cf_units.Unit('1.0e-9')
        return cubes

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
