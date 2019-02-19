# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for BNU ESM model."""
from cf_units import Unit
from dask import array as da

from ..fix import Fix


class fgco2(Fix):
    """Fixes for fgco2."""

    def fix_metadata(self, cubes):
        """
        Fix metadata.

        Fixes cube units

        Parameters
        ----------
        cube: iris.cube.CubeList

        Returns
        -------
        iris.cube.Cube

        """
        self.get_cube_from_list(cubes).units = Unit('kg m-2 s-1')
        return cubes

    def fix_data(self, cube):
        """
        Fix data.

        Fixes cube units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        metadata = cube.metadata
        cube *= 12.0 / 44.0
        cube.metadata = metadata
        return cube


class ch4(Fix):
    """Fixes for ch4."""

    def fix_metadata(self, cubes):
        """
        Fix metadata.

        Fixes cube units

        Parameters
        ----------
        cubes: iris.cube.CubeList

        Returns
        -------
        iris.cube.CubeList

        """
        self.get_cube_from_list(cubes).units = Unit('1e-9')
        return cubes

    def fix_data(self, cube):
        """
        Fix metadata.

        Fixes cube units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        metadata = cube.metadata
        cube *= 29.0 / 16.0 * 1.e9
        cube.metadata = metadata
        return cube


class co2(Fix):
    """Fixes for co2."""

    def fix_metadata(self, cubes):
        """
        Fix metadata.

        Fixes cube units

        Parameters
        ----------
        cubes: iris.cube.CubeList

        Returns
        -------
        iris.cube.CubeList

        """
        self.get_cube_from_list(cubes).units = Unit('1e-6')
        return cubes

    def fix_data(self, cube):
        """
        Fix data.

        Fixes cube units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        metadata = cube.metadata
        cube *= 29.0 / 44.0 * 1.e6
        cube.metadata = metadata
        return cube


class spco2(Fix):
    """Fixes for spco2."""

    def fix_data(self, cube):
        """
        Fix data.

        Fixes cube units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        metadata = cube.metadata
        cube *= 1.e6
        cube.metadata = metadata
        return cube


class od550aer(Fix):
    """Fixes for od550aer."""

    def fix_data(self, cube):
        """
        Fix data.

        Masks invalid values.

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        data = da.ma.masked_equal(cube.core_data(), 1.e36)
        return cube.copy(data)


# No clear way to apply this fix now that we are working with cubes, not files

# class sftlf(Fix):
#
#     def fix_metadata(self):
#         self.cube = self.cube * 1.e6

#   if (name.eq."sftlf") then
#       files = systemfunc("ls " + INFILE)
#       f=addfile(files(0), "r")
#       tmp=f->lat
#       var&lat = tmp
#       delete(tmp)
#       delete(f)
#       ret = 0
#   end if
#
