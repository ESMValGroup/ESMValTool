"""Fixes for BNU ESM model"""
import iris
from cf_units import Unit

from ..fix import Fix


class fgco2(Fix):
    """Fixes for fgco2"""

    def fix_metadata(self, cube):
        """
        Fix metadata

        Fixes cube units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        cube.units = Unit('kg m-2 s-1')
        return cube

    def fix_data(self, cube):
        """
        Fix data

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
    """Fixes for ch4"""

    def fix_metadata(self, cube):
        """
        Fix metadata

        Fixes cube units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        cube.units = Unit('1e-9')
        return cube

    def fix_data(self, cube):
        """
        Fix metadata

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
    """Fixes for co2"""

    def fix_metadata(self, cube):
        """
        Fix metadata

        Fixes cube units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        cube.units = Unit('1e-6')
        return cube

    def fix_data(self, cube):
        """
        Fix data

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
    """Fixes for spco2"""

    def fix_data(self, cube):
        """
        Fix data

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


class o2(Fix):
    """Fixes for o2"""

    def fix_metadata(self, cube):
        """
        Fix metadata

        Fixes cube units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        depth = cube.coord('generic ocean level')
        depth.standard_name = 'depth'
        depth.long_name = 'depth'
        return cube

    def fix_file(self, filepath, output_dir):
        """
        Apply fixes to the files prior to creating the cube.

        Should be used only to fix errors that prevent loading or can
        not be fixed in the cube (i.e. those related with missing_value
        and _FillValue or missing standard_name).
        Parameters
        ----------
        filepath: basestring
            file to fix.
        output_dir: basestring
            path to the folder to store the fix files, if required.
        Returns
        -------
        basestring
            Path to the corrected file. It can be different from the original
            filepath if a fix has been applied, but if not it should be the
            original filepath.
        """
        new_path = Fix.get_fixed_filepath(output_dir, filepath)
        cube = iris.load_cube(filepath)

        std = 'mole_concentration_of_dissolved_molecular_oxygen_in_sea_water'
        long_name = 'Dissolved Oxygen Concentration'

        cube.long_name = long_name
        cube.standard_name = std

        iris.save(cube, new_path)
        return new_path

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
