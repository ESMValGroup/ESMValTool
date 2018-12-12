"""Fixes for GFDL ESM2M."""
from shutil import copyfile

from netCDF4 import Dataset

from ..fix import Fix


class sftof(Fix):
    """Fixes for sftof."""

    def fix_data(self, cube):
        """Fix data.

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


class co2(Fix):
    """Fixes for co2."""

    def fix_data(self, cube):
        """Fix data.

        Fixes discrepancy between declared units and real units.

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        metadata = cube.metadata
        cube *= 1e6
        cube.metadata = metadata

        return cube


class tos(Fix):
    """Fixes for tos."""

    def fix_file(self, filepath, output_dir):
        """Apply fixes to the files prior to creating the cube.

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
        copyfile(filepath, new_path)
        dataset = Dataset(new_path, 'a')
        dataset.variables['tos'].standard_name = 'sea_surface_temperature'
        dataset.close()
        return new_path
