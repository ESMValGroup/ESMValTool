# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for HadGEM2_CC"""
import iris
import numpy as np

from ..fix import Fix


class allvars(Fix):
    """Fixes common to all vars"""

    def fix_metadata(self, cubes):
        """
        Fixes latitude.

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        for cube in cubes:
            lats = cube.coords('latitude')
            if lats:
                lat = cube.coord('latitude')
                lat.points = np.clip(lat.points, -90., 90.)
                lat.bounds = np.clip(lat.bounds, -90., 90.)

        return cubes


class o2(Fix):
    """Fixes for o2"""

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
