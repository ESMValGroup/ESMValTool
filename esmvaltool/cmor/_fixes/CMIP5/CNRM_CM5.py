"""Fixes for CNRM-CM5 model"""
import iris
from ..fix import Fix


class msftmyz(Fix):
    """Fixes for msftmyz"""

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
        cube *= 1e6
        cube.metadata = metadata
        return cube


class msftmyzba(msftmyz):
    """Fixes for msftmyzba"""

    pass


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
        depth = cube.coord('ocean depth coordinate')
        depth.standard_name = 'depth'
        depth.long_name = 'depth'
        depth.attributes['positive'] = 'down'
        return cube
