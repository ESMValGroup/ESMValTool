# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for inmcm4 model."""
import iris

from ..fix import Fix


class gpp(Fix):
    """Fixes for gpp."""

    def fix_data(self, cube):
        """
        Fix data.

        Fixes discrepancy between declared units and real units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        metadata = cube.metadata
        cube *= -1
        cube.metadata = metadata
        return cube


class lai(Fix):
    """Fixes for lai."""

    def fix_data(self, cube):
        """
        Fix data.

        Fixes discrepancy between declared units and real units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        metadata = cube.metadata
        cube *= 0.01
        cube.metadata = metadata
        return cube


class nbp(Fix):
    """Fixes for nbp."""

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
        cube.standard_name = ('surface_net_downward_mass_flux_of_carbon_'
                              'dioxide_expressed_as_carbon_due_to_all_land_'
                              'processes')
        iris.save(cube, new_path)
        return new_path


class baresoilFrac(Fix):
    """Fixes for baresoilFrac."""

    def fix_metadata(self, cubelist):
        """
        Fix missing scalar dimension.

        Parameters
        ----------
        cubelist: iris CubeList
            List of cubes to fix

        Returns
        -------
        iris.cube.CubeList

        """
        typebare = iris.coords.AuxCoord(
            'bare_ground',
            standard_name='area_type',
            long_name='surface type',
            var_name='type',
            units='1',
            bounds=None)
        for cube in cubelist:
            cube.add_aux_coord(typebare)
        return cubelist
