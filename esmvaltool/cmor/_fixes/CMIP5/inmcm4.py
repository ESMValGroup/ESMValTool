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

    def fix_metadata(self, cubes):
        """Fix metadata.

        Fixes cube `standard_name`.

        Parameters
        ----------
        cube: iris.cube.CubeList

        Returns
        -------
        iris.cube.Cube

        """
        cube = self.get_cube_from_list(cubes)
        cube.standard_name = (
            'surface_net_downward_mass_flux_of_carbon_dioxide_expressed_as_'
            'carbon_due_to_all_land_processes')
        return [cube]


class fgco2(Fix):
    """Fixes for fgco2."""

    def fix_metadata(self, cubes):
        """Fix metadata.

        Fixes cube `standard_name`.

        Parameters
        ----------
        cube: iris.cube.CubeList

        Returns
        -------
        iris.cube.Cube

        """
        cube = self.get_cube_from_list(cubes)
        cube.standard_name = (
            'surface_downward_mass_flux_of_carbon_dioxide_expressed_as_carbon')
        return [cube]


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
