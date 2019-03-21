# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for IPSL-CM6A-LR."""
import iris

from ..fix import Fix


class fgco2(Fix):
    """Fixes for fgco2."""

    def fix_metadata(self, cubes):
        """Fix metadata.

        Remove unnecessary variables prohibiting cube concatenation.

        Parameters
        ----------
        cubes: iris.cube.CubeList

        Returns
        -------
        iris.cube.CubeList

        """
        cube = cubes.extract_strict(
            'surface_downward_mass_flux_of_carbon_dioxide_expressed_as_carbon')
        cube.coord('latitude').var_name = 'lat'
        cube.coord('longitude').var_name = 'lon'
        return [cube]
