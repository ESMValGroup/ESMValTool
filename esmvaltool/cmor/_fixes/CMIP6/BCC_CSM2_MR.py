# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for BCC-CSM-MR."""

from ..fix import Fix


class fgco2(Fix):
    """Fixes for fgco2."""

    def fix_metadata(self, cubes):
        """Fix metadata.

        Remove second `latitude` and `longitude` coordinates.

        Parameters
        ----------
        cubes: iris.cube.CubeList

        Returns
        -------
        iris.cube.CubeList

        """
        cube = cubes.extract_strict(
            'surface_downward_mass_flux_of_carbon_dioxide_expressed_as_carbon')
        cube.coord('latitude', dim_coord=True).rename('lat_idx')
        cube.coord('lat_idx').var_name = 'y'
        cube.coord('longitude', dim_coord=True).rename('lon_idx')
        cube.coord('lon_idx').var_name = 'x'
        return [cube]
