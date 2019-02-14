"""Fixes for JMA-TRANSCOM."""
from ..fix import Fix


class fgco2(Fix):
    """Class to fix fgco2."""

    def fix_metadata(self, cubes):
        """Fix metadata for fgco2.

        Longitude is [-180, 180], needs to be [0, 360].

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        cube = self.get_cube_from_list(cubes)
        cube = cube.intersection(longitude=(0.0, 360.0))
        return [cube]


class nbp(Fix):
    """Class to fix nbp."""

    def fix_data(self, cube):
        """
        Fix data for nbp.

        Data lives on lon = [0, 360] grid, but lon = 0 is East Asia, not
        Europe. To fix this, the data needs to be rotated.

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        lon_coord = cube.coord('longitude')
        lon_coord.points = lon_coord.points - 180.0
        cube = cube.intersection(longitude=(0.0, 360.0))
        return cube
