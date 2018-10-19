"""Fixes for EC-Earth3 CMIP6 project data."""
import iris.coords
import iris.util
from ..fix import Fix


class allvars(Fix):
    """Fixes common to all variables"""

    def fix_metadata(self, cube):
        """
        Fixes cube metadata.

        Parameters
        ----------
        cube: Cube
            Cube to fix

        Returns
        -------
        Cube:
            Fixed cube. It is the same instance that was received
        """
        latitude = cube.coord('latitude')
        latitude.var_name = 'lat'

        longitude = cube.coord('longitude')
        longitude.var_name = 'lon'
        return cube


class siconc(Fix):
    """Fixes common to all variables"""

    def fix_metadata(self, cube):
        """
        Fixes cube metadata.

        Add typesi coordinate

        Parameters
        ----------
        cube: Cube
            Cube to fix

        Returns
        -------
        Cube:
            Fixed cube. It is the same instance that was received
        """
        cube.add_aux_coord(iris.coords.AuxCoord(['sea_ice'],
                                                standard_name='area_type',
                                                var_name='type',
                                                long_name='Sea Ice area type'))
        return cube


class zg(Fix):
    """Fixes for geopotential height variable"""

    def fix_metadata(self, cube):
        """
        Fixes cube metadata.

        Simplify lat lon coordinates to make them 1D

        Parameters
        ----------
        cube: Cube
            Cube to fix

        Returns
        -------
        Cube:
            Fixed cube. It is the same instance that was received
        """
        cube.attributes['realm'] = 'atmos'

        lat_2d = cube.coord('latitude')
        lat_1d = lat_2d.copy(lat_2d.points[:, 0], -lat_2d.bounds[:, 0, 1:3])
        cube.remove_coord('latitude')
        cube.add_aux_coord(lat_1d, 2)

        lon_2d = cube.coord('longitude')
        lon_1d = lon_2d.copy(lon_2d.points[0, :], lon_2d.bounds[0, :, 0:2])
        cube.remove_coord('longitude')
        cube.add_aux_coord(lon_1d, 3)

        iris.util.promote_aux_coord_to_dim_coord(cube, lat_1d)
        iris.util.promote_aux_coord_to_dim_coord(cube, lon_1d)
        return cube


class tas(zg):
    """Fixes for surface temperature variable."""

    pass
