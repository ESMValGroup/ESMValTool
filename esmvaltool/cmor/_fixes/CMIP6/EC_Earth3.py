"""Fixes for EC-Earth3-HR PRIMAVERA project data"""
import iris.coords
import iris.util
from ..fix import Fix
from netCDF4 import Dataset


class allvars(Fix):
    """Fixes common to all variables"""

    def fix_metadata(self, cube):
        """
        Fixes cube metadata

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
        Fixes cube metadata

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

    def fix_metadata(self, cube):
        cube.attributes['realm'] = 'atmos'

        lat_2D = cube.coord('latitude')
        lat_1D = lat_2D.copy(lat_2D.points[:, 0], -lat_2D.bounds[:, 0, 1:3])
        cube.remove_coord('latitude')
        cube.add_aux_coord(lat_1D, 2)

        lon_2D = cube.coord('longitude')
        lon_1D = lon_2D.copy(lon_2D.points[0, :], lon_2D.bounds[0, :, 0:2])
        cube.remove_coord('longitude')
        cube.add_aux_coord(lon_1D, 3)

        iris.util.promote_aux_coord_to_dim_coord(cube, lat_1D)
        iris.util.promote_aux_coord_to_dim_coord(cube, lon_1D)
        return cube


class tas(zg):
    pass
