"""Fixes for BDBP"""

import cf_units

from ..fix import Fix


class tro3(Fix):
    """Class to fix tro3"""

    def fix_metadata(self, cube):
        """
        Fix metadata for tro3prof

        Fix air_pressure coordinate
        """
        old = cube.coord('air_pressure')
        dims = cube.coord_dims(old)
        cube.remove_coord(old)
        points = old.points * 100
        if old.bounds is None:
            bounds = None
        else:
            bounds = old.bounds * 100
        plev = old.copy(points, bounds)
        plev.var_name = 'plev'
        plev.standard_name = 'air_pressure'
        plev.long_name = 'Pressure '
        plev.units = cf_units.Unit('Pa')
        cube.add_dim_coord(plev, dims)
        return cube

    def fix_data(self, cube):
        """
        Fix data

        Fixes cube units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        metadata = cube.metadata
        cube *= 1.e9
        cube.metadata = metadata
        return cube


class tro3prof(Fix):
    """Class to fix tro3prof"""

    def fix_metadata(self, cube):
        """
        Fix metadata for tro3prof

        Fix air_pressure coordinate
        """
        old = cube.coord('air_pressure')
        dims = cube.coord_dims(old)
        cube.remove_coord(old)
        points = old.points * 100
        if old.bounds is None:
            bounds = None
        else:
            bounds = old.bounds * 100
        plev = old.copy(points, bounds)
        plev.var_name = 'plev'
        plev.standard_name = 'air_pressure'
        plev.long_name = 'Pressure '
        plev.units = cf_units.Unit('Pa')
        cube.add_dim_coord(plev, dims)
        return cube
