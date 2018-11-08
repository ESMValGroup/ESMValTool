"""Fixes for ESA CCI ozone"""
import cf_units

from ..fix import Fix


class tro3prof(Fix):
    """Fixes for tro3prof"""

    def fix_metadata(self, cube):
        """
        Fix metadata

        Fixes air_pressure coordinate

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

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
