from esmvaltool.interface_scripts.fixes.fix import Fix
import cf_units


class tro3prof(Fix):
    def fix_metadata(self, cube):
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

