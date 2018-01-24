from esmvaltool.interface_scripts.fixes.fix import Fix


class allvars(Fix):
    def fix_metadata(self, cube):
        latitude = cube.coord('latitude')
        latitude.var_name = 'lat'

        longitude = cube.coord('longitude')
        longitude.var_name = 'lon'
        return cube
