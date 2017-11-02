from cf_units import Unit

from esmvaltool.interface_scripts.fixes.fix import Fix


class allvars(Fix):
    def fix_metadata(self, cube):
        time = cube.coord('time')
        time.units = Unit('days since 1850-01-01 00:00:00',
                          time.units.calendar)
        return cube


class sftof(Fix):
    def fix_data(self, cube):
        return cube * 100


class co2(Fix):
    def fix_data(self, cube):
        return cube * 1e6
