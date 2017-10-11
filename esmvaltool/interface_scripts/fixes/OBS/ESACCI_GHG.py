from esmvaltool.interface_scripts.fixes.fix import Fix
import cf_units


class xco2Stderr(Fix):
    def fix_metadata(self, cube):
        cube.units = cf_units.Unit('1.0e-6')
        return cube

    def fix_data(self, cube):
        return cube * 1.0e6


class xco2Stddev(xco2Stderr):
    pass


class xch4Stderr(Fix):
    def fix_metadata(self, cube):
        cube.units = cf_units.Unit('1.0e-9')
        return cube

    def fix_data(self, cube):
        return cube * 1.0e9


class xch4Stddev(xch4Stderr):
    pass
