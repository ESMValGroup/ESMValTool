from orchestrator.interface_scripts.fixes.fix import Fix
import cf_units


class tro3(Fix):

    def fix_data(self, cube):
        return cube * 1000


class co2(Fix):

    def fix_metadata(self, cube):
        cube.units = cf_units.Unit('1.0e-6')
        return cube


class gpp(Fix):

    def fix_metadata(self, cube):
        # Fixing the metadata, automatic unit conversion should do the trick
        cube.units = cf_units.Unit('g m-2 day-1')
        return cube



