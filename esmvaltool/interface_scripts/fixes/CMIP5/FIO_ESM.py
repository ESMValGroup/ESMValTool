from esmvaltool.interface_scripts.fixes.fix import Fix


class co2(Fix):
    def fix_data(self, cube):
        return cube * 29 / 44 * 1.e6


class ch4(Fix):
    def fix_data(self, cube):
        return cube * 29 / 16 * 1.e9
