from esmvaltool.interface_scripts.fixes.fix import Fix


class co2(Fix):

    def fix_data(self, cube):
        return cube * 1e6

