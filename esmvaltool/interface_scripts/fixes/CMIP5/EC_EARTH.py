from esmvaltool.interface_scripts.fixes.fix import Fix


class sic(Fix):

    def fix_data(self, cube):
        return cube * 100


class sftlf(Fix):

    def fix_data(self, cube):
        return cube * 100
