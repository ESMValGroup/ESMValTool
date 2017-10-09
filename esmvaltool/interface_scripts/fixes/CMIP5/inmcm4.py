from esmvaltool.interface_scripts.fixes.fix import Fix


class gpp(Fix):
    def fix_data(self, cube):
        return cube * -1


class lai(Fix):
    def fix_data(self, cube):
        return cube / 100.0
