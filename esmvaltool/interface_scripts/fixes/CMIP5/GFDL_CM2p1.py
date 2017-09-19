from esmvaltool.interface_scripts.fixes.fix import Fix


class sftof(Fix):

    def fix_data(self, cube):
        return cube * 100

