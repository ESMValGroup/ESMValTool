from esmvaltool.interface_scripts.fixes.fix import Fix


# noinspection PyPep8Naming
class fgco2(Fix):

    def fix_data(self, cube):
        return cube * 12.0 / 44.0
