from orchestrator.interface_scripts.fixes.fix import Fix


class fgco2(Fix):

    def fix_data(self, cube):
        return cube * 12.0 / 44.0
