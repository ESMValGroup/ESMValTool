from orchestrator.interface_scripts.fixes.fix import Fix


class pctisccp(Fix):

    def fix_data(self, cube):
        return cube * 100
