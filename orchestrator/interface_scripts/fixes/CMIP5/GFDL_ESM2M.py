from orchestrator.interface_scripts.fixes.fix import Fix


class sftof(Fix):

    def fix_data(self, cube):
        return cube * 100

class co2(Fix):

    def fix_data(self, cube):
        return cube * 1e6

