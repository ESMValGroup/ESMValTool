from orchestrator.interface_scripts.fixes.fix import Fix


class msftmyz(Fix):

    def fix_data(self, cube):
        return cube * 1e6


class msftmyzba(Fix):

    def fix_data(self, cube):
        return cube * 1e6
