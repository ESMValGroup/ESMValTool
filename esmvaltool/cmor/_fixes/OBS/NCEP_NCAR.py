"""Fixes for NCEP-NCAR"""

from iris import Constraint

from ..fix import Fix


class zg(Fix):
    """Class to fix clisccp"""

    def fix_raw_cubes(self, cubes):
        zg_cube = cubes.extract(Constraint(var_name='zg'))[0]
        zg_cube.standard_name = 'geopotential_height'
        return cubes

