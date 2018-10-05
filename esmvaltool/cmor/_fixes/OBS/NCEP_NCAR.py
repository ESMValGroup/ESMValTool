"""Fixes for NCEP-NCAR"""

from iris import Constraint

from ..fix import Fix


class zg(Fix):
    """Class to fix clisccp"""

    def fix_raw_cubes(self, cubes):
        zg_constraint = Constraint(cube_func=(lambda c: c.var_name == 'zg'))
        zg_cube = cubes.extract(zg_constraint)[0]
        zg_cube.standard_name = 'geopotential_height'
        return cubes
