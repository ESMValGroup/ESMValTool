"""Derivation of variable `nbp_grid`."""

import iris
from iris import Constraint

from ._baseclass import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `nbp_grid`."""

    # Required variables
    required = [{
        'short_name': 'nbp',
        'fx_files': ['sftlf'],
    }]

    @staticmethod
    def calculate(cubes):
        """Compute net biome production relative to grid cell area.

        Note
        ----
        By default, `nbp` is defined relative to land area. For easy spatial
        integration, the original quantity is multiplied by the land area
        fraction (`sftlf`), so that the resuting derived variable is defined
        relative to the grid cell area. This correction is only relevant for
        coastal regions.
        """
        nbp_cube = cubes.extract_strict(
            Constraint(name='surface_net_downward_mass_flux_of_carbon_dioxide_'
                       'expressed_as_carbon_due_to_all_land_processes'))
        try:
            sftlf_cube = cubes.extract_strict(
                Constraint(name='land_area_fraction'))
            nbp_cube.data *= sftlf_cube.data / 100.0
        except iris.exceptions.ConstraintMismatchError:
            pass
        return nbp_cube
