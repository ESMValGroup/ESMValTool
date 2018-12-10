"""Derivation of variable `nbp_global`."""

import iris
from iris import Constraint

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `nbp_global`."""

    # Required variables
    _required_variables = {
        'vars': [{
            'short_name': 'nbp',
            'field': 'T2{frequency}s'
        }],
        'fx_files': ['sftlf']
    }

    def calculate(self, cubes):
        """Compute globally averaged net biome production.

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

        # Global average if necessary
        if ('latitude' in nbp_cube.coords()
                and 'longitude' in nbp_cube.coords()):
            for coord in ('latitude', 'longitude'):
                if not nbp_cube.coord(coord).has_bounds():
                    nbp_cube.coord(coord).guess_bounds()
            area_weights = iris.analysis.cartography.area_weights(nbp_cube)
            nbp_cube = nbp_cube.collapsed(['latitude', 'longitude'],
                                          iris.analysis.MEAN,
                                          weights=area_weights)
        return nbp_cube
