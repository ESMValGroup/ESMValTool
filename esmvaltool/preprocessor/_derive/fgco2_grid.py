"""Derivation of variable `fgco2_grid`."""

import logging

import iris
from iris import Constraint

from ._derived_variable_base import DerivedVariableBase

logger = logging.getLogger(__name__)


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `fgco2_grid`."""

    # Required variables
    _required_variables = {
        'vars': [{
            'short_name': 'fgco2',
            'field': 'T2{frequency}s'
        }],
        'fx_files': ['sftof']
    }

    def calculate(self, cubes):
        """Compute gas exchange flux of CO2 relative to grid cell area.

        Note
        ----
        By default, `fgco2` is defined relative to sea area. For easy spatial
        integration, the original quantity is multiplied by the sea area
        fraction (`sftof`), so that the resuting derived variable is defined
        relative to the grid cell area. This correction is only relevant for
        coastal regions.

        """
        fgco2_cube = cubes.extract_strict(
            Constraint(name='surface_downward_mass_flux_of_carbon_dioxide_'
                       'expressed_as_carbon'))
        try:
            sftof_cube = cubes.extract_strict(
                Constraint(name='sea_area_fraction'))
            fgco2_cube.data *= sftof_cube.data / 100.0
        except iris.exceptions.ConstraintMismatchError:
            pass
        return fgco2_cube
