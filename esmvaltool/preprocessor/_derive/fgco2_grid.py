"""Derivation of variable `fgco2_grid`."""

import logging

from ._derived_variable_base import DerivedVariableBase
from ._shared import grid_area_correction

logger = logging.getLogger(__name__)


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `fgco2_grid`."""

    # Required variables
    _required_variables = {
        'vars': [{
            'short_name': 'fgco2',
            'field': 'T2{frequency}s'
        }],
        'fx_files': ['areacello', 'sftof']
    }

    def calculate(self, cubes):
        """Compute gas exchange flux of CO2 per grid cell.

        Note
        ----
        By default, `fgco2` is defined relative to sea area. For spatial
        integration, the original quantity is multiplied by the sea area
        fraction (`sftof`), so that the resuting derived variable is defined
        relative to the grid cell area. This correction is only relevant for
        coastal regions.

        """
        return grid_area_correction(
            cubes,
            'surface_downward_mass_flux_of_carbon_dioxide_expressed_as_carbon',
            fraction_var='sea_area_fraction')
