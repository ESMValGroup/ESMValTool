"""Derivation of variable `nbp_grid`."""
from ._baseclass import DerivedVariableBase
from ._shared import grid_area_correction


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
        By default, `nbp` is defined relative to land area. For spatial
        integration, the original quantity is multiplied by the land area
        fraction (`sftlf`), so that the resuting derived variable is defined
        relative to the grid cell area. This correction is only relevant for
        coastal regions.

        """
        return grid_area_correction(
            cubes,
            'surface_net_downward_mass_flux_of_carbon_dioxide_expressed_as_'
            'carbon_due_to_all_land_processes')
