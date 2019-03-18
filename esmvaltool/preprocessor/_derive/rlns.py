"""Derivation of variable `rlns`."""

from iris import Constraint

from ._baseclass import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `rlns`."""

    # Required variables
    required = [
        {
            'short_name': 'rlds'
        },
        {
            'short_name': 'rlus'
        },
    ]

    @staticmethod
    def calculate(cubes):
        """Compute surface net downward longwave radiation."""
        rlds_cube = cubes.extract_strict(
            Constraint(name='surface_downwelling_longwave_flux_in_air'))
        rlus_cube = cubes.extract_strict(
            Constraint(name='surface_upwelling_longwave_flux_in_air'))

        rlns_cube = rlds_cube - rlus_cube

        return rlns_cube
