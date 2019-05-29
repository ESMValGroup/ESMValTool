"""Derivation of variable `rsns`."""

from iris import Constraint

from ._baseclass import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `rsns`."""

    # Required variables
    required = [
        {
            'short_name': 'rsds'
        },
        {
            'short_name': 'rsus'
        },
    ]

    @staticmethod
    def calculate(cubes):
        """Compute surface net downward shortwave radiation."""
        rsds_cube = cubes.extract_strict(
            Constraint(name='surface_downwelling_shortwave_flux_in_air'))
        rsus_cube = cubes.extract_strict(
            Constraint(name='surface_upwelling_shortwave_flux_in_air'))

        rsns_cube = rsds_cube - rsus_cube

        return rsns_cube
