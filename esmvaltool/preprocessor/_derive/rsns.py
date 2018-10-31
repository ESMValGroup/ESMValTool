"""Derivation of variable `rsns`."""


from iris import Constraint

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `rsns`."""

    # Required variables
    _required_variables = {'vars': [{'short_name': 'rsds',
                                     'field': 'T2{frequency}s'},
                                    {'short_name': 'rsus',
                                     'field': 'T2{frequency}s'}]}

    def calculate(self, cubes):
        """Compute surface net downward shortwave radiation."""
        rsds_cube = cubes.extract_strict(
            Constraint(name='surface_downwelling_shortwave_flux_in_air'))
        rsus_cube = cubes.extract_strict(
            Constraint(name='surface_upwelling_shortwave_flux_in_air'))

        rsns_cube = rsds_cube - rsus_cube

        return rsns_cube
