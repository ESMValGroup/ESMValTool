"""Derivation of variable `rtnt`."""


from iris import Constraint

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `rtnt`."""

    # Required variables
    _required_variables = {'vars': [{'short_name': 'rsdt',
                                     'field': 'T2{frequency}s'},
                                    {'short_name': 'rsut',
                                     'field': 'T2{frequency}s'},
                                    {'short_name': 'rlut',
                                     'field': 'T2{frequency}s'}]}

    def calculate(self, cubes):
        """Compute toa net downward total radiation."""
        rsdt_cube = cubes.extract_strict(
            Constraint(name='toa_incoming_shortwave_flux'))
        rsut_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_shortwave_flux'))
        rlut_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_longwave_flux'))

        rtnt_cube = rsdt_cube - rsut_cube - rlut_cube

        return rtnt_cube
