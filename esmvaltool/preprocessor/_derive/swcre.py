"""Derivation of variable `swcre`."""


from iris import Constraint

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `swcre`."""

    # Required variables
    _required_variables = {'vars': [('rsut', 'T2{frequency}s'),
                                    ('rsutcs', 'T2{frequency}s')]}

    def calculate(self, cubes):
        """Compute shortwave cloud radiative effect."""
        rsut_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_shortwave_flux'))
        rsutcs_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_shortwave_flux_assuming_clear_sky'))

        swcre_cube = rsutcs_cube - rsut_cube

        return swcre_cube
