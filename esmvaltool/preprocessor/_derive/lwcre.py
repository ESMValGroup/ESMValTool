"""Derivation of variable `lwcre`."""


from iris import Constraint

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `lwcre`."""

    # Required variables
    _required_variables = {'vars': [{'short_name': 'rlut',
                                     'field': 'T2{frequency}s'},
                                    {'short_name': 'rlutcs',
                                     'field': 'T2{frequency}s'}]}

    def calculate(self, cubes):
        """Compute longwave cloud radiative effect."""
        rlut_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_longwave_flux'))
        rlutcs_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_longwave_flux_assuming_clear_sky'))

        lwcre_cube = rlutcs_cube - rlut_cube
        lwcre_cube.units = rlut_cube.units

        return lwcre_cube
