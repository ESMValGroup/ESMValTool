"""Derivation of variable `rsnt`."""


from iris import Constraint

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `rsnt`."""

    # Required variables
    _required_variables = {'vars': [('rsdt', 'T2{frequency}s'),
                                    ('rsut', 'T2{frequency}s')]}

    def calculate(self, cubes):
        """Compute toa net downward shortwave radiation."""
        rsdt_cube = cubes.extract_strict(
            Constraint(name='toa_incoming_shortwave_flux'))
        rsut_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_shortwave_flux'))

        rsnt_cube = rsdt_cube - rsut_cube

        return rsnt_cube
