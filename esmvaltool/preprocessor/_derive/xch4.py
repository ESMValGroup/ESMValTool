"""Derivation of variable `xch4`."""

from iris import Constraint

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `xch4`."""

    # Required variables
    _required_variables = {
        'vars': [{
            'short_name': 'xch4',
            'field': 'T2{frequency}s'
        }]
    }

    def calculate(self, cubes):
        """Compute/forward xch4."""
        xch4_cube = cubes.extract_strict(
            Constraint(name='column_averaged_methane'))

        return xch4_cube
