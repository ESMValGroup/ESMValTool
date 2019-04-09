"""Derivation of variable `xch4`."""

from iris import Constraint

from ._baseclass import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `xch4`."""

    # Required variables
    required = [{
                'short_name': 'xch4'
                }]

    @staticmethod
    def calculate(cubes):
        """Compute/forward xch4."""
        xch4_cube = cubes.extract_strict(
            Constraint(name='column_averaged_methane'))

        return xch4_cube
