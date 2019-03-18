"""Derivation of variable `rsnt`."""

from iris import Constraint

from ._baseclass import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `rsnt`."""

    # Required variables
    required = [
        {
            'short_name': 'rsdt'
        },
        {
            'short_name': 'rsut'
        },
    ]

    @staticmethod
    def calculate(cubes):
        """Compute toa net downward shortwave radiation."""
        rsdt_cube = cubes.extract_strict(
            Constraint(name='toa_incoming_shortwave_flux'))
        rsut_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_shortwave_flux'))

        rsnt_cube = rsdt_cube - rsut_cube

        return rsnt_cube
