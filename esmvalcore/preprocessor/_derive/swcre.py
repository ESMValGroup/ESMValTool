"""Derivation of variable `swcre`."""

from iris import Constraint

from ._baseclass import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `swcre`."""

    # Required variables
    required = [
        {
            'short_name': 'rsut'
        },
        {
            'short_name': 'rsutcs'
        },
    ]

    @staticmethod
    def calculate(cubes):
        """Compute shortwave cloud radiative effect."""
        rsut_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_shortwave_flux'))
        rsutcs_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_shortwave_flux_assuming_clear_sky'))

        swcre_cube = rsutcs_cube - rsut_cube

        return swcre_cube
