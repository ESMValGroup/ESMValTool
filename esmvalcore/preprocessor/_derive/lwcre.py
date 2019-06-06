"""Derivation of variable `lwcre`."""

from iris import Constraint

from ._baseclass import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `lwcre`."""

    # Required variables
    required = [
        {
            'short_name': 'rlut'
        },
        {
            'short_name': 'rlutcs'
        },
    ]

    @staticmethod
    def calculate(cubes):
        """Compute longwave cloud radiative effect."""
        rlut_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_longwave_flux'))
        rlutcs_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_longwave_flux_assuming_clear_sky'))

        lwcre_cube = rlutcs_cube - rlut_cube
        lwcre_cube.units = rlut_cube.units

        return lwcre_cube
