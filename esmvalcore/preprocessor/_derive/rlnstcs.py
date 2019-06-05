"""Derivation of variable `rlnstcs`.

authors:
    - weig_ka

"""
from iris import Constraint

from ._baseclass import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `rlnstcs`."""

    # Required variables
    required = [
        {'short_name': 'rldscs'},
        {'short_name': 'rlus'},
        {'short_name': 'rlutcs'},
    ]

    @staticmethod
    def calculate(cubes):
        """
        Compute variable `rlnstcs`.

        Compute Net Atmospheric Longwave Cooling
        to surface and outer space assuming clear sky.
        """
        rldscs_cube = cubes.extract_strict(
            Constraint(name='surface_downwelling_longwave_flux_in_air_' +
                       'assuming_clear_sky'))
        rlus_cube = cubes.extract_strict(
            Constraint(name='surface_upwelling_longwave_flux_in_air'))
        rlutcs_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_longwave_flux_assuming_clear_sky'))

        rlnstcs_cube = rlutcs_cube + (rldscs_cube - rlus_cube)

        return rlnstcs_cube
