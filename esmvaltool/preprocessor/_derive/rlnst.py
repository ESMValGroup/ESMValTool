"""Derivation of variable `rlnst`.

authors:
    - weig_ka

"""
from iris import Constraint

from ._baseclass import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `rlnst`."""

    # Required variables
    required = [
        {'short_name': 'rlds'},
        {'short_name': 'rlus'},
        {'short_name': 'rlut'},
    ]

    @staticmethod
    def calculate(cubes):
        """
        Compute variable `rlnst`.

        Compute Net Atmospheric Longwave Cooling
        to surface and outer space.
        """
        rlds_cube = cubes.extract_strict(
            Constraint(name='surface_downwelling_longwave_flux_in_air'))
        rlus_cube = cubes.extract_strict(
            Constraint(name='surface_upwelling_longwave_flux_in_air'))
        rlut_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_longwave_flux'))

        rlnst_cube = rlut_cube + (rlds_cube - rlus_cube)

        return rlnst_cube
