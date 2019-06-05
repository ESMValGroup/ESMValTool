"""Derivation of variable `rsnstcs`.

authors:
    - weig_ka

"""
from iris import Constraint

from ._baseclass import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `rsnstcs`."""

    # Required variables
    required = [
        {'short_name': 'rsdscs'},
        {'short_name': 'rsdt'},
        {'short_name': 'rsuscs'},
        {'short_name': 'rsutcs'},
    ]

    @staticmethod
    def calculate(cubes):
        """Compute Heating from Shortwave Absorption assuming clear sky."""
        rsdscs_cube = cubes.extract_strict(
            Constraint(name='surface_downwelling_shortwave_flux_in_air_' +
                       'assuming_clear_sky'))
        rsdt_cube = cubes.extract_strict(
            Constraint(name='toa_incoming_shortwave_flux'))
        rsuscs_cube = cubes.extract_strict(
            Constraint(name='surface_upwelling_shortwave_flux_in_air_' +
                       'assuming_clear_sky'))
        rsutcs_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_shortwave_flux_assuming_clear_sky'))

        rsnstcs_cube = (rsdt_cube - rsutcs_cube) - (rsdscs_cube - rsuscs_cube)

        return rsnstcs_cube
