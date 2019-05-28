"""Derivation of variable `lvp`.

authors:
    - weig_ka

"""
from iris import Constraint

from ._baseclass import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `lvp`."""

    # Required variables
    required = [
        {'short_name': 'hfls'},
        {'short_name': 'pr'},
        {'short_name': 'evspsbl'},
    ]

    @staticmethod
    def calculate(cubes):
        """Compute Latent Heat Release from Precipitation."""
        hfls_cube = cubes.extract_strict(
            Constraint(name='surface_upward_latent_heat_flux'))
        pr_cube = cubes.extract_strict(Constraint(name='precipitation_flux'))
        evspsbl_cube = cubes.extract_strict(
            Constraint(name='water_evaporation_flux'))

        lvp_cube = hfls_cube * (pr_cube / evspsbl_cube)

        return lvp_cube
