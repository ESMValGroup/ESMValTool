"""Derivation of variable `alb`.

authors:
    - crez_ba

"""

from iris import Constraint

from ._baseclass import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `alb`."""

    # Required variables
    required = [
        {
            'short_name': 'rsds'
        },
        {
            'short_name': 'rsus'
        },
    ]

    @staticmethod
    def calculate(cubes):
        """Compute surface albedo."""
        rsds_cube = cubes.extract_strict(
            Constraint(name='surface_downwelling_shortwave_flux_in_air'))
        rsus_cube = cubes.extract_strict(
            Constraint(name='surface_upwelling_shortwave_flux_in_air'))

        rsns_cube = rsus_cube / rsds_cube

        return rsns_cube
