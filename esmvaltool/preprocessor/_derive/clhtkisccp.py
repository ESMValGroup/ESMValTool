"""Derivation of variable `clhtkisccp`."""


from iris import Constraint

from ._shared import cloud_area_fraction
from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `clhtkisccp`."""

    # Required variables
    _required_variables = {'vars': [{'short_name': 'clisccp',
                                     'field': 'T4{frequency}'}]}

    def calculate(self, cubes):
        """Compute ISCCP high level thick cloud area fraction."""
        tau = Constraint(
            atmosphere_optical_thickness_due_to_cloud=lambda t: t > 23.)
        plev = Constraint(air_pressure=lambda p: p <= 44000.)

        return cloud_area_fraction(cubes, tau, plev)
