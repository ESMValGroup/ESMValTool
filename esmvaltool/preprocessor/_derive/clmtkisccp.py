"""Derivation of variable `clmtkisccp`."""


from iris import Constraint

from ._shared import cloud_area_fraction
from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `clmtkisccp`."""

    # Required variables
    _required_variables = {'vars': [('clisccp', 'T4{frequency}')]}

    def calculate(self, cubes):
        """Compute ISCCP middle level thick cloud area fraction."""
        tau = Constraint(
            atmosphere_optical_thickness_due_to_cloud=lambda t: t > 23.)
        plev = Constraint(air_pressure=lambda p: 44000. < p <= 68000.)

        return cloud_area_fraction(cubes, tau, plev)
