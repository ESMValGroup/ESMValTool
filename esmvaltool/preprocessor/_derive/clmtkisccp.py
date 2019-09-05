"""Derivation of variable `clmtkisccp`."""

from iris import Constraint

from ._baseclass import DerivedVariableBase
from ._shared import cloud_area_fraction


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `clmtkisccp`."""

    # Required variables
    required = [{'short_name': 'clisccp'}]

    @staticmethod
    def calculate(cubes):
        """Compute ISCCP middle level thick cloud area fraction."""
        tau = Constraint(
            atmosphere_optical_thickness_due_to_cloud=lambda t: t > 23.)
        plev = Constraint(air_pressure=lambda p: 44000. < p <= 68000.)

        return cloud_area_fraction(cubes, tau, plev)
