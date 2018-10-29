"""Derivation of variable `rlns`."""


from iris import Constraint

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `rlns`."""

    # Required variables
    _required_variables = {'vars': [('rlds', 'T2{frequency}s'),
                                    ('rlus', 'T2{frequency}s')]}

    def calculate(self, cubes):
        """Compute surface net downward longwave radiation."""
        rlds_cube = cubes.extract_strict(
            Constraint(name='surface_downwelling_longwave_flux_in_air'))
        rlus_cube = cubes.extract_strict(
            Constraint(name='surface_upwelling_longwave_flux_in_air'))

        rlns_cube = rlds_cube - rlus_cube

        return rlns_cube
