"""Derivation of variable `netcre`."""


from ._derived_variable_base import DerivedVariableBase
from .lwcre import DerivedVariable as Lwcre
from .swcre import DerivedVariable as Swcre


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `netcre`."""

    # Required variables
    _required_variables = {'vars': [('rlut', 'T2{frequency}s'),
                                    ('rlutcs', 'T2{frequency}s'),
                                    ('rsut', 'T2{frequency}s'),
                                    ('rsutcs', 'T2{frequency}s')]}

    def calculate(self, cubes):
        """Compute net cloud radiative effect.

        Note
        ----
        Calculate net cloud radiative effect as sum of longwave and shortwave
        cloud radiative effects.
        """
        lwcre_var = Lwcre()
        swcre_var = Swcre()
        lwcre_cube = lwcre_var.calculate(cubes)
        swcre_cube = swcre_var.calculate(cubes)

        netcre_cube = lwcre_cube + swcre_cube
        netcre_cube.units = lwcre_cube.units

        return netcre_cube
