"""Derivation of variable `netcre`."""


from ._derived_variable_base import DerivedVariableBase
from .lwcre import DerivedVariable as Lwcre
from .swcre import DerivedVariable as Swcre


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `netcre`."""

    # Required variables
    _required_variables = {'vars': [{'short_name': 'rlut',
                                     'field': 'T2{frequency}s'},
                                    {'short_name': 'rlutcs',
                                     'field': 'T2{frequency}s'},
                                    {'short_name': 'rsut',
                                     'field': 'T2{frequency}s'},
                                    {'short_name': 'rsutcs',
                                     'field': 'T2{frequency}s'}]}

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
