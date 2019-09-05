"""Derivation of variable `netcre`."""

from ._baseclass import DerivedVariableBase
from .lwcre import DerivedVariable as Lwcre
from .swcre import DerivedVariable as Swcre


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `netcre`."""

    # Required variables
    required = [
        {
            'short_name': 'rlut'
        },
        {
            'short_name': 'rlutcs'
        },
        {
            'short_name': 'rsut'
        },
        {
            'short_name': 'rsutcs'
        },
    ]

    @staticmethod
    def calculate(cubes):
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
