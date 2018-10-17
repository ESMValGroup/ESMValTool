"""Derivation of variable `netcre`."""


from ._derived_variable import DerivedVariable
from .lwcre import lwcre
from .swcre import swcre


class netcre(DerivedVariable):  # noqa
    """Derivation of variable `netcre`."""

    def get_required(self, frequency):
        """Get variable `short_name` and `field` pairs required for derivation.

        Parameters
        ----------
        frequency : str
            Frequency of the desired derived variable.

        Returns
        -------
        list of tuples
            List of tuples (`short_name`, `field`) of all variables required
            for derivation.

        """
        return [('rlut', 'T2' + frequency + 's'),
                ('rlutcs', 'T2' + frequency + 's'),
                ('rsut', 'T2' + frequency + 's'),
                ('rsutcs', 'T2' + frequency + 's')]

    def calculate(self, cubes):
        """Compute net cloud radiative effect.

        Calculate net cloud radiative effect as sum of longwave and shortwave
        cloud radiative effects.

        Parameters
        ----------
        cubes : iris.cube.CubeList
            `CubeList` containing `rlut` (`toa_outgoing_longwave_flux`),
            `rlutcs` (`toa_outgoing_longwave_flux_assuming_clear_sky`),
            `rsut` (`toa_outgoing_shortwave_flux`) and `rsutcs`
            (`toa_outgoing_shortwave_flux_assuming_clear_sky`).

        Returns
        -------
        iris.cube.Cube
            `Cube` containing net cloud radiative effect.

        """
        lwcre_var = lwcre()
        swcre_var = swcre()
        lwcre_cube = lwcre_var.calculate(cubes)
        swcre_cube = swcre_var.calculate(cubes)

        netcre_cube = lwcre_cube + swcre_cube
        netcre_cube.units = lwcre_cube.units

        return netcre_cube
