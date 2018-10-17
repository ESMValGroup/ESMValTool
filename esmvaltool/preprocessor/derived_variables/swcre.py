"""Derivation of variable `swcre`."""


from iris import Constraint

from ._derived_variable import DerivedVariable


class swcre(DerivedVariable):  # noqa
    """Derivation of variable `swcre`."""

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
        return [('rsut', 'T2' + frequency + 's'),
                ('rsutcs', 'T2' + frequency + 's')]

    def calculate(self, cubes):
        """Compute shortwave cloud radiative effect.

        Parameters
        ----------
        cubes : iris.cube.CubeList
            `CubeList` containing `rsut` (`toa_outgoing_shortwave_flux`) and
            `rsutcs` (`toa_outgoing_shortwave_flux_assuming_clear_sky`).

        Returns
        -------
        iris.cube.Cube
            `Cube` containing shortwave cloud radiative effect.

        """
        rsut_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_shortwave_flux'))
        rsutcs_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_shortwave_flux_assuming_clear_sky'))

        swcre_cube = rsutcs_cube - rsut_cube

        return swcre_cube
