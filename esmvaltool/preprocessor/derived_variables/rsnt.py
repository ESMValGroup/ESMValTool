"""Derivation of variable `rsnt`."""


from iris import Constraint

from .derived_variable import DerivedVariable


class rsnt(DerivedVariable):  # noqa
    """Derivation of variable `rsnt`."""

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
        return [('rsdt', 'T2' + frequency + 's'),
                ('rsut', 'T2' + frequency + 's')]

    def calculate(self, cubes):
        """Compute toa net downward shortwave radiation.

        Parameters
        ----------
        cubes : iris.cube.CubeList
            `CubeList` containing `rsut` (`toa_outgoing_shortwave_flux`) and
            `rsdt` (`toa_incoming_shortwave_flux`).

        Returns
        -------
        iris.cube.Cube
            `Cube` containing toa net downward shortwave radiation.

        """
        rsdt_cube = cubes.extract_strict(
            Constraint(name='toa_incoming_shortwave_flux'))
        rsut_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_shortwave_flux'))

        rsnt_cube = rsdt_cube - rsut_cube

        return rsnt_cube
