"""Derivation of variable `rtnt`."""


from iris import Constraint

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `rtnt`."""

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
                ('rsut', 'T2' + frequency + 's'),
                ('rlut', 'T2' + frequency + 's')]

    def calculate(self, cubes):
        """Compute toa net downward total radiation.

        Parameters
        ----------
        cubes : iris.cube.CubeList
            `CubeList` containing `rsut` (`toa_outgoing_shortwave_flux`),
            `rsdt` (`toa_incoming_shortwave_flux`) and `rlut`
            (`toa_outgoing_longwave_flux`).

        Returns
        -------
        iris.cube.Cube
            `Cube` containing toa net downward total radiation.

        """
        rsdt_cube = cubes.extract_strict(
            Constraint(name='toa_incoming_shortwave_flux'))
        rsut_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_shortwave_flux'))
        rlut_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_longwave_flux'))

        rtnt_cube = rsdt_cube - rsut_cube - rlut_cube

        return rtnt_cube
