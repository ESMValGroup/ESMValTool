"""Derivation of variable `lwcre`."""


from iris import Constraint

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `lwcre`."""

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
                ('rlutcs', 'T2' + frequency + 's')]

    def calculate(self, cubes):
        """Compute longwave cloud radiative effect.

        Parameters
        ----------
        cubes : iris.cube.CubeList
            `CubeList` containing `rlut` (`toa_outgoing_longwave_flux`) and
            `rlutcs` (`toa_outgoing_longwave_flux_assuming_clear_sky`).

        Returns
        -------
        iris.cube.Cube
            `Cube` containing longwave cloud radiative effect.

        """
        rlut_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_longwave_flux'))
        rlutcs_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_longwave_flux_assuming_clear_sky'))

        lwcre_cube = rlutcs_cube - rlut_cube
        lwcre_cube.units = rlut_cube.units

        return lwcre_cube
