"""Derivation of variable `swcre`."""


from iris import Constraint

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
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

    def calculate(self, cubes, fx_files=None):
        """Compute shortwave cloud radiative effect.

        Parameters
        ----------
        cubes : iris.cube.CubeList
            `CubeList` containing `rsut` (`toa_outgoing_shortwave_flux`) and
            `rsutcs` (`toa_outgoing_shortwave_flux_assuming_clear_sky`).
        fx_files : dict, optional
            If required, dictionary containing fx files  with `short_name`
            (key) and path (value) of the fx variable.

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
