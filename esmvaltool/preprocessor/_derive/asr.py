"""Derivation of variable `asr`."""


from iris import Constraint

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `asr`."""

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

    def calculate(self, cubes, fx_files=None):
        """Compute absorbed shortwave radiation.

        Parameters
        ----------
        cubes : iris.cube.CubeList
            `CubeList` containing `rsdt` (`toa_incoming_shortwave_flux`) and
            `rsut` (`toa_outgoing_shortwave_flux`).
        fx_files : dict, optional
            If required, dictionary containing fx files  with `short_name`
            (key) and path (value) of the fx variable.

        Returns
        -------
        iris.cube.Cube
            `Cube` containing absorbed shortwave radiation.

        """
        rsdt_cube = cubes.extract_strict(
            Constraint(name='toa_incoming_shortwave_flux'))
        rsut_cube = cubes.extract_strict(
            Constraint(name='toa_outgoing_shortwave_flux'))

        asr_cube = rsdt_cube - rsut_cube
        asr_cube.units = rsdt_cube.units

        return asr_cube
