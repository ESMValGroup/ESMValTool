"""Derivation of variable `rsns`."""


from iris import Constraint

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `rsns`."""

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
        return [('rsds', 'T2' + frequency + 's'),
                ('rsus', 'T2' + frequency + 's')]

    def calculate(self, cubes):
        """Compute surface net downward shortwave radiation.

        Parameters
        ----------
        cubes : iris.cube.CubeList
            `CubeList` containing `rsus`
            (`surface_upwelling_shortwave_flux_in_air`) and `rsds`
            (`surface_downwelling_shortwave_flux_in_air`).

        Returns
        -------
        iris.cube.Cube
            `Cube` containing surface net downward shortwave radiation.

        """
        rsds_cube = cubes.extract_strict(
            Constraint(name='surface_downwelling_shortwave_flux_in_air'))
        rsus_cube = cubes.extract_strict(
            Constraint(name='surface_upwelling_shortwave_flux_in_air'))

        rsns_cube = rsds_cube - rsus_cube

        return rsns_cube
