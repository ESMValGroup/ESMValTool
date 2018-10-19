"""Derivation of variable `rlns`."""


from iris import Constraint

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `rlns`."""

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
        return [('rlds', 'T2' + frequency + 's'),
                ('rlus', 'T2' + frequency + 's')]

    def calculate(self, cubes):
        """Compute surface net downward longwave radiation.

        Parameters
        ----------
        cubes : iris.cube.CubeList
            `CubeList` containing `rlds`
            (`surface_downwelling_longwave_flux_in_air`) and `rlus`
            (`surface_upwelling_longwave_flux_in_air`).

        Returns
        -------
        iris.cube.Cube
            `Cube` containing surface net downward longwave radiation.

        """
        rlds_cube = cubes.extract_strict(
            Constraint(name='surface_downwelling_longwave_flux_in_air'))
        rlus_cube = cubes.extract_strict(
            Constraint(name='surface_upwelling_longwave_flux_in_air'))

        rlns_cube = rlds_cube - rlus_cube

        return rlns_cube
