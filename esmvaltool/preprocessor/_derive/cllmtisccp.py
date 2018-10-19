"""Derivation of variable `cllmtisccp`."""


import iris
from iris import Constraint

from ._derived_variable_base import DerivedVariableBase


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `cllmtisccp`."""

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
        return [('clisccp', 'T4' + frequency)]

    def calculate(self, cubes):
        """Compute ISCCP low level medium-thickness cloud area fraction.

        Parameters
        ----------
        cubes : iris.cube.CubeList
            `CubeList` containing `clisccp` (`isccp_cloud_area_fraction`).

        Returns
        -------
        iris.cube.Cube
            `Cube` containing ISCCP low level medium-thickness cloud area
            fraction.

        """
        clisccp_cube = cubes.extract_strict(
            Constraint(name='isccp_cloud_area_fraction'))

        tau = iris.Constraint(
            atmosphere_optical_thickness_due_to_cloud=lambda t: 3.6 < t <= 23.)
        plev = iris.Constraint(air_pressure=lambda p: p > 68000.)
        cllmtisccp_cube = clisccp_cube
        cllmtisccp_cube = cllmtisccp_cube.extract(tau & plev)
        coord_names = [
            coord.standard_name for coord in cllmtisccp_cube.coords()
            if len(coord.points) > 1
        ]
        if 'atmosphere_optical_thickness_due_to_cloud' in coord_names:
            cllmtisccp_cube = cllmtisccp_cube.collapsed(
                'atmosphere_optical_thickness_due_to_cloud', iris.analysis.SUM)
        if 'air_pressure' in coord_names:
            cllmtisccp_cube = cllmtisccp_cube.collapsed('air_pressure',
                                                        iris.analysis.SUM)

        return cllmtisccp_cube
