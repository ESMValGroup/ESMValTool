"""Derivation of variable `fgco2_grid`."""


import logging

import iris
from iris import Constraint

from ._derived_variable_base import DerivedVariableBase

logger = logging.getLogger(__name__)


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `fgco2_grid`."""

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
        return [('fgco2', 'T2' + frequency + 's'),
                ('fx_files', ['sftof'])]

    def calculate(self, cubes, fx_files=None):
        """Compute gas exchange flux of CO2 relative to grid cell area.

        By default, `fgco2` is defined relative to land area. For easy spatial
        integration, the original quantity is multiplied by the land area
        fraction (`sftlf`), so that the resuting derived variable is defined
        relative to the grid cell area. This correction is only relevant for
        coastal regions.

        Parameters
        ----------
        cubes : iris.cube.CubeList
            `CubeList` containing `fgco2`(`surface_downward_mass_flux_of_
            carbon_dioxide_expressed_as_carbon`).
        fx_files : dict, optional
            If required, dictionary containing fx files  with `short_name`
            (key) and path (value) of the fx variable.

        Returns
        -------
        iris.cube.Cube
            `Cube` containing gas exchange flux of CO2 relative to grid cell
            area.

        """
        fgco2_cube = cubes.extract_strict(
            Constraint(name='surface_downward_mass_flux_of_carbon_dioxide_'
                            'expressed_as_carbon'))
        if fx_files.get('sftof'):
            sftof_cube = iris.load_cube(fx_files['sftof'])
            fgco2_cube.data = fgco2_cube.data * sftof_cube.data / 100.0
        return fgco2_cube
