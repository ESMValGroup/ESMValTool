"""Derivation of variable `lai_grid`."""


import logging

import iris
from iris import Constraint

from ._derived_variable_base import DerivedVariableBase

logger = logging.getLogger(__name__)


class DerivedVariable(DerivedVariableBase):
    """Derivation of variable `lai_grid`."""

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
        return [('lai', 'T2' + frequency + 's'),
                ('fx_files', ['sftlf'])]

    def calculate(self, cubes, fx_files=None):
        """Compute leaf area index relative to grid cell area.

        By default, `lai` is defined relative to land area. For easy spatial
        integration, the original quantity is multiplied by the land area
        fraction (`sftlf`), so that the resuting derived variable is defined
        relative to the grid cell area. This correction is only relevant for
        coastal regions.

        Parameters
        ----------
        cubes : iris.cube.CubeList
            `CubeList` containing `lai` (`leaf_area_index`).
        fx_files : dict, optional
            If required, dictionary containing fx files  with `short_name`
            (key) and path (value) of the fx variable.

        Returns
        -------
        iris.cube.Cube
            `Cube` containing leaf area index relative to grid cell area.

        """
        lai_cube = cubes.extract_strict(
            Constraint(name='leaf_area_index'))
        if fx_files.get('sftlf'):
            sftlf_cube = iris.load_cube(fx_files['sftlf'])
            lai_cube.data = lai_cube.data * sftlf_cube.data / 100.0
        return lai_cube
