"""Fixes for ESA-CCI Landcover."""

import iris

from ..fix import Fix


class baresoilFrac(Fix):
    """Fixes for baresoilFrac."""

    def fix_metadata(self, cube):
        """
        Fix missing scalar dimension.

        Parameters
        ----------
        cube: iris.cube.Cube
            Cube to fix

        Returns
        -------
        iris.cube.Cube

        """
        typebare = iris.coords.AuxCoord(
            'bare_ground',
            standard_name='area_type',
            long_name='surface type',
            var_name='type',
            units='1',
            bounds=None)
        cube.add_aux_coord(typebare)
        return cube
