"""Fixes for IPSL-CM6A-LR model."""
# pylint: disable=invalid-name, no-self-use, too-few-public-methods
import iris
from iris.cube import CubeList
from iris.coords import AuxCoord
from iris.exceptions import ConstraintMismatchError

from ..fix import Fix


class allvars(Fix):
    """Fixes for thetao."""

    def fix_metadata(self, cubelist):
        """
        Fix cell_area coordinate.

        Parameters
        ----------
        cubelist: iris CubeList
            List of cubes to fix

        Returns
        -------
        iris.cube.CubeList

        """
        try:
            cell_area = cubelist.extract_strict('cell_area')
        except ConstraintMismatchError:
            return cubelist

        cell_area = AuxCoord(
            cell_area.data,
            standard_name=cell_area.standard_name,
            long_name=cell_area.long_name,
            var_name=cell_area.var_name,
            units=cell_area.units,
        )
        new_list = CubeList()
        for cube in cubelist
            if cube.name == 'cell_area':
                continue
            cube.add_aux_coord(cell_area, (2, 3))
            cube.coord('latitude').var_name = 'lat'
            cube.coord('longitude').var_name = 'lon'
            new_list.append()
        return CubeList(new_list)
