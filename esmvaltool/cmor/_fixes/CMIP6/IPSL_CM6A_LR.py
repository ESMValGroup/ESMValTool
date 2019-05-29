"""Fixes for IPSL-CM6A-LR model."""
# pylint: disable=invalid-name, no-self-use, too-few-public-methods
import iris
from iris.cube import CubeList
from iris.coords import AuxCoord

from ..fix import Fix


class thetao(Fix):
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
        thetao = cubelist.extract_strict('sea_water_potential_temperature')
        cell_area = cubelist.extract_strict('cell_area')
        cell_area = AuxCoord(
            cell_area.data,
            standard_name=cell_area.standard_name,
            long_name=cell_area.long_name,
            var_name=cell_area.var_name,
            units=cell_area.units,
        )
        thetao.add_aux_coord(cell_area, (2, 3))
        thetao.coord('latitude').var_name = 'lat'
        thetao.coord('longitude').var_name = 'lon'
        return CubeList([thetao])
