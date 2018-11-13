"""Fixes for GFDL ESM2G"""
import iris
from iris.coords import AuxCoord
from ..fix import Fix



class allvars(Fix):
    """Common fixes"""

    def fix_metadata(self, cubes):
        """
        Fix metadata

        Fixes bad standard names

        Parameters
        ----------
        cubes: iris.cube.CubeList

        Returns
        -------
        iris.cube.Cube

        """


        start_time = self._get_and_remove(
            cubes, 'Start time for average period'
        )
        end_time = self._get_and_remove(cubes, 'End time for average period')
        length = self._get_and_remove(cubes, 'Length of average period')

        self._add_aux_coord(cubes[0], start_time)
        self._add_aux_coord(cubes[0], end_time)
        self._add_aux_coord(cubes[0], length)

        return cubes

    def _get_and_remove(self, cubes, long_name):
        try:
            cube = cubes.extract_strict(long_name)
            cubes.remove(cube)
            return cube
        except iris.exceptions.ConstraintMismatchError:
            return None


    def _add_aux_coord(self, target, cube):
        if cube is None:
            return
        coordinate = AuxCoord(
            cube.data,
            standard_name=cube.standard_name,
            long_name=cube.long_name,
            var_name=cube.var_name,
            units=cube.units,
        )
        target.add_aux_coord(coordinate, target.coord_dims('time'))

class co2(Fix):
    """Fixes for co2"""

    def fix_data(self, cube):
        """
        Fix data

        Fixes discrepancy between declared units and real units

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        metadata = cube.metadata
        cube *= 1e6
        cube.metadata = metadata
        return cube
