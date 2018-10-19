"""Fixes for MPI-ESM1-2-HR CMIP6 project data."""
from iris.cube import CubeList
from ..fix import Fix


class zg(Fix):
    """Fixes common to all variables."""

    def fix_metadata(self, cube):
        """
        Fix cube metadata.

        Parameters
        ----------
        cube: Cube
            Cube to fix

        Returns
        -------
        Cube:
            Fixed cube. It is the same instance that was received
        """
        plev = cube.coord('air_pressure')
        plev.var_name = 'plev'
        slices = CubeList(reversed(
            [lat_slice for lat_slice in cube.slices_over('latitude')]
        ))
        cube = slices.merge_cube()
        return cube
