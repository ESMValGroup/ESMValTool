"""Fixes for IPSL-CM6A-LR."""

from ..fix import Fix


class allvars(Fix):
    """Fixes for all vars."""

    def fix_metadata(self, cubelist):
        """
        Removes addition cube in cube list (Cell area).
        It has been created as a cell_area_cube which can be re-added here,
        if needed.

        Also forced the latitude and longitude to be 'lat' & 'lon',
        instead of 'nav_lat' & 'nav_lon'.

        Parameters
        ----------
        cubelist: iris CubeList
            List of cubes to fix

        Returns
        -------
        iris.cube.CubeList
        """
        for cube in cubelist:
            if cube.standard_name == 'cell_area':
                # If you need to add the cell area back in, uncomment this line
                # cell_area_cube = cube.copy()
                #
                cubelist.remove(cube)

        for cube in cubelist:
            lats = cube.coord('latitude')
            lons = cube.coord('longitude')
            lats.var_name = 'lat'
            lons.var_name = 'lon'
        return cubelist
