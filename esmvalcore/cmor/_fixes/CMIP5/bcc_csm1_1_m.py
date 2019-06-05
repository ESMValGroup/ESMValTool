# pylint: disable=invalid-name, no-self-use, too-few-public-methods
"""Fixes for bcc-csm1-1-m."""
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.ndimage import map_coordinates

from ..fix import Fix


class tos(Fix):
    """Fixes for tos."""

    def fix_data(self, cube):
        """Fix data.

        Calculate missing latitude/longitude boundaries using interpolation.

        Parameters
        ----------
        cube: iris.cube.Cube

        Returns
        -------
        iris.cube.Cube

        """
        rlat = cube.coord('grid_latitude').points
        rlon = cube.coord('grid_longitude').points

        # Transform grid latitude/longitude to array indices [0, 1, 2, ...]
        rlat_to_idx = InterpolatedUnivariateSpline(
            rlat, np.arange(len(rlat)), k=1)
        rlon_to_idx = InterpolatedUnivariateSpline(
            rlon, np.arange(len(rlon)), k=1)
        rlat_idx_bnds = rlat_to_idx(cube.coord('grid_latitude').bounds)
        rlon_idx_bnds = rlon_to_idx(cube.coord('grid_longitude').bounds)

        # Calculate latitude/longitude vertices by interpolation
        lat_vertices = []
        lon_vertices = []
        for (i, j) in [(0, 0), (0, 1), (1, 1), (1, 0)]:
            (rlat_v, rlon_v) = np.meshgrid(
                rlat_idx_bnds[:, i], rlon_idx_bnds[:, j], indexing='ij')
            lat_vertices.append(
                map_coordinates(
                    cube.coord('latitude').points, [rlat_v, rlon_v],
                    mode='nearest'))
            lon_vertices.append(
                map_coordinates(
                    cube.coord('longitude').points, [rlat_v, rlon_v],
                    mode='wrap'))
        lat_vertices = np.array(lat_vertices)
        lon_vertices = np.array(lon_vertices)
        lat_vertices = np.moveaxis(lat_vertices, 0, -1)
        lon_vertices = np.moveaxis(lon_vertices, 0, -1)

        # Copy vertices to cube
        cube.coord('latitude').bounds = lat_vertices
        cube.coord('longitude').bounds = lon_vertices

        return cube
