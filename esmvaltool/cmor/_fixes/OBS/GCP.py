"""Fixes for GCP."""
import numpy as np

import iris

from ..fix import Fix


class nbp(Fix):
    """Class to fix nbp."""

    def fix_metadata(self, cube):
        """
        Fix metadata for npb.

        Add latitude and longitude coordinates as GCP data is actually T0M.

        FIXME:
            Better introduce a new variable (e.g. 'nbp_t0m'), but at the
            moment NCL is not able to process multiple variables (#531).
            Best solution right now: rename "T0M" -> "T2Ms" in the filename.
        """
        time_coord = cube.coord('time')
        lat_coord = iris.coords.DimCoord([0],
                                         standard_name='latitude',
                                         long_name='latitude coordinate',
                                         var_name='lat',
                                         units='degrees_north',
                                         bounds=[-90, 90])
        lon_coord = iris.coords.DimCoord([180],
                                         standard_name='longitude',
                                         long_name='longitude coordinate',
                                         var_name='lon',
                                         units='degrees_east',
                                         bounds=[0, 360])
        metadata = cube.metadata

        # Create new cube with latitude and longitude coordinate
        new_data = cube.data[:, np.newaxis, np.newaxis]
        cube = iris.cube.Cube(new_data,
                              dim_coords_and_dims=[(time_coord, 0),
                                                   (lat_coord, 1),
                                                   (lon_coord, 2)],
                              **metadata._asdict())
        cube.attributes['field'] = "T2M"

        return cube
