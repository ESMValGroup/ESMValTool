"""Fixes for GCP."""
import iris
import numpy as np

from ..fix import Fix


class nbp(Fix):
    """Class to fix nbp."""

    def fix_metadata(self, cubes):
        """Fix metadata for npb.

        Add latitude and longitude coordinates as GCP data is actually T0M.

        """
        cube = self.get_cube_from_list(cubes)
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
        cube = iris.cube.Cube(
            new_data,
            dim_coords_and_dims=[
                (time_coord, 0),
                (lat_coord, 1),
                (lon_coord, 2),
            ],
            **metadata._asdict(),
        )
        cube.attributes['field'] = 'T2Ms'
        return [cube]


class fgco2(Fix):
    """Class to fix fgco2."""

    def fix_metadata(self, cubes):
        """Fix metadata for fgco2.

        Add latitude and longitude coordinates as GCP data is actually T0M.

        """
        cube = self.get_cube_from_list(cubes)
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
        cube = iris.cube.Cube(
            new_data,
            dim_coords_and_dims=[
                (time_coord, 0),
                (lat_coord, 1),
                (lon_coord, 2),
            ],
            **metadata._asdict(),
        )
        cube.attributes['field'] = 'T2Ms'
        return [cube]
