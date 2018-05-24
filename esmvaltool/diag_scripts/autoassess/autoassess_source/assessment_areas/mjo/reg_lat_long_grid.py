"Create cube with regular lat-long grid"

from __future__ import division
import numpy as np
import iris


def create_cube(latitudes=(-90, 90), longitudes=(0, 360), spacing=1):
    """
    Create empty cube with regular lat-long grid.

    Args:
        latitudes: Range of latitudes in degree
        longitudes: Range of longitudes in degree
        spacing: lat-long grid spacing in degree

    Returns:
        Empty cube with specified grid spacing, covering specifed
        latitude and longitude intervals.
    """
    # TODO most used coord system?
    cs = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)

    # lat coord
    lower_bound = latitudes[0] + spacing / 2.
    upper_bound = latitudes[1] - spacing / 2.
    number_of_points = (upper_bound - lower_bound) / spacing + 1
    lat_coord = iris.coords.DimCoord(np.linspace(lower_bound, upper_bound, number_of_points),
                                     standard_name='latitude',
                                     units='degrees',
                                     coord_system=cs)
    lat_coord.guess_bounds()

    # long coord
    lower_bound = longitudes[0] + spacing / 2.
    upper_bound = longitudes[1] - spacing / 2.
    number_of_points = (upper_bound - lower_bound) / spacing + 1
    lon_coord = iris.coords.DimCoord(np.linspace(lower_bound, upper_bound, number_of_points),
                                     standard_name='longitude',
                                     units='degrees',
                                     coord_system=cs)
    lon_coord.guess_bounds()

    # data
    data = np.zeros((len(lat_coord.points), len(lon_coord.points)))

    # build cube
    cube = iris.cube.Cube(
        data,
        long_name='zeros',
        units='',
        attributes=None,
        dim_coords_and_dims=[(lat_coord, 0), (lon_coord, 1)]
        )

    return cube


if __name__ == '__main__':
    cube = create_cube()
    print cube
    print cube.coord_system()
    print cube.coord('latitude')
    print cube.coord('longitude')
