#! /usr/local/sci/bin/python
"""
Port for ESMValTool v2 from v1.

Uses: ESMValTool v2, Python 3.x
Valeriu Predoi, UREAD, July 2018

The script is well different than the v1 vresion but executes the
same set of functionalities; script name kept the same as in v1
for historical purposes.
"""

import iris
import iris.analysis.cartography


def get_cube_ready(cube):
    """Remve unwanted coords and check bounds."""
    to_remove_list = [
        'forecast_reference_time', 'forecast_period', 'source', 'season',
        'time'
    ]
    for coord in cube.coords():
        if coord.name() in to_remove_list:
            cube.remove_coord(coord)
    if not cube.coord(axis='x').has_bounds():
        cube.coord(axis='x').guess_bounds()
    if not cube.coord(axis='y').has_bounds():
        cube.coord(axis='y').guess_bounds()

    return cube


def area_avg(cube, coord1=None, coord2=None):
    """
    Get area average.

    Perform an area average of a cube using weights to account for
    changes in latitude.
    """
    for coord in (coord1, coord2):
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(cube)
    result = cube.collapsed(
        [coord1, coord2], iris.analysis.MEAN, weights=grid_areas)

    return result


def perform_equation(dataset_1, dataset_2, analysis_type):
    """
    Perform a simple cube operation.

    analysis_type = type of analysis (zonal_mean, vertical_mean,...)
    This can be easily adapted for more than one type of operation
    by passing an argument e.g. 'sum_of_squares' etc.
    """
    # Make sure all the fields have correct units
    dataset_1_ready = get_cube_ready(dataset_1)
    dataset_2_ready = get_cube_ready(dataset_2)

    if analysis_type == 'zonal_mean':
        dataset_1_mean = dataset_1_ready.collapsed('longitude',
                                                   iris.analysis.MEAN)
        dataset_2_mean = dataset_2_ready.collapsed('longitude',
                                                   iris.analysis.MEAN)

    elif analysis_type == 'vertical_mean':
        dataset_1_mean = dataset_1_ready.collapsed('pressure',
                                                   iris.analysis.MEAN)
        dataset_2_mean = dataset_2_ready.collapsed('pressure',
                                                   iris.analysis.MEAN)
    elif analysis_type == 'lat_lon':
        dataset_1_mean = dataset_1_ready
        dataset_2_mean = dataset_2_ready

    # Perform simple difference
    toplot_cube = dataset_1_mean - dataset_2_mean

    return toplot_cube
