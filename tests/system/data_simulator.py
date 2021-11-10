"""Simulate test data for `esmvaltool`."""
import os
import sys
import time

import cf_units
import iris
import numpy as np
from esmvalcore._config import read_config_user_file
from esmvalcore._recipe import read_recipe_file


def get_input_filename(variable, rootpath, drs):
    """Get a valid input filename."""
    # TODO PROTOTYPE: implement this according to esmvalcore._data_finder.py
    # or patch get_input_filelist there.
    try:
        filename = (
            f"{variable['short_name']}_{variable['mip']}_"
            f"{variable['dataset']}_{variable['exp']}_{variable['ensemble']}_"
            f"{variable['grid']}_{variable['start_year']}-"
            f"{variable['end_year']}.nc")
        filepath = os.path.join(rootpath['CMIP6'][0], filename)
    except KeyError:
        filename = (
            f"{variable['project']}_{variable['dataset']}_{variable['type']}_"
            f"{variable['version']}_{variable['mip']}_"
            f"{variable['short_name']}_"
            f"{variable['start_year']}-{variable['end_year']}.nc")

        filepath = os.path.join(rootpath['CMIP6'][0], 'Tier3/ERA-Interim',
                                filename)

    return filepath


def dummy_cube(short_name, start_year, end_year):
    """Return dummy cube."""
    # TODO PROTOTYPE: Replace this with code that looks in the MIP table to
    # create the cube.

    # Time coordinate.
    unit_string = f"days since {start_year}-01-01"
    time_unit = cf_units.Unit(unit_string, calendar="360_day")
    number_of_days = np.arange(0, (end_year - start_year) * 30, step=30)
    time_coord = iris.coords.DimCoord(number_of_days,
                                      standard_name="time",
                                      units=time_unit)
    time_len = len(time_coord.points)
    time_coord.guess_bounds()

    # Latitude coordinate.
    lat_coord = iris.coords.DimCoord(np.arange(-87.5, 88, 45),
                                     standard_name="latitude",
                                     var_name="lat",
                                     units="degrees")
    lat_len = len(lat_coord.points)
    lat_coord.guess_bounds()

    # Longitude coordinate.
    lon_coord = iris.coords.DimCoord(np.arange(-175, 176, 90),
                                     standard_name="longitude",
                                     var_name="lon",
                                     units="degrees")
    lon_len = len(lon_coord.points)
    lon_coord.guess_bounds()

    # Dealing with error tas: does not match coordinate rank

    # height_coord
    height_coord = iris.coords.DimCoord(np.array(2.5),
                                        standard_name="height",
                                        var_name="height",
                                        units="m")

    # Data.
    data = np.ones((time_len, lat_len, lon_len))

    # Create cube.
    dim_coords_and_dims = [(time_coord, 0), (lat_coord, 1), (lon_coord, 2)]
    aux_coords_and_dims = [(height_coord, None)]
    cube = iris.cube.Cube(data,
                          dim_coords_and_dims=dim_coords_and_dims,
                          aux_coords_and_dims=aux_coords_and_dims)
    cube.var_name = short_name
    cube.standard_name = 'air_temperature'

    cube.attributes['parent_time_units'] = unit_string
    cube.units = 'K'

    return cube


def write_data_file(short_name, filename, start_year, end_year):
    """Write a file containing simulated data."""
    cube = dummy_cube(short_name, start_year, end_year)
    iris.save(cube, filename)


def simulate_input_data(recipe_file, config_user_file=None):
    """Simulate data for variables defined in recipe."""
    if config_user_file:
        user_config = read_config_user_file(config_file=config_user_file,
                                            folder_name='')
    else:
        user_config = {
            'rootpath': {
                'default': '.',
            },
            'drs': {},
        }

    recipe = read_recipe_file(recipe_file, user_config, initialize_tasks=False)

    start_time = time.time()

    for diagnostic in recipe.diagnostics.values():
        np.random.seed(0)
        for variables in diagnostic['preprocessor_output'].values():
            for variable in variables:
                filename = get_input_filename(variable=variable,
                                              rootpath=user_config['rootpath'],
                                              drs=user_config['drs'])
                dirname = os.path.dirname(filename)
                if not os.path.exists(dirname):
                    print("Creating {}".format(dirname))
                    os.makedirs(dirname)

                print("Writing {}".format(filename))
                write_data_file(
                    short_name=variable['short_name'],
                    filename=filename,
                    start_year=variable['start_year'],
                    end_year=variable['end_year'],
                )

    print("Simulating data took {:.0f} seconds".format(time.time() -
                                                       start_time))


if __name__ == '__main__':
    for path in sys.argv[1:]:
        simulate_input_data(recipe_file=path, config_user_file=None)
