"""Simulate test data for `esmvaltool`."""
import os
import sys
import tempfile
import time

import numpy as np

from esmvalcore._config import read_config_user_file
from esmvalcore._recipe import read_recipe_file


def get_input_filename(variable, rootpath, drs):
    """Get a valid input filename."""
    # TODO: implement this according to esmvalcore._data_finder.py
    # or patch get_input_filelist there.
    return tempfile.NamedTemporaryFile().name + '.nc'


def write_data_file(short_name, filename, field, start_year, end_year):
    """Write a file containing simulated data."""
    from dummydata.model2 import Model2
    from dummydata.model3 import Model3

    if 'T2M' in field:
        writer = Model2
    elif 'T3M' in field:
        writer = Model3
    else:
        raise NotImplementedError(
            "Cannot create a model from field {}".format(field))

    # TODO: Maybe this should be made configurable per diagnostic or model
    cfg = {
        'ta': {
            'method': 'gaussian_blobs',
            'low': 223,
            'high': 303,
        },
        'pr': {
            'method': 'gaussian_blobs',
            'low': 1e-7,
            'high': 2e-4,
        }
    }

    kwargs = cfg[short_name] if short_name in cfg else {}

    writer(
        var=short_name,
        oname=filename,
        start_year=start_year,
        stop_year=end_year,
        **kwargs)


def simulate_input_data(recipe_file, config_user_file=None):
    """Simulate data for variables defined in recipe"""
    if config_user_file:
        user_config = read_config_user_file(
            config_file=config_user_file, recipe_name='')
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
        for variables in diagnostic['variables'].values():
            for variable in variables:
                filename = get_input_filename(
                    variable=variable,
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
                    field=variable['field'],
                    start_year=variable['start_year'],
                    end_year=variable['end_year'],
                )

    print(
        "Simulating data took {:.0f} seconds".format(time.time() - start_time))


if __name__ == '__main__':
    for path in sys.argv[1:]:
        simulate_input_data(recipe_file=path, config_user_file=None)
