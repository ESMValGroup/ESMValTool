"""Simulate test data for `esmvaltool`."""
from __future__ import print_function

import os
import sys
import time

import numpy as np
from dummydata.model2 import Model2
from dummydata.model3 import Model3

from esmvaltool.interface_scripts.data_finder import get_input_filename
from esmvaltool.interface_scripts.yaml_parser import load_namelist
from esmvaltool.main import read_config_file


def write_data_file(short_name, filename, field, start_year, end_year):
    """Write a file containing simulated data."""
    if 'T2M' in field:
        writer = Model2
    elif 'T3M' in field:
        writer = Model3
    else:
        raise NotImplementedError("Cannot create a model from field {}"
                                  .format(field))

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


def simulate_input_data(namelist_file, config_user_file=None):
    """Simulate data for variables defined in namelist"""
    if config_user_file:
        user_config = read_config_file(
            config_file=config_user_file, namelist_name='')
    else:
        user_config = {
            'rootpath': {
                'default': '.',
            },
            'drs': {},
        }

    namelist = load_namelist(namelist_file)

    start_time = time.time()

    for diagnostic in namelist.DIAGNOSTICS.values():
        np.random.seed(0)
        for model in namelist.MODELS + diagnostic.additional_models:
            for variable in diagnostic.variables:
                filename = get_input_filename(
                    variable=variable,
                    model=model,
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
                    start_year=model['start_year'],
                    end_year=model['end_year'],
                )

    print("Simulating data took {:.0f} seconds"
          .format(time.time() - start_time))


if __name__ == '__main__':
    for path in sys.argv[1:]:
        simulate_input_data(namelist_file=path, config_user_file=None)
