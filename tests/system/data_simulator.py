"""Simulate test data for `esmvaltool`."""
from __future__ import print_function

import os
import sys
import time

import numpy as np

from dummydata.model2 import Model2
from dummydata.model3 import Model3
from esmvaltool.interface_scripts import data_finder
from esmvaltool.interface_scripts.yaml_parser import load_namelist
from esmvaltool.main import read_config_file


def write_data_file(var_name, filename, field, start_year, end_year):
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

    kwargs = cfg[var_name] if var_name in cfg else {}

    writer(
        var=var_name,
        oname=filename,
        start_year=start_year,
        stop_year=end_year,
        **kwargs)


def get_input_filename(project_info, model, var):
    """Return the path to input file that is expected by `esmvaltool`.

    This function should match the ESMValTool function
    esmvaltool.interface_scripts.data_finder.get_input_filelist
    """
    project = model['project']

    cfg = data_finder.read_config_file(project)

    # Apply variable-dependent model keys
    for key in 'mip', 'ensemble', 'exp':
        if key in var:
            model[key] = var[key]

    # Set the rootpath
    if project in project_info['GLOBAL']['rootpath']:
        dir1 = project_info['GLOBAL']['rootpath'][project]
    elif 'default' in project_info['GLOBAL']['rootpath']:
        dir1 = project_info['GLOBAL']['rootpath']['default']
    else:
        raise KeyError(
            'default rootpath must be specified in config-user file')

    # Set the drs
    if project in project_info['GLOBAL']['drs']:
        drs = project_info['GLOBAL']['drs'][project]
    else:
        drs = 'default'

    if drs in cfg['input_dir']:
        dir2 = data_finder.replace_tags(cfg['input_dir'][drs], model, var)
    else:
        raise KeyError(
            'drs {} for {} project not specified in config-developer file'
            .format(drs, project))

    dirname = os.path.join(dir1, dir2)

    # Find latest version if required
    if '[latestversion]' in dirname:
        part1, part2 = dirname.split('[latestversion]')
        dirname = os.path.join(part1, 'dummy', part2)

    # Set the filename
    filename = data_finder.replace_tags(cfg['input_file'], model, var)
    if filename.endswith('*'):
        filename = filename.rstrip(
            '*') + "{start_year}01-{end_year}12.nc".format(**model)

    # Full path to files
    return os.path.join(dirname, filename)


def simulate_input_data(namelist_file, config_user_file=None):
    """Simulate data for variables defined in namelist"""
    if config_user_file:
        user_cfg = read_config_file(
            config_file=config_user_file, namelist_name='')
    else:
        user_cfg = {
            'rootpath': {
                'default': '.',
            },
            'drs': {},
        }

    project_info = {'GLOBAL': user_cfg}

    namelist = load_namelist(namelist_file)

    start_time = time.time()

    for diagnostic in namelist.DIAGNOSTICS.values():
        np.random.seed(0)
        for model in namelist.MODELS + diagnostic.additional_models:
            for variable in diagnostic.variables:
                filename = get_input_filename(
                    project_info=project_info, model=model, var=variable)
                dirname = os.path.dirname(filename)
                if not os.path.exists(dirname):
                    print("Creating {}".format(dirname))
                    os.makedirs(dirname)

                print("Writing {}".format(filename))
                write_data_file(
                    var_name=variable['name'],
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
