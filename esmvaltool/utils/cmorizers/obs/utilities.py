"""Utils module for Python cmorizers."""
import logging
import os

import iris
import numpy as np

import yaml

logger = logging.getLogger(__name__)


def _read_cmor_config(cmor_config):
    """Read cmor configuration in a dict."""
    reg_path = os.path.join(
        os.path.dirname(__file__), 'cmor_config', cmor_config)
    with open(reg_path, 'r') as file:
        cfg = yaml.safe_load(file)
    return cfg


def _convert_timeunits(cube, start_year):
    """Convert time axis from malformed Year 0."""
    # TODO any more weird cases?
    if cube.coord('time').units == 'months since 0000-01-01 00:00:00':
        real_unit = 'months since {}-01-01 00:00:00'.format(str(start_year))
    if cube.coord('time').units == 'days since 0000-01-01 00:00:00':
        real_unit = 'days since {}-01-01 00:00:00'.format(str(start_year))
    cube.coord('time').units = real_unit
    return cube


def _add_metadata(cube, proj):
    """Complete the cmorized file with useful metadata."""
    for att in proj['metadata_attributes']:
        if att not in cube.metadata.attributes:
            cube.metadata.attributes[att] = proj['metadata_attributes'][att]


def _roll_cube_data(data, shift, axis):
    """Roll a cube data on specified axis."""
    data = np.roll(data, shift, axis=axis)
    return data


def _save_variable(cube, var, outdir, yr, proj):
    """Saver function."""
    # CMOR standard
    time_suffix = '-'.join([str(yr) + '01', str(yr) + '12'])
    cmor_prefix = '_'.join([
        'OBS', proj['dataset'], proj['realm'], proj['version'],
        proj['frequency'][var], var
    ])
    file_name = cmor_prefix + '_' + time_suffix + '.nc'
    file_path = os.path.join(outdir, file_name)
    logger.info('Saving: %s', file_path)
    iris.save(cube, file_path)
