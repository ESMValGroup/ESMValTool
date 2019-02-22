"""
# #############################################################################
# ESMValTool CMORizer for CSIRO data
# #############################################################################
#
# Tier
#    Tier 2: other freely-available dataset.
#
# Source
#    https://
#
# Last access
#    20190222
#
# Download and processing instructions
#    Download the following files:
#
# Modification history
#    20190131-A_pred_va: written for v2.
#
# #############################################################################
"""

import datetime
import logging
import os

import iris
from cf_units import Unit

from .utilities import (_add_metadata,
                        _convert_timeunits,
                        _fix_coords,
                        _read_cmor_config,
                        _roll_cube_data,
                        _save_variable)

logger = logging.getLogger(__name__)

# used vars
ALL_VARS = ['so']

# all years to be analyzed
ALL_YEARS = [
    1950, 2000
]

# read in CMOR configuration
cfg = _read_cmor_config('CSIRO.yml')
proj = cfg['proj']
timestamp = datetime.datetime.utcnow()
timestamp_format = "%Y-%m-%d %H:%M:%S"
now_time = timestamp.strftime(timestamp_format)
proj['metadata_attributes']['CMORcreated'] = now_time
VAR_TO_FILENAME = cfg['VAR_TO_FILENAME']
FIELDS = cfg['FIELDS']
STANDARD_NAMES = cfg['STANDARD_NAMES']
LONG_NAMES = cfg['LONG_NAMES']


def _fix_metadata(cube, var):
    """Fix all aspects of metadata for different vars."""
    return cube


def _fix_data(cube, var):
    """Specific data fixes for different variables."""
    return cube


def extract_variable(var, raw_file, out_dir, yr):
    """Extract to all vars."""
    cubes = iris.load(raw_file)
    field = FIELDS[var]
    for cube in cubes:
        if cube.long_name == field:
            cube.standard_name = STANDARD_NAMES[var]
            cube.long_name = LONG_NAMES[var]
            cube.var_name = var
            _convert_timeunits(cube, yr)
            _fix_coords(cube)
            #  _roll_cube_data(cube, 180, -1)
            _fix_data(cube, var)
            _fix_metadata(cube, var)
            _add_metadata(cube, proj)
            _save_variable(cube, var, out_dir, yr, proj)


def cmorization(in_dir, out_dir):
    """Cmorization func call."""
    logger.info("Starting cmorization for CSIRO OBS files: Tier2")
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)

    # run the cmorization
    for var in ALL_VARS:
        if not os.path.exists(out_dir):
            os.path.makedirs(out_dir)
        yr1, yr2 = ALL_YEARS
        file_suffix = str(yr1) + '-' + str(yr2) + '.nc'
        raw_file = os.path.join(in_dir, VAR_TO_FILENAME[var] + file_suffix)
        logger.info("CMORizing var %s in file %s", var, raw_file)
        extract_variable(var, raw_file, out_dir, ALL_YEARS)
