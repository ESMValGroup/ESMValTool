"""
# #############################################################################
# ESMValTool CMORizer for WOA data
# #############################################################################
#
# Tier
#    Tier 2
#
# Source
#    https://cds.climate.copernicus.eu/cdsapp#!/dataset/
#    satellite-albedo
#
# Last access
#    20190311
#
# Caveats: 
#    -Fix coordinates, not good yet! 
#
# Download and processing instructions
#    Download the following files:
#
# Modification history
#    20190320-A_crez_ba: adapted from cmorize_obs_WOA.py
#
# #############################################################################
"""


#TODO fix coordinates
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
ALL_VARS = ['bhalb']

# all years to be analyzed
START_YYYYMM = '200301'
END_YYYYMM = '200312'

# read in CMOR configuration
cfg = _read_cmor_config('CDS-SATELLITE-ALBEDO.yml')
proj = cfg['proj']
timestamp = datetime.datetime.utcnow()
timestamp_format = "%Y-%m-%d %H:%M:%S"
now_time = timestamp.strftime(timestamp_format)
proj['metadata_attributes']['CMORcreated'] = now_time
VAR_TO_FILENAME = cfg['VAR_TO_FILENAME']
FIELDS = cfg['FIELDS']
STANDARD_NAMES = cfg['STANDARD_NAMES']
LONG_NAMES = cfg['LONG_NAMES']


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
            _roll_cube_data(cube, 180, -1)
            _add_metadata(cube, proj)
            _save_variable(cube, var, out_dir, yr, proj)

def cmorization(in_dir, out_dir):
    """Cmorization func call."""
    logger.info("Starting cmorization for CDS-SATELLITE-ALBEDO files: Tier2")
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)

    # run the cmorization
    for var in ALL_VARS:
        if not os.path.exists(out_dir):
            os.path.makedirs(out_dir)
        raw_file = os.path.join(in_dir,'OBS_CDS-SATELLITE-ALBEDO-025deg_sat_L3_' + VAR_TO_FILENAME[var] + '_' + START_YYYYMM + '-' + END_YYYYMM+'.nc')
        logger.info("CMORizing var %s in file %s", var, raw_file)
        extract_variable(var, raw_file, out_dir, '2003')
