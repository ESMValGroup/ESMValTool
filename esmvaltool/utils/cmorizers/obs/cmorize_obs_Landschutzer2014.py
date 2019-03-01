"""
# #############################################################################
# ESMValTool CMORizer for Landschutzer2014 data
# #############################################################################
#
# Tier
#    Tier 2: other freely-available dataset.
#
# Source
#    ftp://ftp.nodc.noaa.gov/nodc/archive/arc0105/0160558/1.1/data/0-data/
#
# Last access
#    20190227
#
# Download and processing instructions
#    Download the following file in ${RAWOBS}/Tier2/Landschutzer2014:
#     pco2_1998-2011_ETH_SOM-FFN_CDIAC_G05.nc
#
#
# Modification history
#    20190227-A_lova_to: written.
#
# #############################################################################
"""

import datetime
import logging
import os

import iris
from cf_units import Unit
import numpy as np

from .utilities import (_add_metadata,
                        _convert_timeunits,
                        _fix_coords,
                        _read_cmor_config,
                        _roll_cube_data,
                        _save_variable)

logger = logging.getLogger(__name__)

# read in CMOR configuration
cfg = _read_cmor_config('Landschutzer2014.yml')
proj = cfg['proj']
timestamp = datetime.datetime.utcnow()
timestamp_format = "%Y-%m-%d %H:%M:%S"
now_time = timestamp.strftime(timestamp_format)
proj['metadata_attributes']['CMORcreated'] = now_time
VAR_TO_CMOR = cfg['VAR_TO_CMOR']
VAR_TO_FILENAME = cfg['VAR_TO_FILENAME']
STANDARD_NAMES = cfg['STANDARD_NAMES']
LONG_NAMES = cfg['LONG_NAMES']
NC4_ZIP = cfg['NC4_ZIP']
ALLVARS = list(VAR_TO_CMOR.keys())


def _fix_metadata(cube, var):
    """Fix all aspects of metadata for different vars."""
    logger.info("Fixing units...")
    if var in ['fgco2', ]:
        cube.units = Unit('kg m-2 s-1')
    if var in ['spco2', 'dpco2', ]:
        cube.units = Unit('Pa')
    return cube


def _fix_data(cube, var):
    """Specific data fixes for different variables."""
    logger.info("Fixing data ...")
    # fix for bad missing value definition
    cube.data = np.ma.masked_values(cube.data, cube.data.fill_value)
    if var in ['fgco2', ]:
        # Assume standard year 365_day
        cube.data = cube.data * -12.01 / 1000. / (86400. * 365.)
        cube.attributes['positive'] = 'down'
    if var in ['spco2', 'dpco2', ]:
        cube.data = cube.data * 101325. / 1.e06
    return cube


def extract_variable(var, raw_file, out_dir, yr):
    """Extract to all vars."""
    cubes = iris.load(raw_file)
    rawvar = VAR_TO_CMOR[var]
    for cube in cubes:
        if cube.var_name == rawvar:
            cube.standard_name = STANDARD_NAMES[var]
            cube.long_name = LONG_NAMES[var]
            cube.var_name = var
            _convert_timeunits(cube, yr)
            _fix_coords(cube)
            _roll_cube_data(cube, 180, -1)
            _fix_data(cube, var)
            _fix_metadata(cube, var)
            _add_metadata(cube, proj)
            yr1 = cube.coord('time').cell(0).point.strftime('%Y')
            yr2 = cube.coord('time').cell(-1).point.strftime('%Y')
            fillvalue = cube.data.fill_value
            _save_variable(cube, var, out_dir,
                           [yr1, yr2], proj, fillvalue=fillvalue,
                           local_keys=['positive'])


def cmorization(in_dir, out_dir):
    """Cmorization func call."""
    logger.info("Starting cmorization for WOA OBS files: Tier2")
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)

    yrdummy = 2000

    # run the cmorization
    for var in ALLVARS:
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        raw_file = os.path.join(in_dir, VAR_TO_FILENAME[var] + '.nc')
        logger.info("CMORizing var %s in file %s", var, raw_file)
        extract_variable(var, raw_file, out_dir, yrdummy)
