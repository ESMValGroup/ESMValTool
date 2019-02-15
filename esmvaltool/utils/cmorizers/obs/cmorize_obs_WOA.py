"""
# #############################################################################
# ESMValTool CMORizer for WOA data
# #############################################################################
#
# Tier
#    Tier 2: other freely-available dataset.
#
# Source
#    https://data.nodc.noaa.gov/woa/WOA13/DATAv2/
#
# Last access
#    20190131
#
# Download and processing instructions
#    Download the following files:
#      temperature/netcdf/decav81B0/1.00/woa13_decav81B0_t00_01.nc
#      salinity/netcdf/decav81B0/1.00/woa13_decav81B0_s00_01.nc
#      oxygen/netcdf/all/1.00/woa13_all_o00_01.nc
#      nitrate/netcdf/all/1.00/woa13_all_n00_01.nc
#      phosphate/netcdf/all/1.00/woa13_all_p00_01.nc
#      silicate/netcdf/all/1.00/woa13_all_i00_01.nc
#
# Modification history
#    20190131-A_pred_va: adapted to v2.
#    20190131-A_demo_le: written.
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
ALL_VARS = ['thetao', 'so', 'no3', 'po4', 'si', 'o2']

# all years to be analyzed
ALL_YEARS = [
    2000,
]

# read in CMOR configuration
cfg = _read_cmor_config('WOA.yml')
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
    mol_m3 = ['si', 'po4', 'no3', 'o2']
    if var in mol_m3:
        cube.units = Unit('mol m-3')
    if var == 'thetao':
        cube.convert_units(Unit('kelvin'))
    if var == 'so':
        cube.units = Unit('Unknown')
    return cube


def _fix_data(cube, var):
    """Specific data fixes for different variables."""
    mll_to_mol = ['po4', 'si', 'no3']
    if var in mll_to_mol:
        cube.data = cube.data / 1000.  # Convert from ml/l to mol/m^3
    if var == 'o2':
        cube.data = cube.data * 44.661 / 1000.  # Convert from ml/l to mol/m^3
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
            _roll_cube_data(cube, 180, -1)
            _fix_data(cube, var)
            _fix_metadata(cube, var)
            _add_metadata(cube, proj)
            _save_variable(cube, var, out_dir, yr, proj)


def cmorization(in_dir, out_dir):
    """Cmorization func call."""
    logger.info("Starting cmorization for WOA OBS files: Tier2")
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)

    # run the cmorization
    for var in ALL_VARS:
        if not os.path.exists(out_dir):
            os.path.makedirs(out_dir)
        for yr in ALL_YEARS:
            file_suffix = str(yr)[-2:] + '_' + str(yr + 1)[-2:] + '.nc'
            raw_file = os.path.join(in_dir, VAR_TO_FILENAME[var] + file_suffix)
            logger.info("CMORizing var %s in file %s", var, raw_file)
            extract_variable(var, raw_file, out_dir, yr)
