"""
# #############################################################################
# ESMValTool CMORizer for CDS-SATELLITE-ALBEDO dataset
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
#
#
# Instructions:
#    - Download the raw datafiles from the source above
#    - Run the non-ESMValTool regridding script (cmorize_obs_CDS-SATELLITE-ALBEDO_with_gdal.py)
#
# Caveats
#    - This dataset is at 1-km  resolution originally, therefore a first preprocessing step
#      is taken outside ESMValTool making use of gdalwarp to regrid from 1 km to 0.25 deg resolution
#    - Potentially, this can be implemented within Python, using the gdal bindings. import gdal; gdal.Translate; gdal.Warp etc. 
#
# Modification history
#    20190320-A_crez_ba: adapted from cmorize_obs_WOA.py
#
# #############################################################################
"""

import datetime
import logging
import os
import iris
from cf_units import Unit

from .utilities import (_add_metadata, _convert_timeunits, _fix_coords,
                        _read_cmor_config, _roll_cube_data, _save_variable)

logger = logging.getLogger(__name__)

# Set these parameters
ALL_VARS = ['bhalb', 'dhalb']
yr_start, yr_end = '2003', '2003'
delta_lon = 0.25  # The resolution of data in the longitudinal direction.

# read in CMOR configuration
cfg = _read_cmor_config('CDS-SATELLITE-ALBEDO.yml')
proj = cfg['proj']
timestamp = datetime.datetime.utcnow()
timestamp_format = "%Y-%m-%d %H:%M:%S"
now_time = timestamp.strftime(timestamp_format)
proj['metadata_attributes']['CMORcreated'] = now_time

RAWFILENAMES = cfg['RAWFILENAMES']
FIELDS = cfg['FIELDS']
STANDARD_NAMES = cfg['STANDARD_NAMES']
LONG_NAMES = cfg['LONG_NAMES']


def extract_variable(var, raw_file, out_dir, yr):
    """Extract to all vars."""
    cubes = iris.load(raw_file)
    field = FIELDS[var]
    for cube in cubes:
        print(cube.long_name, field)
        if cube.long_name == field:
            cube.standard_name = STANDARD_NAMES[var]
            cube.long_name = LONG_NAMES[var]
            cube.var_name = var
            _fix_coords(cube)
            # Shift the data in the longitudinal direction (needed because coords are shifted in _fix_coords)
            _roll_cube_data(cube, int(180 / delta_lon), -1)
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
        rawfilename = RAWFILENAMES[var]
        raw_file = os.path.join(in_dir, rawfilename)
        logger.info("CMORizing var %s in file %s", var, raw_file)
        extract_variable(var, raw_file, out_dir, [yr_start, yr_end])
