# pylint: disable=invalid-name
"""ESMValTool CMORizer script.

# #############################################################################
  Landschuetzer 2016
# #############################################################################

  Tier
     Tier 2: other freely-available dataset.

  Source
     https://www.nodc.noaa.gov/archive/arc0105/0160558/3.3/data/0-data/

  Last access
     20190308

  Download and processing instructions
     Download the following file in ${RAWOBS}/Tier2/Landschuetzer2016:
      spco2_1982-2015_MPI_SOM-FFN_v2016.nc


  Modification history
     20190227-A_lova_to: written.

# #############################################################################

"""

import logging
import os

import iris
from cf_units import Unit
import numpy as np

from .utilities import (_add_metadata,
                        _fix_coords,
                        _read_cmor_config,
                        _save_variable)

logger = logging.getLogger(__name__)

# read in CMOR configuration
CFG = _read_cmor_config('Landschuetzer2016.yml')
PROJ = CFG['proj']
VAR_TO_CMOR = CFG['VAR_TO_CMOR']
VAR_TO_FILENAME = CFG['VAR_TO_FILENAME']
STANDARD_NAMES = CFG['STANDARD_NAMES']
LONG_NAMES = CFG['LONG_NAMES']

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
    if var in ['dpco2', ]:
        cube.data = cube.data * -1.0 * 101325. / 1.e06
    if var in ['spco2', ]:
        cube.data = cube.data * 101325. / 1.e06
    return cube


def extract_variable(var, raw_file, out_dir):
    """Extract to all vars."""
    cubes = iris.load(raw_file)
    rawvar = VAR_TO_CMOR[var]
    for cube in cubes:
        if cube.var_name == rawvar:
            cube.standard_name = STANDARD_NAMES[var]
            cube.long_name = LONG_NAMES[var]
            cube.var_name = var
            _fix_coords(cube)
            _fix_data(cube, var)
            _fix_metadata(cube, var)
            _add_metadata(cube, PROJ)
            yr1 = cube.coord('time').cell(0).point.strftime('%Y')
            yr2 = cube.coord('time').cell(-1).point.strftime('%Y')
            fillvalue = cube.data.fill_value
            _save_variable(cube, var, out_dir,
                           [yr1, yr2], PROJ, fill_value=fillvalue,
                           local_keys=['positive'],
                           unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir):
    """Cmorization func call."""
    logger.info("Starting cmorization for WOA OBS files: Tier2")
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)

    # run the cmorization
    for var in ALLVARS:
        raw_file = os.path.join(in_dir, VAR_TO_FILENAME[var] + '.nc')
        logger.info("CMORizing var %s in file %s", var, raw_file)
        extract_variable(var, raw_file, out_dir)
