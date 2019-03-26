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
                        _fix_format,
                        _fix_coords,
                        _read_cmor_config,
                        _save_variable)

logger = logging.getLogger(__name__)

# read in CMOR configuration
CFG = _read_cmor_config('Landschuetzer2016.yml')
PROJ = CFG['proj']
RAW_VAR = CFG['RAW_VAR']
RAW_FILE = CFG['RAW_FILE']


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


def extract_variable(var_info, raw_file, out_dir):
    """Extract to all vars."""
    var = var_info.short_name
    cubes = iris.load(raw_file)
    rawvar = RAW_VAR[var]
    
    for cube in cubes:
        if cube.var_name == rawvar:
            _fix_format(cube, var_info)
            _fix_coords(cube)
            _fix_data(cube, var)
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

    cmor_table = PROJ['cmip_table']

    # run the cmorization
    for var in PROJ['vars']:
        raw_file = os.path.join(in_dir, RAW_FILE[var] + '.nc')
        logger.info("CMORizing var %s in file %s", var, raw_file)
        var_info = cmor_table.get_variable(PROJ['frequency'][var], var)
        extract_variable(var_info, raw_file, out_dir)
