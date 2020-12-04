"""ESMValTool CMORizer for ESACCI-SST data.

Tier
   Tier 2: other freely-available dataset.

Source
   https:<insert_here>

Last access
   20202203

Download and processing instructions
   Download the following files:
     Go to http://surftemp.net/regridding/index.html and request regridded data with the default options
     except for the time resolution, which must be monthly. At the time of writing the script the period
     of the data was from 1982 to 2019 inclusive.

Modification history
   20202203-predoi_valeriu: written.
   <insert_revision>
   <insert_approval>

"""

import logging
import os

import iris

from esmvalcore.preprocessor._io import concatenate
from .utilities import (constant_metadata, convert_timeunits, fix_coords,
                        fix_var_metadata, save_variable, set_global_atts)

logger = logging.getLogger(__name__)


def _fix_data(cube, var):
    """Specific data fixes for different variables."""
    logger.info("Fixing data ...")
    with constant_metadata(cube):
        logger.info("Pending")
    return cube


def extract_variable(var_info, raw_info, out_dir, attrs, year):
    """Extract to all vars."""
    var = var_info.short_name
    cubes = iris.load(raw_info['file'])
    rawvar = raw_info['name']

    for cube in cubes:
        if cube.var_name == rawvar:
            fix_var_metadata(cube, var_info)
            convert_timeunits(cube, year)
            fix_coords(cube)
            _fix_data(cube, var)
            set_global_atts(cube, attrs)
            return cube


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

    # run the cmorization
    for var, vals in cfg['variables'].items():
        var_info = cmor_table.get_variable(vals['mip'], var)
        glob_attrs['mip'] = vals['mip']
        raw_info = {'name': vals['raw'], 'file': vals['file']}
        inpfile = os.path.join(in_dir, cfg['filename'])
        logger.info("CMORizing var %s from file type %s", var, inpfile)
        years = range(1982, 2020)
        months = ["0" + str(mo) for mo in range(1, 10)] + ["10", "11", "12"]
        for year in years:
            monthly_cubes = []
            for month in months:
                raw_info['file'] = inpfile.format(year=year, month=month)
                logger.info("CMORizing var %s from file type %s", var, raw_info['file'])
                cube = extract_variable(var_info, raw_info, out_dir, glob_attrs, year)
                monthly_cubes.append(cube)
            yearly_cube = concatenate(monthly_cubes)
            save_variable(
                yearly_cube, var, out_dir, glob_attrs, unlimited_dimensions=['time'])
