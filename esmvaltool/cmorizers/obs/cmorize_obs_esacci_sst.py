"""ESMValTool CMORizer for WOA data.

Tier
   Tier 2: other freely-available dataset.

Source
   https://data.nodc.noaa.gov/woa/WOA13/DATAv2/

Last access
   20190131

Download and processing instructions
   Download the following files:
     temperature/netcdf/decav81B0/1.00/woa13_decav81B0_t00_01.nc
     salinity/netcdf/decav81B0/1.00/woa13_decav81B0_s00_01.nc
     oxygen/netcdf/all/1.00/woa13_all_o00_01.nc
     nitrate/netcdf/all/1.00/woa13_all_n00_01.nc
     phosphate/netcdf/all/1.00/woa13_all_p00_01.nc
     silicate/netcdf/all/1.00/woa13_all_i00_01.nc

Modification history
   20130328-lovato_tomas: cmorizer revision
   20190131-predoi_valeriu: adapted to v2.
   20190131-demora_lee: written.

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
        years = range(1993, 1996)
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
