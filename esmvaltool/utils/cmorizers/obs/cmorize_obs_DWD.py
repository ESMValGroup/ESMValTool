# pylint: disable=invalid-name
"""ESMValTool CMORizer for DWD CDC data.

Tier
   Tier 2: other freely-available dataset.

Source
   https://cdc.dwd.de/portal/

Last access
   20190605

Download and processing instructions
   Download any files:
      - tested on temperature/2m/hourly

Modification history
   20190605-A_muel_bn: written.

"""

import logging
import os

import iris
import dask.array as da
import dask.dataframe as dd

from .utilities import (constant_metadata, convert_timeunits, fix_coords,
                        fix_var_metadata, read_cmor_config, save_variable,
                        set_global_atts)

logger = logging.getLogger(__name__)

# read in CMOR configuration
CFG = read_cmor_config('DWD.yml')


def _fix_data(cube, var):
    """Specific data fixes for different variables."""
    logger.info("Fixing data ...")
    with constant_metadata(cube):
        mll_to_mol = ['po4', 'si', 'no3']
        if var in mll_to_mol:
            cube /= 1000.  # Convert from ml/l to mol/m^3
        elif var == 'thetao':
            cube += 273.15  # Convert to Kelvin
        elif var == 'o2':
            cube *= 44.661 / 1000.  # Convert from ml/l to mol/m^3
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
            save_variable(
                cube, var, out_dir, attrs, unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir):
    """Cmorization func call."""
    
    cmor_table = CFG['cmor_table']
    glob_attrs = CFG['attributes']
    
    logger.info("Starting cmorization for Tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    df = dd.read_csv(in_dir + os.sep + 'data_*.csv')
    logger.info("Output will be written to: %s", out_dir)
    
    
    logger.info(cmor_table)
    logger.info(glob_attrs)
    logger.info(df)
#    

#

#    logger.info("Input data from: %s", in_dir)
#    logger.info("Output will be written to: %s", out_dir)
#
#    # run the cmorization
#    for var, vals in CFG['variables'].items():
#        yr = None
#        for yr in CFG['custom']['years']:
#            file_suffix = str(yr)[-2:] + '_' + str(yr + 1)[-2:] + '.nc'
#            inpfile = os.path.join(in_dir, vals['file'] + file_suffix)
#            logger.info("CMORizing var %s from file %s", var, inpfile)
#            var_info = cmor_table.get_variable(vals['mip'], var)
#            raw_info = {'name': vals['raw'], 'file': inpfile}
#            glob_attrs['mip'] = vals['mip']
#            extract_variable(var_info, raw_info, out_dir, glob_attrs, yr)
