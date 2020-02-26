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

from .utilities import (constant_metadata, convert_timeunits, fix_coords,
                        fix_var_metadata, save_variable, set_global_atts)

logger = logging.getLogger(__name__)


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


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

    # run the cmorization
    for var, vals in cfg['variables'].items():
        for yr in cfg['custom']['years']:
            file_suffix = str(yr)[-2:] + '_' + str(yr + 1)[-2:] + '.nc'
            inpfile = os.path.join(in_dir, vals['file'] + file_suffix)
            logger.info("CMORizing var %s from file %s", var, inpfile)
            var_info = cmor_table.get_variable(vals['mip'], var)
            raw_info = {'name': vals['raw'], 'file': inpfile}
            glob_attrs['mip'] = vals['mip']
            extract_variable(var_info, raw_info, out_dir, glob_attrs, yr)
