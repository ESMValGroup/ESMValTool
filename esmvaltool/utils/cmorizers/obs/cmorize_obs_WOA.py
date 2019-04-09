# pylint: disable=invalid-name
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
   20130328-A_lova_to: cmorizer revision
   20190131-A_pred_va: adapted to v2.
   20190131-A_demo_le: written.

"""

import logging
import os

import iris

from .utilities import (_set_global_atts, _convert_timeunits, _fix_coords,
                        _fix_var_metadata, _read_cmor_config, _save_variable,
                        constant_metadata)

logger = logging.getLogger(__name__)

# read in CMOR configuration
CFG = _read_cmor_config('WOA.yml')


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
            _fix_var_metadata(cube, var_info)
            _convert_timeunits(cube, year)
            _fix_coords(cube)
            _fix_data(cube, var)
            _set_global_atts(cube, attrs)
            _save_variable(
                cube, var, out_dir, attrs, unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir):
    """Cmorization func call."""
    cmor_table = CFG['cmor_table']
    glob_attrs = CFG['attributes']

    logger.info("Starting cmorization for Tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)

    # run the cmorization
    for var, vals in CFG['variables'].items():
        yr = None
        for yr in CFG['custom']['years']:
            file_suffix = str(yr)[-2:] + '_' + str(yr + 1)[-2:] + '.nc'
            inpfile = os.path.join(in_dir, vals['file'] + file_suffix)
            logger.info("CMORizing var %s from file %s", var, inpfile)
            var_info = cmor_table.get_variable(vals['mip'], var)
            raw_info = {'name': vals['raw'], 'file': inpfile}
            glob_attrs['mip'] = vals['mip']
            extract_variable(var_info, raw_info, out_dir, glob_attrs, yr)
