"""ESMValTool CMORizer for WOA data.

Tier
   Tier 2: other freely-available dataset.

Source
   WOA18: https://www.ncei.noaa.gov/data/oceans/woa/WOA18/DATA
   WOA13: https://www.ncei.noaa.gov/data/oceans/woa/WOA13/DATAv2

Last access
   WOA18: 20210311

Download and processing instructions
   All handled by the script (download only if local data are missing)

   Alternatively, download the following files:
     temperature/netcdf/decav81B0/1.00/woa18_decav81B0_t00_01.nc
     salinity/netcdf/decav81B0/1.00/woa18_decav81B0_s00_01.nc
     oxygen/netcdf/all/1.00/woa18_all_o00_01.nc
     nitrate/netcdf/all/1.00/woa18_all_n00_01.nc
     phosphate/netcdf/all/1.00/woa18_all_p00_01.nc
     silicate/netcdf/all/1.00/woa18_all_i00_01.nc
   (To get WOA13, replace filenames prefix woa18 with woa13)


Modification history
   20210311-lovato_tomas: handle WOA18/WOA13, raw data download, use OBS6
   20200911-bock_lisa: extend to WOA18
   20190328-lovato_tomas: cmorizer revision
   20190131-predoi_valeriu: adapted to v2.
   20190131-demora_lee: written.
"""

import logging
import os
from warnings import catch_warnings, filterwarnings

import iris
from cf_units import Unit

from esmvaltool.cmorizers.data.utilities import (
    constant_metadata,
    fix_coords,
    fix_var_metadata,
    save_variable,
    set_global_atts,
)

logger = logging.getLogger(__name__)


def _fix_data(cube, var, version):
    """Specific data fixes for different variables."""
    logger.info("Fixing data ...")

    if version == '2018':
        with constant_metadata(cube):
            if var in ['o2', 'po4', 'si', 'no3']:
                cube /= 1000.  # Convert from umol/kg to mol/m^3

    if version == '2013v2':
        with constant_metadata(cube):
            mll_to_mol = ['po4', 'si', 'no3']
            if var in mll_to_mol:
                cube /= 1000.  # Convert from ml/l to mol/m^3
            elif var == 'thetao':
                cube += 273.15  # Convert to Kelvin
            elif var == 'o2':
                cube *= 44.661 / 1000.  # Convert from ml/l to mol/m^3

    return cube


def collect_files(in_dir, var, cfg):
    """Compose input file list and download if missing."""
    file_list = []
    var_dict = cfg['variables'][var]
    in_dir = os.path.join(in_dir, var_dict['name'])

    fname = cfg['attributes']['short_name'].lower(
    ) + '_' + var_dict['file'] + '00_01.nc'
    in_file = os.path.join(in_dir, fname)
    file_list.append(in_file)

    return file_list


def extract_variable(in_files, out_dir, attrs, raw_info, cmor_table):
    """Extract variables and create OBS dataset."""
    var = raw_info['var']
    var_info = cmor_table.get_variable(raw_info['mip'], var)
    rawvar = raw_info['raw_var']
    with catch_warnings():
        filterwarnings(
            action='ignore',
            message='Ignoring netCDF variable .* invalid units .*',
            category=UserWarning,
            module='iris',
        )
        cubes = iris.load(in_files, rawvar)
    iris.util.equalise_attributes(cubes)
    cube = cubes.concatenate_cube()

    # set reference time
    year = raw_info['reference_year']
    cube.coord('time').climatological = False
    cube.coord('time').points = 6.5
    cube.coord('time').units = Unit('months since ' + str(year) +
                                    '-01-01 00:00:00',
                                    calendar='gregorian')

    fix_var_metadata(cube, var_info)
    fix_coords(cube)
    _fix_data(cube, var, attrs['version'])
    set_global_atts(cube, attrs)
    save_variable(cube, var, out_dir, attrs, unlimited_dimensions=['time'])

    # derive ocean surface
    if 'srf_var' in raw_info:
        var_info = cmor_table.get_variable(raw_info['mip'],
                                           raw_info['srf_var'])
        logger.info("Extract surface OBS for %s", raw_info['srf_var'])
        level_constraint = iris.Constraint(cube.var_name, depth=0)
        cube_os = cube.extract(level_constraint)
        fix_var_metadata(cube_os, var_info)
        save_variable(cube_os,
                      raw_info['srf_var'],
                      out_dir,
                      attrs,
                      unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

    # run the cmorization
    for var, vals in cfg['variables'].items():
        in_files = collect_files(in_dir, var, cfg)
        logger.info("CMORizing var %s from input set %s", var, vals['name'])
        raw_info = cfg['variables'][var]
        raw_info.update({
            'var': var,
            'reference_year': cfg['custom']['reference_year'],
        })
        glob_attrs['mip'] = vals['mip']
        extract_variable(in_files, out_dir, glob_attrs, raw_info, cmor_table)
