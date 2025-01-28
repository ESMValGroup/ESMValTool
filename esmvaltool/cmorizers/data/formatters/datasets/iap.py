"""ESMValTool CMORizer for IAP data.

Tier
   Tier 2: other freely-available dataset.

Source
   IAPv4.2: http://www.ocean.iap.ac.cn/ftp/cheng/IAPv4.2_IAP_Temperature_gridded_1month_netcdf/Monthly/

Last access: 20250220

Download and processing instructions
   All handled by the script (download only if local data are missing)

   Alternatively, download the following files:
     Temperature_IAPv4.2_gridded_data_1940_1949.zip
     Temperature_IAPv4.2_gridded_data_1950_1959.zip
     Temperature_IAPv4.2_gridded_data_1960_1969.zip
     Temperature_IAPv4.2_gridded_data_1970_1979.zip
     Temperature_IAPv4.2_gridded_data_1980_1989.zip
     Temperature_IAPv4.2_gridded_data_1990_1999.zip
     Temperature_IAPv4.2_gridded_data_2000_2009.zip
     Temperature_IAPv4.2_gridded_data_2010_2019.zip
     Temperature_IAPv4.2_gridded_data_2020_2023.zip
"""

import logging
import os
from warnings import catch_warnings, filterwarnings
from datetime import datetime
from dateutil import relativedelta

import iris
from cf_units import Unit

from esmvaltool.cmorizers.data.utilities import (
    fix_coords,
    fix_var_metadata,
    save_variable,
    set_global_atts,
)

logger = logging.getLogger(__name__)


def collect_files(in_dir, var, cfg, start_date, end_date):
    file_list = []
    #var_dict = cfg['variables'][var]
    #in_dir = os.path.join(in_dir, var_dict['name'])

    if start_date is None:
        start_date = datetime(year=1940, month=1, day=1)
    if end_date is None:
        end_date = datetime(year=2024, month=12, day=31)

    loop_date = start_date

    while loop_date <= end_date:
        fname = f"IAPv4_Temp_monthly_1_6000m_year_{loop_date.year}_month_{loop_date.month:02d}.nc"
        in_file = os.path.join(in_dir, fname)
        file_list.append(in_file)
        loop_date += relativedelta.relativedelta(months=1)

    return file_list


def extract_variable(in_files, out_dir, attrs, raw_info, cmor_table):
    """Extract variables and create OBS dataset."""
    print(raw_info)
    var = raw_info['var']
    var_info = cmor_table.get_variable(raw_info['mip'], var)
    rawvar = raw_info['raw_var']
    #with catch_warnings():
    #    filterwarnings(
    #        action='ignore',
    #        message='Ignoring netCDF variable .* invalid units .*',
    #        category=UserWarning,
    #        module='iris',
    #    )
    #    cubes = iris.load(in_files, rawvar)
    cubes = iris.load(in_files, rawvar)

    print(in_files)

    print(cubes)
    iris.util.equalise_attributes(cubes)
    cube = cubes.concatenate_cube()

    # set reference time
    #year = raw_info['reference_year']
    #cube.coord('time').climatological = False
    #cube.coord('time').points = 6.5
    #cube.coord('time').units = Unit('months since ' + str(year) +
    #                                '-01-01 00:00:00',
    #                                calendar='gregorian')

    fix_var_metadata(cube, var_info)
    fix_coords(cube)
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
        in_files = collect_files(in_dir, var, cfg, start_date, end_date)
        logger.info("CMORizing var %s from input set %s", var, vals['name'])
        raw_info = cfg['variables'][var]
        raw_info.update({
            'var': var,
            #'reference_year': cfg['custom']['reference_year'],
        })
        glob_attrs['mip'] = vals['mip']
        extract_variable(in_files, out_dir, glob_attrs, raw_info, cmor_table)
