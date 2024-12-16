"""ESMValTool CMORizer for ESACCI-CLOUD data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://public.satproj.klima.dwd.de/data/ESA_Cloud_CCI/CLD_PRODUCTS/v3.0/L3U/AVHRR-PM/

Last access
    20230619

Download and processing instructions
    see downloading script
"""

import copy
import glob
import logging
import os

import cf_units
import iris
import numpy as np
from cf_units import Unit
from dask import array as da
from iris import NameConstraint
from iris.exceptions import ConstraintMismatchError, MergeError
from calendar import monthrange
from datetime import datetime

from esmvaltool.cmorizers.data import utilities as utils
from esmvalcore.preprocessor import regrid, daily_statistics

logger = logging.getLogger(__name__)

def _create_nan_cube(cube, year, month, day):
    """Create cube containing only NaN from existing cube."""
    nan_cube = cube.copy()
    nan_cube.data = da.ma.masked_greater(cube.core_data(), -1e20)

    # Read dataset time unit and calendar from file
    dataset_time_unit = str(nan_cube.coord('time').units)
    dataset_time_calender = nan_cube.coord('time').units.calendar
    # Convert datetime
    newtime = datetime(year=year, month=month, day=day)
    newtime = cf_units.date2num(newtime, dataset_time_unit,
                                dataset_time_calender)
    nan_cube.coord('time').points = float(newtime) + 0.025390625

    return nan_cube

def _extract_variable_daily(short_name, var, cfg, in_dir,
                            out_dir, start_date, end_date):
    """Extract daily variable."""

    fill_cube = None
    glob_attrs = cfg['attributes']
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)

    if not start_date:
        start_date = datetime(glob_attrs['start_year'], 1, 1)
    if not end_date:
        end_date = datetime(glob_attrs['end_year'], 12, 31)

    for year in range(start_date.year, end_date.year + 1):
        for month in range(start_date.month, end_date.month + 1):
            cubes = iris.cube.CubeList()
            num_days = monthrange(year, month)[1]
            for iday in range(1, num_days + 1):

                filelist = glob.glob(os.path.join(in_dir, f'{year}{month:02}' + f'{iday:02}' + var['file']))

                if filelist:

                    for inum, ifile in enumerate(filelist):
                        logger.info("CMORizing file %s", ifile)

                        # load data
                        raw_var = var.get('raw', short_name)

                        for ivar, raw_name in enumerate(raw_var):
                            daily_cube = iris.load_cube(ifile, NameConstraint(var_name=raw_name))

                            # set arbitrary time of day
                            daily_cube.coord('time').points = (daily_cube.coord('time').points
                                                               + (inum + 0.5 * ivar) * 0.1)
                            daily_cube.attributes.clear()
                            daily_cube.coord('time').long_name = 'time'

                            # Fix coordinates
                            daily_cube = utils.fix_coords(daily_cube)
                            #Fix dtype
                            utils.fix_dtype(daily_cube)
                            #Fix metadata
                            utils.fix_var_metadata(daily_cube, cmor_info)

                            if fill_cube is None:
                                fill_cube = daily_cube

                            cubes.append(daily_cube)

                else:

                    logger.info("Fill missing day %s in month %s and year %s", iday, month, year)
                    daily_cube = _create_nan_cube(fill_cube, year, month, iday)
                    cubes.append(daily_cube)

            cube = cubes.concatenate_cube()

            # regridding from 0.05x0.05 to 0.5x0.5
            cube = regrid(cube, target_grid='0.5x0.5', scheme='area_weighted')

            logger.info("Building daily means")
            # Calc daily
            cube = daily_statistics(cube)

            # Fix units
            if short_name == 'clt':
                cube.data = 100 * cube.core_data()
            else:
                if 'raw_units' in var:
                    cube.units = var['raw_units']
                cube.convert_units(cmor_info.units)

            # Fix metadata and  update version information
            attrs = copy.deepcopy(cfg['attributes'])
            attrs['mip'] = var['mip']
            utils.set_global_atts(cube, attrs)

            # Save variable
            utils.save_variable(cube,
                                short_name,
                                out_dir,
                                attrs,
                                unlimited_dimensions=['time'])


# def _extract_variable_monthly(short_name, var, cfg, in_dir, out_dir, start_date, end_date):
#     """Extract monthly variable with improved handling for multiple cubes."""

#     glob_attrs = cfg['attributes']
#     cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)

#     if not start_date:
#         start_date = datetime(glob_attrs['start_year'], 1, 1)
#     if not end_date:
#         end_date = datetime(glob_attrs['end_year'], 12, 31)

#     for year in range(start_date.year, end_date.year + 1):
#         for month in range(1, 13):  # Cover all months
#             # Search for files for the given year and month
#             filelist = glob.glob(os.path.join(in_dir, f"{year}{month:02}" + var['file']))

#             if not filelist:
#                 logger.warning("No monthly file found for %s-%02d", year, month)
#                 continue

#             for ifile in filelist:
#                 logger.info("CMORizing file %s", ifile)
#                 try:
#                     # Attempt to load the cube using a constraint
#                     constraint = iris.Constraint(var_name=short_name)
#                     cube = iris.load_cube(ifile, constraint)

#                     if cube is None:
#                         logger.warning("Cube could not be loaded for file '%s'", ifile)
#                         continue  # Skip this file and move to the next

#                 except (ConstraintMismatchError, MergeError) as e:
#                     logger.warning("Constraint mismatch in file '%s': %s", ifile, e)
#                     cubes = iris.load(ifile)
#                     matching_cubes = [c for c in cubes if c.var_name == short_name]

#                     if not matching_cubes:
#                         logger.error("No cube found with var_name '%s' in file '%s'", short_name, ifile)
#                         continue  # Skip this file

#                     if len(matching_cubes) > 1:
#                         logger.warning(
#                             "Multiple cubes found with var_name '%s' in file '%s'. Using the first one.",
#                             short_name, ifile
#                         )
#                     cube = matching_cubes[0]  # Use the first matching cube

#                 except Exception as e:
#                     logger.error("Unexpected error while loading file '%s': %s", ifile, e)
#                     continue

#                 try:
#                     # Fix coordinates
#                     logger.info("Fixing coordinates for cube '%s'", cube)
#                     cube = utils.fix_coords(cube)

#                     # Regrid to target grid
#                     cube = regrid(cube, target_grid='0.5x0.5', scheme='area_weighted')

#                     # Fix units
#                     if 'raw_units' in var:
#                         cube.units = var['raw_units']
#                     cube.convert_units(cmor_info.units)

#                     # Fix metadata and update global attributes
#                     attrs = copy.deepcopy(cfg['attributes'])
#                     attrs['mip'] = var['mip']
#                     utils.set_global_atts(cube, attrs)

#                     # Save the processed variable
#                     utils.save_variable(
#                         cube,
#                         short_name,
#                         out_dir,
#                         attrs,
#                         unlimited_dimensions=['time']
#                     )
#                 except Exception as e:
#                     logger.error("Error processing cube for file '%s': %s", ifile, e)


def _extract_variable_monthly(short_name, var, cfg, in_dir, out_dir, start_date, end_date):
    """Extract monthly variable with improved handling for multiple cubes."""
    
    glob_attrs = cfg['attributes']
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)

    if not start_date:
        start_date = datetime(glob_attrs['start_year'], 1, 1)
    if not end_date:
        end_date = datetime(glob_attrs['end_year'], 12, 31)

    for year in range(start_date.year, end_date.year + 1):
        for month in range(1, 13):  # Loop through all months (1-12)
            # Construct the file list for the current month
            filelist = glob.glob(os.path.join(in_dir, f"{year}{month:02}" + var['file']))

            if not filelist:
                logger.warning("No monthly file found for %s-%02d", year, month)
                continue

            cubes = iris.cube.CubeList()

            for ifile in filelist:
                logger.info("CMORizing file %s", ifile)

                try:
                    # Extract raw names from the variable dictionary, like in the daily function
                    raw_var = var.get('raw', short_name)
                    for ivar, raw_name in enumerate(raw_var):
                        # Try to load the cube using a constraint based on the raw_name
                        cube = iris.load_cube(ifile, NameConstraint(var_name=raw_name))

                        if cube is None:
                            logger.warning("Cube could not be loaded for file '%s'", ifile)
                            continue  # Skip this file and move to the next

                        # Set an arbitrary time of day for the monthly data (similar to daily)
                        cube.coord('time').points = (cube.coord('time').points + (ivar + 0.5) * 0.1)
                        cube.attributes.clear()
                        cube.coord('time').long_name = 'time'

                        # Fix coordinates
                        cube = utils.fix_coords(cube)
                        
                        # Fix data type
                        utils.fix_dtype(cube)
                        
                        # Fix metadata
                        utils.fix_var_metadata(cube, cmor_info)

                        # Add the cube to the list
                        cubes.append(cube)
                
                except Exception as e:
                    logger.error("Error processing file '%s': %s", ifile, e)

            # After gathering all cubes for the month, concatenate them
            if cubes:
                cube = cubes.concatenate_cube()

                # Regrid the cube to the target grid (e.g., 0.5x0.5)
                cube = regrid(cube, target_grid='0.5x0.5', scheme='area_weighted')

                # Fix units and handle any special cases like 'clt'
                if short_name == 'clt':
                    cube.data = 100 * cube.core_data()  # Example conversion
                else:
                    if 'raw_units' in var:
                        cube.units = var['raw_units']
                    cube.convert_units(cmor_info.units)

                # Set global attributes and fix metadata
                attrs = copy.deepcopy(cfg['attributes'])
                attrs['mip'] = var['mip']
                utils.set_global_atts(cube, attrs)

                # Save the processed variable
                utils.save_variable(
                    cube,
                    short_name,
                    out_dir,
                    attrs,
                    unlimited_dimensions=['time']
                )

            else:
                logger.warning("No valid cubes processed for %s-%02d", year, month)


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """CMORization function call."""
    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", short_name)
        if 'L3U' in var['file']:
            _extract_variable_daily(short_name, var, cfg, in_dir, out_dir,
                                    start_date, end_date)
        elif 'L3C' in var['file']:
            _extract_variable_monthly(short_name, var, cfg, in_dir, out_dir,
                                       start_date, end_date)





