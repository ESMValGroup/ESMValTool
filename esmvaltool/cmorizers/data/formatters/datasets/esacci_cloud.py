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

import os
import glob
import copy
import logging
from datetime import datetime
from calendar import monthrange
import iris
import cf_units
from iris.cube import CubeList
from iris.util import equalise_attributes
from esmvalcore.cmor.table import CMOR_TABLES
import esmvaltool.cmorizers.obs.utilities as utils
from iris.constraints import NameConstraint

# Set up logging

logger = logging.getLogger(__name__)

def regrid(cube, target_grid='0.5x0.5', scheme='area_weighted'):
    """Dummy regridding function (replace with actual regridding)."""
    logger.info(f"Regridding cube to {target_grid} with {scheme}")
    return cube

def daily_statistics(cube):
    """Dummy daily statistics function (replace with actual statistics computation)."""
    logger.info("Computing daily statistics")
    return cube

def _create_nan_cube(template_cube, year, month, day):
    """Create a cube filled with NaN values as a placeholder for missing data."""
    logger.info(f"Creating NaN cube for missing data on {year}-{month:02}-{day:02}")
    time = datetime(year=year, month=month, day=day)
    newtime = cf_units.date2num(time, template_cube.coord('time').units, template_cube.coord('time').units.calendar)
    new_cube = template_cube.copy(data=None)
    new_cube.coord('time').points = [newtime]
    new_cube.data = float('nan')
    return new_cube

def _extract_variable(short_name, var, cfg, in_dir, out_dir, start_date, end_date, frequency='daily'):
    """Extract variable (supports daily and monthly)."""

    fill_cube = None
    glob_attrs = cfg['attributes']
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)

    if not start_date:
        start_date = datetime(glob_attrs['start_year'], 1, 1)
    if not end_date:
        end_date = datetime(glob_attrs['end_year'], 12, 31)

    # Loop through each year
    for year in range(start_date.year, end_date.year + 1):
        # Loop through each month
        for month in range(start_date.month, end_date.month + 1):
            cubes = CubeList()

            # Handle daily data
            if frequency == 'daily':
                num_days = monthrange(year, month)[1]
                for iday in range(1, num_days + 1):
                    filelist = glob.glob(os.path.join(in_dir, f'{year}{month:02}{iday:02}' + var['file']))

                    if filelist:
                        for inum, ifile in enumerate(filelist):
                            logger.info("CMORizing daily file %s", ifile)

                            raw_var = var.get('raw', short_name)

                            for ivar, raw_name in enumerate(raw_var):
                                daily_cube = iris.load_cube(ifile, NameConstraint(var_name=raw_name))

                                # Adjust time for daily data
                                daily_cube.coord('time').points = (daily_cube.coord('time').points
                                                                   + (inum + 0.5 * ivar) * 0.1)
                                daily_cube.attributes.clear()
                                daily_cube.coord('time').long_name = 'time'

                                daily_cube = utils.fix_coords(daily_cube)
                                utils.fix_dtype(daily_cube)
                                utils.fix_var_metadata(daily_cube, cmor_info)

                                if fill_cube is None:
                                    fill_cube = daily_cube

                                cubes.append(daily_cube)

                    else:
                        logger.info("Fill missing day %s in month %s and year %s", iday, month, year)
                        daily_cube = _create_nan_cube(fill_cube, year, month, iday)
                        cubes.append(daily_cube)

            # Handle monthly data
            elif frequency == 'monthly':
                filelist = glob.glob(os.path.join(in_dir, f'{year}{month:02}' + var['file']))

                if filelist:
                    for inum, ifile in enumerate(filelist):
                        logger.info("CMORizing monthly file %s", ifile)

                        raw_var = var.get('raw', short_name)

                        for ivar, raw_name in enumerate(raw_var):
                            monthly_cube = iris.load_cube(ifile, NameConstraint(var_name=raw_name))

                            # Set time for monthly data (e.g., the middle of the month)
                            newtime = datetime(year=year, month=month, day=15)
                            dataset_time_unit = str(monthly_cube.coord('time').units)
                            dataset_time_calendar = monthly_cube.coord('time').units.calendar
                            newtime = cf_units.date2num(newtime, dataset_time_unit, dataset_time_calendar)
                            monthly_cube.coord('time').points = float(newtime)

                            monthly_cube.attributes.clear()
                            monthly_cube.coord('time').long_name = 'time'

                            monthly_cube = utils.fix_coords(monthly_cube)
                            utils.fix_dtype(monthly_cube)
                            utils.fix_var_metadata(monthly_cube, cmor_info)

                            if fill_cube is None:
                                fill_cube = monthly_cube

                            cubes.append(monthly_cube)

                else:
                    logger.warning("No monthly file found for month %s and year %s", month, year)

            # Combine the cubes into a single one
            if cubes:
                cube = cubes.concatenate_cube()

                # Apply regridding if needed
                cube = regrid(cube, target_grid='0.5x0.5', scheme='area_weighted')

                if frequency == 'daily':
                    logger.info("Building daily means")
                    cube = daily_statistics(cube)

                # Fix units
                if short_name == 'clt':
                    cube.data = 100 * cube.core_data()
                else:
                    if 'raw_units' in var:
                        cube.units = var['raw_units']
                    cube.convert_units(cmor_info.units)

                # Fix metadata and update version information
                attrs = copy.deepcopy(cfg['attributes'])
                attrs['mip'] = var['mip']
                utils.set_global_atts(cube, attrs)

                # Save the variable
                utils.save_variable(cube,
                                    short_name,
                                    out_dir,
                                    attrs,
                                    unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date, frequency='daily'):
    """Main cmorization function."""
    # Run the cmorization for each variable
    for (short_name, var) in cfg['variables'].items():
        logger.info("CMORizing variable '%s' with frequency '%s'", short_name, frequency)
        _extract_variable(short_name, var, cfg, in_dir, out_dir, start_date, end_date, frequency=frequency)