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
from calendar import monthrange
from datetime import datetime

from esmvaltool.cmorizers.data import utilities as utils
from esmvalcore.preprocessor import regrid, daily_statistics

logger = logging.getLogger(__name__)


def _create_nan_cube(cube, year, month, day):
    """Create cube containing only NaN values from an existing cube."""
    nan_cube = cube.copy()
    nan_cube.data = da.ma.masked_greater(cube.core_data(), -1e20)

    # Read dataset time unit and calendar from file
    dataset_time_unit = str(nan_cube.coord('time').units)
    dataset_time_calendar = nan_cube.coord('time').units.calendar

    # Convert datetime to numeric time
    new_time = cf_units.date2num(datetime(year=year, month=month, day=day),
                                 dataset_time_unit, dataset_time_calendar)
    nan_cube.coord('time').points = float(new_time) + 0.025390625

    return nan_cube


def _extract_variable(short_name, var, cfg, in_dir, out_dir, start_date,
                      end_date):
    """Extract and process a variable."""

    fill_cube = None
    glob_attrs = cfg['attributes']
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)

    if not start_date:
        start_date = datetime(glob_attrs['start_year'], 1, 1)
    if not end_date:
        end_date = datetime(glob_attrs['end_year'], 12, 31)

    for year in range(start_date.year, end_date.year + 1):
        for month in range(1, 13):
            cubes = iris.cube.CubeList()
            num_days = monthrange(year, month)[1]

            for day in range(1, num_days + 1):
                file_pattern = f"{year}{month:02}{day:02}" + var['file']
                file_list = glob.glob(os.path.join(in_dir, file_pattern))

                if file_list:
                    for file_index, file_path in enumerate(file_list):
                        logger.info("Processing file: %s", file_path)

                        raw_var = var.get('raw', [short_name])

                        for var_index, raw_name in enumerate(raw_var):
                            daily_cube = iris.load_cube(
                                file_path, NameConstraint(var_name=raw_name))

                            # Adjust time
                            daily_cube.coord('time').points += (
                                file_index + 0.5 * var_index) * 0.1
                            daily_cube.attributes.clear()
                            daily_cube.coord('time').long_name = 'time'

                            # Fix coordinates, dtype, and metadata
                            daily_cube = utils.fix_coords(daily_cube)
                            utils.fix_dtype(daily_cube)
                            utils.fix_var_metadata(daily_cube, cmor_info)

                            if fill_cube is None:
                                fill_cube = daily_cube

                            cubes.append(daily_cube)
                else:
                    logger.info("Filling missing day: %s-%s-%s",
                                year, month, day)
                    daily_cube = _create_nan_cube(fill_cube, year, month, day)
                    cubes.append(daily_cube)

            # Combine cubes and process
            cube = cubes.concatenate_cube()
            cube = regrid(cube, target_grid='0.5x0.5', scheme='area_weighted')

            logger.info("Calculating daily statistics")
            cube = daily_statistics(cube)

            # Adjust units
            if short_name == 'clt':
                cube.data = 100 * cube.core_data()
            else:
                if 'raw_units' in var:
                    cube.units = var['raw_units']
                cube.convert_units(cmor_info.units)

            # Finalize metadata
            attrs = copy.deepcopy(cfg['attributes'])
            attrs['mip'] = var['mip']
            utils.set_global_atts(cube, attrs)

            # Save the variable
            utils.save_variable(cube, short_name, out_dir, attrs,
                                unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date=None,
                end_date=None):
    """Main cmorization function."""
    for short_name, var in cfg['variables'].items():
        logger.info("CMORizing variable: %s", short_name)
        _extract_variable(short_name, var, cfg, in_dir, out_dir,
                          start_date, end_date)
