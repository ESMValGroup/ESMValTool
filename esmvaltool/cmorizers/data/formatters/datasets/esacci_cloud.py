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

    """Create cube containing only nan from existing cube."""
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


def _extract_variable(short_name, var, cfg, in_dir,
                      out_dir, start_date, end_date):
    """Extract variable."""

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
            #filename = f'{year}{month:02}*' + var['file']
            #filelist = glob.glob(os.path.join(in_dir, filename))
            num_days = monthrange(year, month)[1]
            for iday in range(1, num_days+1):

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
                            print(daily_cube.coords)
                            daily_cube = utils.fix_coords(daily_cube)
                            print(daily_cube.coords)
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

            #cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)

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
            #utils.fix_var_metadata(cube, cmor_info)
            utils.set_global_atts(cube, attrs)

            print(cube.coords)

            # Save variable
            utils.save_variable(cube,
                                short_name,
                                out_dir,
                                attrs,
                                unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", short_name)
        _extract_variable(short_name, var, cfg, in_dir, out_dir,
                          start_date, end_date)
