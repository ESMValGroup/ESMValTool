"""ESMValTool CMORizer for ESACCI-WATERVAPOUR data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://wui.cmsaf.eu/safira/action/viewDoiDetails?acronym=COMBI_V001

Last access
    20240125

Download and processing instructions
    To download the data you need to register on the EUMETSAT website.
    After login you could choose at the bottom of the website
    "https://wui.cmsaf.eu/safira/action/viewDoiDetails?acronym=COMBI_V001"
    between the monthly and daily dataset. We chose for both datasets the 
    0.5x0.5 resolution. Press "Order to carts" and then set the time
    period of the dataset and choose " HTTPS/SFTP" as type of data 
    transfer. Afterwards you receive an email with the 'sftp' infos. 
"""

import copy
import logging
import os

import cf_units
import iris
import numpy as np
from iris import NameConstraint
from calendar import monthrange
from datetime import datetime

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _create_nan_cube(cube, year, month, day):

    nan_cube = cube.copy(
        np.ma.masked_all(cube.shape, dtype=cube.dtype))

    # Read dataset time unit and calendar from file
    dataset_time_unit = str(nan_cube.coord('time').units)
    dataset_time_calender = nan_cube.coord('time').units.calendar
    # Convert datetime
    newtime = datetime(year=year, month=month, day=day)
    newtime = cf_units.date2num(newtime, dataset_time_unit,
                                dataset_time_calender)
    nan_cube.coord('time').points = float(newtime)

    return nan_cube


def save_data(var, cfg, cube, out_dir):
    # fix and save data

    # load cmor table
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], var['short_name'])

    if 'raw_units' in var:
        cube.units = var['raw_units']
    cube.convert_units(cmor_info.units)

    # Fix metadata and  update version information
    attrs = copy.deepcopy(cfg['attributes'])
    attrs['mip'] = var['mip']
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(cube,
                        var['short_name'],
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def load_cube(var, ifile):
    # load data
    raw_var = var['raw']
    cube = iris.load_cube(ifile, NameConstraint(var_name=raw_var))

    cube.attributes.clear()
    cube.coord('time').long_name = 'time'

    # Fix coordinates
    cube = utils.fix_coords(cube)
    # Fix dtype
    utils.fix_dtype(cube)

    return cube


def _extract_daily_variable(var, cfg, inpfile, years, months, out_dir):
    """Extract variable."""

    fill_cube = None

    for year in years:
        for month in months:
            cubes = iris.cube.CubeList()
            num_days = monthrange(year, int(month))[1]
            print(num_days)
            days = ["0" + str(da) for da in range(1, 10)] + [str(da) for da in range(10, num_days+1)]
            for day in days:
                ifile = inpfile.format(year=year, month=month, day=day)
                logger.info("CMORizing var %s from file type %s", var, ifile)

                if os.path.exists(ifile):
                    logger.info("CMORizing file %s", ifile)

                    # load data
                    daily_cube = load_cube(var, ifile)

                    if fill_cube is None:
                        fill_cube = daily_cube

                else:
                    logger.info("Fill missing day %s in month %s and year %s", day, month, year)
                    daily_cube = _create_nan_cube(fill_cube, year, month, day)

                cubes.append(daily_cube)

            cube = cubes.concatenate_cube()

            save_data(var, cfg, cube, out_dir)


def _extract_monthly_variable(var, cfg, inpfile, years, months, out_dir):
    """Extract variable."""

    for year in years:
        cubes = iris.cube.CubeList()
        for month in months:
            ifile = inpfile.format(year=year, month=month)
            logger.info("CMORizing var %s from file type %s", var, ifile)

            logger.info("CMORizing file %s", ifile)

            # load data
            monthly_cube = load_cube(var, ifile)

            cubes.append(monthly_cube)

        cube = cubes.concatenate_cube()

        save_data(var, cfg, cube, out_dir)


def _cmorize_variable(short_name, var, cfg, in_dir,
                      out_dir):
    """Extract variable."""

    glob_attrs = cfg['attributes']

    years = range(var['start_year'], var['end_year'] + 1)
    months = ["0" + str(mo) for mo in range(1, 10)] + ["10", "11", "12"]

    inpfile = os.path.join(in_dir, var['filename'])

    if 'day' in var['mip']:
        _extract_daily_variable(var, cfg, inpfile, years, months, out_dir)
    else:
        _extract_monthly_variable(var, cfg, inpfile, years, months, out_dir)
       


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        if 'short_name' not in var:
            var['short_name'] = short_name
        logger.info("CMORizing variable '%s'", short_name)
        _cmorize_variable(short_name, var, cfg, in_dir, out_dir)
