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
from iris import NameConstraint
from calendar import monthrange
from datetime import datetime

from esmvaltool.cmorizers.data import utilities as utils
from esmvalcore.preprocessor import regrid

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


def _extract_variable(short_name, var, cfg, in_dir,
                      out_dir):
    """Extract variable."""

    fill_cube = None
    glob_attrs = cfg['attributes']

    for year in range(glob_attrs['start_year'],
                      glob_attrs['end_year'] + 1):
        for month in range(1,13):
            cubes = iris.cube.CubeList()
            filename = f'{year}{month:02}*' + var['file']
            filelist = glob.glob(os.path.join(in_dir, filename))
            num_days = monthrange(year, month)[1]
            for iday in range(1, num_days+1):

                ifile = glob.glob(os.path.join(in_dir, f'{year}{month:02}' + f'{iday:02}' + var['file']))

                if ifile:
                    logger.info("CMORizing file %s", ifile)

                    # load data
                    raw_var = var.get('raw', short_name)
                    daily_cube = iris.load_cube(ifile, NameConstraint(var_name=raw_var))

                    daily_cube.attributes.clear()
                    daily_cube.coord('time').long_name = 'time'

                    # Fix coordinates
                    daily_cube = utils.fix_coords(daily_cube)
                    #Fix dtype
                    utils.fix_dtype(daily_cube)

                    if fill_cube is None:
                        fill_cube = daily_cube

                else:

                    logger.info("Fill missing day %s in month %s and year %s", iday, month, year)

                    daily_cube = _create_nan_cube(fill_cube, year, month, iday)

                cubes.append(daily_cube)

            cube = cubes.concatenate_cube()

            # regridding from 0.05x0.05 to 0.5x0.5
            cube = regrid(cube, target_grid='0.5x0.5', scheme='area_weighted')

            cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)

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
            utils.fix_var_metadata(cube, cmor_info)
            utils.set_global_atts(cube, attrs)

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
        _extract_variable(short_name, var, cfg, in_dir, out_dir)
