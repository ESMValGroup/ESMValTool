"""ESMValTool CMORizer for ESACCI-WATERVAPOUR data.

Tier
<<<<<<< HEAD
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
=======
   Tier 3: CDR2 requires registration at EUMETSAT CM SAF.

Source
   https://wui.cmsaf.eu/safira/action/viewDoiDetails?acronym=COMBI_V001

Last access
   20240221

Download and processing instructions
   CDR2 requires registration at EUMETSAT CM SAF, the information on how to
        download the order will be emailed once the order is ready.
   All files need to be in one directory, not in yearly subdirectories.

Modification history
   20240221-malinina_elizaveta: Adjust for daily cmorization and updated
                                filenames, remove CDR1 due to irrelevance.
   20210607-weigel_katja: Fix for monthly time bounds.
   20210408-weigel_katja: written.
"""

import glob
>>>>>>> public/main
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


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize data."""
    glob_attrs = cfg['attributes']

    # run the cmorization
    for var_name, vals in cfg['variables'].items():
        var = vals['short_name']
        var_info = cfg['cmor_table'].get_variable(vals['mip'], var)
        glob_attrs['mip'] = vals['mip']
        raw_info = {'name': vals['raw']}
        inpfile_pattern = os.path.join(in_dir, vals['filename'])
        logger.info("CMORizing var %s from file type %s", var, inpfile_pattern)
        for year in range(vals['start_year'], vals['end_year'] + 1):
            data_cubes = []
            year_inpfile_pattern = inpfile_pattern.format(year=year)
            inpfiles = sorted(glob.glob(year_inpfile_pattern))
            for inpfile in inpfiles:
                raw_info['file'] = inpfile
                logger.info("CMORizing var %s from file type %s", var,
                            raw_info['file'])
                data_cubes.append(
                    extract_variable(var_info, raw_info, glob_attrs, year))
            yearly_cube = concatenate(data_cubes)
            # Fix monthly time bounds
            time = yearly_cube.coord('time')
            time.bounds = get_time_bounds(time, vals['frequency'])
            save_variable(yearly_cube,
                          var,
                          out_dir,
                          glob_attrs,
                          unlimited_dimensions=['time'])
