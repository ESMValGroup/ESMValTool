"""ESMValTool CMORizer for GRACE.

Tier
   Tier 3
Source
   https://podaac.jpl.nasa.gov/dataset/TELLUS_GRAC-GRFO_MASCON_CRI_GRID_RL06_V2
Last access
   20201127

Download and processing instructions
 - Go to the above link
 - Click the tab "Data Access"
 - Log in with Earthdata account
 - Download the following files:
     - CLM4.SCALE_FACTOR.JPL.MSCNv02CRI.nc
     - GRCTellus.JPL.200204_202004.GLO.RL06M.MSCNv02CRI.nc
     - LAND_MASK.CRI.nc
 - Download the grace months table  which holds important information
   on data coverage. Save it in the RAWOBSDIR.
      https://podaac-tools.jpl.nasa.gov/drive/files/allData/tellus/L3/docs/GRACE_GRACE-FO_Months_RL06.csv
 - Manually inspect and check the months table


Modification history
   20200630-crezee_bas: written.
   20201127-kazeroni_remi: updated for latest dataset
"""

import logging
import os
from copy import deepcopy
from datetime import datetime

import iris
import numpy as np
import pandas as pd
import xarray as xr
from cf_units import Unit
from dateutil import relativedelta
from esmvalcore.preprocessor import regrid_time

from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)


def _make_monthly_data_contiguous(in_file, out_file, cfg):

    original = xr.open_dataset(in_file)[cfg['variables']['lweGrace']['raw']]

    months_table_file = os.path.join(cfg['in_dir'], cfg['grace_table'])
    # Read CSV file if available
    if os.path.isfile(months_table_file):
        grace_months_table = pd.read_csv(months_table_file)
    else:
        logger.error("CSV file %s does not exist", months_table_file)
    # Construct the time axis
    time_axis = []
    # read the first and last years and months from the csv table
    time_grace = [[], []]  # [start time], [end time]
    time_grace[0].append(grace_months_table['YEAR'].iloc[0])
    time_grace[1].append(grace_months_table['YEAR'].iloc[-1])
    time_grace[0].append(
        datetime.strptime(grace_months_table['MONTH'].iloc[0], '%b').month)
    time_grace[1].append(
        datetime.strptime(grace_months_table['MONTH'].iloc[-1], '%b').month)
    time_grace[0].append(15)
    time_grace[1].append(15)
    start_date = datetime(*time_grace[0])
    end_date = datetime(*time_grace[1])
    while start_date <= end_date:
        time_axis.append(start_date)
        start_date += relativedelta.relativedelta(months=1)

    # Initialize data array with nan
    data = np.ones((len(time_axis), ) + original.shape[1:])
    data[:] = np.nan

    # Now fill the array with grace data
    for nmonth, recindex in enumerate(
            grace_months_table['GRACE/GRACE-FO record index']):
        if not np.isnan(recindex):
            data[nmonth, :, :] = original[int(recindex - 1), :, :].data
    data_array = xr.DataArray(data,
                              coords={
                                  'time': time_axis,
                                  'lat': original.lat,
                                  'lon': original.lon
                              },
                              dims=['time', 'lat', 'lon'])

    dataset = data_array.to_dataset(name=cfg['variables']['lweGrace']['raw'])
    dataset.to_netcdf(out_file)


def _apply_gain_and_land_sea_mask(in_file, out_file, cfg):

    gain_file = os.path.join(cfg['in_dir'], cfg['auxfiles']['scale_factor'])
    lsm_file = os.path.join(cfg['in_dir'], cfg['auxfiles']['land_mask'])

    gain = xr.open_dataset(gain_file)
    lsm = xr.open_dataset(lsm_file)
    data = xr.open_dataset(in_file)

    data = data['lwe_thickness']
    gain = gain['scale_factor']
    lsm = lsm['land_mask']
    data = gain * data
    data = data.transpose('time', 'lat', 'lon')
    data = data.where(lsm)

    # Specify that unit is cm here (will be converted later)
    data.attrs['units'] = 'cm'
    data = data.to_dataset(name='lwe_thickness')
    data.to_netcdf(out_file)


def _cmorize_dataset(in_file, var, cfg, out_dir):
    logger.info("CMORizing variable '%s' from input file '%s'",
                var['short_name'], in_file)
    attributes = deepcopy(cfg['attributes'])
    attributes['mip'] = var['mip']

    cmor_table = cfg['cmor_table']
    definition = cmor_table.get_variable(var['mip'], var['short_name'])

    cube = iris.load_cube(str(in_file),
                          constraint=utils.var_name_constraint(var['raw']))

    # Set correct names
    cube.var_name = definition.short_name
    if definition.standard_name:
        cube.standard_name = definition.standard_name

    cube.long_name = definition.long_name

    # Convert units if required
    cube.convert_units(definition.units)

    # Set global attributes
    utils.set_global_atts(cube, attributes)

    # Setting time right
    cube = regrid_time(cube, 'mon')

    # Set calendar to gregorian instead of proleptic gregorian
    # matplotlib does not correctly format years in proleptic gregorian
    old_unit = cube.coord('time').units
    new_unit = Unit(old_unit.origin, calendar='gregorian')
    cube.coord('time').units = new_unit

    logger.info("Saving CMORized cube for variable %s", cube.var_name)
    utils.save_variable(cube, cube.var_name, out_dir, attributes)

    return in_file


def cmorization(in_dir, out_dir, cfg, cfg_user):
    """Cmorization func call."""
    cfg['work_dir'] = cfg_user['work_dir']
    # Pass on some parameters to cfg file
    cfg['rawobsdir'] = cfg_user['rootpath']['RAWOBS'][0]
    cfg['in_dir'] = in_dir
    # If it doesn't exist, create it
    if not os.path.isdir(cfg['work_dir']):
        logger.info("Creating working directory for resampling: %s",
                    cfg['work_dir'])
        os.mkdir(cfg['work_dir'])

    # run the cmorization
    for short_name, var in cfg['variables'].items():
        var['short_name'] = short_name
        logger.info("Processing var %s", short_name)
        in_file = os.path.join(in_dir, var['file'])
        logger.info("Structure monthly data")
        out_file = os.path.join(cfg['work_dir'],
                                'grace_monthly_data_contiguous.nc')
        _make_monthly_data_contiguous(in_file, out_file, cfg)
        in_file = out_file
        out_file = os.path.join(
            cfg['work_dir'],
            'grace_monthly_data_contiguous_gain_lsm_applied.nc')
        _apply_gain_and_land_sea_mask(in_file, out_file, cfg)
        in_file = out_file
        logger.info("Start CMORization of file %s", in_file)
        _cmorize_dataset(in_file, var, cfg, out_dir)
        logger.info("Finished regridding and CMORizing.")
