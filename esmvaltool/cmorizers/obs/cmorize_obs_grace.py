"""ESMValTool CMORizer for GRACE.

Tier
   Tier 3
Source
   https://podaac.jpl.nasa.gov/dataset/TELLUS_GRAC-GRFO_MASCON_CRI_GRID_RL06_V2
Last access
   20200612

Download and processing instructions
 - Go to the above link
 - Click the tab "Data Access"
 - Log in with Earthdata account
 - Downloaded the following files:
     - CLM4.SCALE_FACTOR.JPL.MSCNv02CRI.nc
     - GRCTellus.JPL.200204_202004.GLO.RL06M.MSCNv02CRI.nc
     - LAND_MASK.CRI.nc


Modification history
   20200630-crezee_bas: written.
"""

import logging
import os
import urllib
from copy import deepcopy
from datetime import datetime

import iris
import numpy as np
import pandas as pd
import xarray as xr
from dateutil import relativedelta

from esmvalcore.preprocessor import regrid_time
from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)


def _make_monthly_data_contiguous(in_file, out_file, raw_varname):
    original = xr.open_dataset(in_file)
    original = original[raw_varname]

    response = urllib.request.urlopen(
        'https://podaac-tools.jpl.nasa.gov/drive/files/\
       allData/tellus/L3/docs/GRACE_GRACE-FO_Months_RL06.csv')
    grace_months_table = pd.read_csv(response)

    # Construct the time axis
    time_axis = []
    start_date = datetime(2002, 4, 15)
    end_date = datetime(2019, 12, 15)
    while start_date <= end_date:
        time_axis.append(start_date)
        start_date += relativedelta.relativedelta(months=1)

    # Determine grid shape from original file
    gridshape = list(original.shape[1:])
    # Add the time dimension to the shape
    gridshape.insert(0, len(time_axis))
    gridshape = tuple(gridshape)

    # Initialize data array with nan
    data = np.ones(gridshape)
    data[:] = np.nan

    # Now fill the array with grace data
    for nmonth, recindex in enumerate(
            grace_months_table['GRACE/GRACE-FO record index']):
        if np.isnan(recindex):  # no grace data so skip
            continue
        else:
            data[nmonth, :, :] = original[int(recindex - 1), :, :].data
    data_array = xr.DataArray(data,
                              coords={
                                  'time': time_axis,
                                  'lat': original.lat,
                                  'lon': original.lon
                              },
                              dims=['time', 'lat', 'lon'])

    dataset = data_array.to_dataset(name=raw_varname)
    dataset.to_netcdf(out_file)


def _apply_gain_and_land_sea_mask(in_file, out_file, cfg):

    gain_file = os.path.join(cfg['rawobsdir'], 'Tier3/GRACE/',
                             'CLM4.SCALE_FACTOR.JPL.MSCNv02CRI.nc')
    lsm_file = os.path.join(cfg['rawobsdir'], 'Tier3/GRACE/',
                            'LAND_MASK.CRI.nc')

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

    logger.info("Saving CMORized cube for variable %s", cube.var_name)
    utils.save_variable(cube, cube.var_name, out_dir, attributes)

    return in_file


def cmorization(in_dir, out_dir, cfg, cfg_user):
    """Cmorization func call."""
    cfg['work_dir'] = cfg_user['work_dir']
    # Pass on rawobsdir to cfg file
    cfg['rawobsdir'] = cfg_user['rootpath']['RAWOBS'][0]
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
        _make_monthly_data_contiguous(in_file, out_file, var['raw'])
        in_file = out_file
        out_file = os.path.join(
            cfg['work_dir'],
            'grace_monthly_data_contiguous_gain_lsm_applied.nc')
        _apply_gain_and_land_sea_mask(in_file, out_file, cfg)
        in_file = out_file
        logger.info("Start CMORization of file %s", in_file)
        _cmorize_dataset(in_file, var, cfg, out_dir)
        logger.info("Finished regridding and CMORizing.")
