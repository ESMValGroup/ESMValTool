# pylint: disable=invalid-name
"""ESMValTool CMORizer for Eppley-VGPM-MODIS data from Oregon State University.

Tier

Source
   http://orca.science.oregonstate.edu/data/1x2/monthly/eppley.r2018.m.chl.m.sst/hdf

Last access
   20190515

Download and processing instructions
   Download and unpack all files under a single directory
   (no subdirectories with years) in ${RAWOBS}/Tier2/Eppley-VGPM-MODIS

Modification history
   20190515-A_lova_to: written.

"""

import logging
import os
import glob
from datetime import datetime as dt
import xarray as xr
import numpy as np

import iris

from .utilities import (constant_metadata, fix_coords, fix_var_metadata,
                        save_variable, set_global_atts)

logger = logging.getLogger(__name__)


def _fix_data(cube, var):
    """Specific data fixes for different variables."""
    logger.info("Fixing data ...")
    with constant_metadata(cube):
        if var == 'intpp':
            cube /= 1000. * 12.01 * 86400.
    return cube


def extract_variable(var_info, raw_info, out_dir, attrs):
    """Extract to all vars."""
    var = var_info.short_name
    cubes = iris.load(raw_info['file'])
    rawvar = raw_info['name']

    for cube in cubes:
        if cube.var_name == rawvar:
            fix_var_metadata(cube, var_info)
            fix_coords(cube)
            _fix_data(cube, var)
            set_global_atts(cube, attrs)
            save_variable(
                cube,
                var,
                out_dir,
                attrs,
                local_keys=['coordinates'],
                unlimited_dimensions=['time'],
            )


def merge_data(in_dir, out_dir, raw_info):
    """Merge all data into a single file."""
    var = raw_info['name']
    datafile = sorted(glob.glob(in_dir + '/' + raw_info['file'] + '*.hdf'))
    for x in datafile:
        ds = xr.open_dataset(x).rename({'fakeDim0': 'lat', 'fakeDim1': 'lon'})
        # create coordinates
        ds = ds.assign_coords(
            time=dt.strptime(ds.attrs['Start Time String'],
                             '%m/%d/%Y %H:%M:%S'))
        ds = ds.expand_dims(dim='time', axis=0)
        dx = 90. / ds.dims['lat']
        ds = ds.assign_coords(
            lat=np.linspace(-90. + dx, 90. - dx, ds.dims['lat']))
        ds.lat.attrs = {'long_name': 'Latitude', 'units': 'degrees_north'}
        ds = ds.assign_coords(
            lon=np.linspace(-180. + dx, 180. - dx, ds.dims['lon']))
        ds.lon.attrs = {'long_name': 'Longitude', 'units': 'degrees_east'}
        # get data
        da = ds[var]
        if x == datafile[0]:
            newda = da
            continue
        newda = xr.concat((newda, da), dim='time')

    # need data flip to match coordinates
    newda.data = np.fliplr(newda.data)

    # save to file
    da = newda.to_dataset(name=var)
    thekeys = {
        'lat': {
            '_FillValue': False
        },
        'lon': {
            '_FillValue': False
        },
        'time': {
            'calendar': 'gregorian'
        },
        var: {
            '_FillValue': da[var].attrs['Hole Value']
        }
    }
    datafile = os.path.join(out_dir, raw_info['file'] + '_merged.nc')
    da.to_netcdf(datafile, encoding=thekeys, unlimited_dims='time')

    logger.info("Merged data written to: %s", datafile)

    return datafile


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

    logger.info("Starting cmorization for Tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)

    # run the cmorization
    for var, vals in cfg['variables'].items():
        var_info = cmor_table.get_variable(vals['mip'], var)
        glob_attrs['mip'] = vals['mip']
        raw_info = {'name': vals['raw'], 'file': vals['file']}

        # merge yearly data and apply binning
        inpfile = merge_data(in_dir, out_dir, raw_info)

        logger.info("CMORizing var %s from file %s", var, inpfile)
        raw_info['file'] = inpfile
        extract_variable(var_info, raw_info, out_dir, glob_attrs)

    # Remove temporary input file
    os.remove(inpfile)
