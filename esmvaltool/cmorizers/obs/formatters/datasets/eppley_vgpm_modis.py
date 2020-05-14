# pylint: disable=invalid-name
"""ESMValTool CMORizer for Eppley-VGPM-MODIS data from Oregon State University.

Tier

Source
   http://orca.science.oregonstate.edu/data/1x2/monthly/eppley.r2018.m.chl.m.sst/hdf

Last access
   20190515

Download and processing instructions
   Download and unpack all the *.tar files under a single directory
   (no subdirectories with years) in ${RAWOBS}/Tier2/Eppley-VGPM-MODIS

Modification history
   20190515-lovato_tomas: written.

"""

import logging
import os
import glob
from datetime import datetime as dt
import xarray as xr
import numpy as np

import iris

from esmvaltool.cmorizers.obs.utilities import (
    constant_metadata, fix_coords, fix_var_metadata, save_variable,
    set_global_atts)

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
    da = []
    var = raw_info['name']
    filelist = sorted(glob.glob(in_dir + '/' + raw_info['file'] + '*.hdf'))
    for filename in filelist:
        ds = xr.open_dataset(filename, engine='pynio').rename({
            'fakeDim0': 'lat',
            'fakeDim1': 'lon'
        })
        # create coordinates
        ds = ds.assign_coords(
            time=dt.strptime(ds.attrs['Start_Time_String'],
                             '%m/%d/%Y %H:%M:%S'))
        ds = ds.expand_dims(dim='time', axis=0)
        dx = 90. / ds.dims['lat']
        ds = ds.assign_coords(
            lat=np.linspace(-90. + dx, 90. - dx, ds.dims['lat']))
        ds.lat.attrs = {'long_name': 'Latitude', 'units': 'degrees_north'}
        ds = ds.assign_coords(
            lon=np.linspace(-180. + dx, 180. - dx, ds.dims['lon']))
        ds.lon.attrs = {'long_name': 'Longitude', 'units': 'degrees_east'}
        # get current file data
        da.append(ds[var])
    damerge = xr.concat(da, dim='time')

    # need data flip to match coordinates
    damerge.data = np.fliplr(damerge.data)

    # save to file
    ds = damerge.to_dataset(name=var)
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
            '_FillValue': ds[var].attrs['Hole_Value']
        }
    }
    filename = os.path.join(out_dir, raw_info['file'] + '_merged.nc')
    ds.to_netcdf(filename, encoding=thekeys, unlimited_dims='time')

    logger.info("Merged data written to: %s", filename)

    return filename


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

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
