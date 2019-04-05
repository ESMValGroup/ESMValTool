# pylint: disable=invalid-name
"""ESMValTool CMORizer for ESACCI-OC data.

Tier

Source
   ftp://oceancolour.org/occci-v3.1/geographic/netcdf/monthly/chlor_a/

Last access
   20190227

Download and processing instructions
   Put all files under a single directory (no subdirectories with years)
   in ${RAWOBS}/Tier2/ESACCI-OC

Modification history
   20190227-A_lova_to: written.

"""

import logging
import os
import glob
import xarray as xr

import iris
import numpy as np

from .utilities import (_set_global_atts, _fix_coords, _fix_var_metadata,
                        _read_cmor_config, _save_variable, constant_metadata)

logger = logging.getLogger(__name__)

# read in CMOR configuration
CFG = _read_cmor_config('ESACCI-OC.yml')


def _fix_data(cube, var):
    """Specific data fixes for different variables."""
    logger.info("Fixing data ...")
    with constant_metadata(cube):
        if var == 'chl':
            cube *= 1.e-06
    return cube


def _fix_attr_cmip5(out_dir, var):
    """Adjust for CMIP5 standard"""
    logger.info("Fix to CMIP5 standard...")
    in_file = glob.glob(out_dir + '/OBS*' + var + '*.nc')[0]
    DS = xr.open_dataset(in_file)
    DS[var].attrs['coordinates'] = 'depth'
    datt = {
        'standard_name': 'depth',
        'long_name': 'depth',
        'units': 'm',
        'axis': 'Z',
        'positive': 'down',
        '_FillValue': False
    }

    DS['depth'] = xr.DataArray(1., name='depth', attrs=datt)
    DS.close()
    DS.to_netcdf(in_file, mode='a')
    return


def extract_variable(var_info, raw_info, out_dir, attrs):
    """Extract to all vars."""
    var = var_info.short_name
    cubes = iris.load(raw_info['file'])
    rawvar = raw_info['name']

    for cube in cubes:
        if cube.var_name == rawvar:
            _fix_var_metadata(cube, var_info)
            _fix_coords(cube)
            _fix_data(cube, var)
            _set_global_atts(cube, attrs)
            _save_variable(
                cube,
                var,
                out_dir,
                attrs,
                local_keys=['coordinates'],
                unlimited_dimensions=['time'],
            )
            _fix_attr_cmip5(out_dir, var)


def da_coarsen(da, bin):
    """ Rebin DataArray with user specified value"""
    data = np.ma.masked_invalid(da.values)
    lat = da.coords['lat'].values
    lon = da.coords['lon'].values
    # rebin
    dd = data.shape
    newda = da[:, :dd[1] // bin, :dd[2] // bin]
    newda.values = data.reshape([dd[0], dd[1] // bin, bin, dd[2] // bin,
                                 bin]).mean(-1).mean(-2)
    newda.coords['lat'].values = lat.reshape([dd[1] // bin, bin]).mean(-1)
    newda.coords['lon'].values = lon.reshape([dd[2] // bin, bin]).mean(-1)

    return newda


def merge_data(in_dir, out_dir, raw_info, bin):
    """Merge all data into a single (regridded) file"""

    var = raw_info['name']
    file = raw_info['file']
    do_bin = True if (bin % 2 == 0) & (bin != 0) else False
    comment = ''
    merged_file = os.path.join(out_dir, file + '_merged.nc')
    #TODO remove 1997* here below before final publication
    in_files = glob.glob(in_dir + '/' + file + '*1997*.nc')
    for ff in in_files:
        DS = xr.open_dataset(ff)
        da = DS[var].sel(lat=slice(None, None, -1))
        # remove inconsitent attributes
        for delkey in [
                'grid_mapping', 'ancillary_variables', 'parameter_vocab_uri'
        ]:
            del da.attrs[delkey]
        #TODO test with xarray coarsen at v0.12+
        if (do_bin):
            da = da_coarsen(da, bin)
        if ff == in_files[0]:
            newda = da
            selkey = ['creator_name','creator_url',\
                    'license','sensor','processing_level']
            DSmeta = dict((k, DS.attrs[k]) for k in selkey)
            if (do_bin):
                comment = ' '.join([
                    'Data binned using ',
                    str(bin), 'by',
                    str(bin), 'cells average'
                ])
                DSmeta['BINNING'] = comment
            continue
        newda = xr.concat((newda, da), dim='time')

    # save to file
    DS = newda.to_dataset()
    for x, y in DSmeta.items():
        DS.attrs[x] = y
    #TODO test encoding with xarray at v0.12+
    encoding = {
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
            '_FillValue': 1.e20
        }
    }
    DS.to_netcdf(merged_file, encoding=encoding, unlimited_dims='time')

    logger.info("Merged data written to: %s", merged_file)

    return (merged_file, comment)


def cmorization(in_dir, out_dir):
    """Cmorization func call."""
    cmor_table = CFG['cmor_table']
    glob_attrs = CFG['attributes']

    logger.info("Starting cmorization for Tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)

    # run the cmorization
    for var, vals in CFG['variables'].items():
        var_info = cmor_table.get_variable(vals['mip'], var)
        glob_attrs['mip'] = vals['mip']
        raw_info = {'name': vals['raw'], 'file': vals['file']}

        # merge yearly data and apply binning
        inpfile, addinfo = merge_data(in_dir, out_dir, raw_info,
                                      CFG['custom']['bin_size'])

        logger.info("CMORizing var %s from file %s", var, inpfile)
        raw_info['file'] = inpfile
        glob_attrs['comment'] = addinfo + glob_attrs['comment']
        extract_variable(var_info, raw_info, out_dir, glob_attrs)
