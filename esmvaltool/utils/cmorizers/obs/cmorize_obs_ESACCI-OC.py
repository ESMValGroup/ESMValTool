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

from .utilities import (constant_metadata, fix_coords, fix_var_metadata,
                        read_cmor_config, save_variable, set_global_atts)

logger = logging.getLogger(__name__)

# read in CMOR configuration
CFG = read_cmor_config('ESACCI-OC.yml')


def _fix_data(cube, var):
    """Specific data fixes for different variables."""
    logger.info("Fixing data ...")
    with constant_metadata(cube):
        if var == 'chl':
            cube *= 1.e-06
    return cube


# pylint: disable=unused-argument
def _fix_auxcoord(cube, field, filename):
    """Add depth auxiliary coordinate for CMIP5 standard."""
    if not cube.coords('depth'):
        depth = 1.
        depth_coord = iris.coords.AuxCoord(
            depth,
            standard_name='depth',
            long_name='depth',
            var_name='depth',
            units='m',
            attributes={'positive': 'down'})
        cube.add_aux_coord(depth_coord)
        cube.coordinates = 'depth'


def extract_variable(var_info, raw_info, out_dir, attrs):
    """Extract to all vars."""
    var = var_info.short_name
    cubes = iris.load(raw_info['file'], callback=_fix_auxcoord)
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


def merge_data(in_dir, out_dir, raw_info, bin):
    """Merge all data into a single (regridded) file."""
    var = raw_info['name']
    do_bin = True if (bin % 2 == 0) & (bin != 0) else False
    comment = ''
    #TODO remove 1997* here below before final publication
    thefiles = glob.glob(in_dir + '/' + raw_info['file'] + '*1997*.nc')
    for x in thefiles:
        ds = xr.open_dataset(x)
        da = ds[var].sel(lat=slice(None, None, -1))
        # remove inconsistent attributes
        for thekeys in [
                'grid_mapping', 'ancillary_variables', 'parameter_vocab_uri'
        ]:
            del da.attrs[thekeys]

        if do_bin:
            da = da.coarsen(lat=bin, boundary='exact').mean()
            da = da.coarsen(lon=bin, boundary='exact').mean()

        if x == thefiles[0]:
            newda = da
            thekeys = [
                'creator_name', 'creator_url', 'license', 'sensor',
                'processing_level'
            ]
            dsmeta = dict((y, ds.attrs[y]) for y in thekeys)
            if do_bin:
                comment = ' '.join([
                    'Data binned using ', "{}".format(bin), 'by',
                    "{}".format(bin), 'cells average'
                ])
                dsmeta['BINNING'] = comment
            continue
        newda = xr.concat((newda, da), dim='time')

    # save to file
    ds = newda.to_dataset(name=var)
    for x, y in dsmeta.items():
        ds.attrs[x] = y
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
    merged_file = os.path.join(out_dir, raw_info['file'] + '_merged.nc')
    ds.to_netcdf(merged_file, encoding=encoding, unlimited_dims='time')

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
