"""ESMValTool CMORizer for ESACCI-OC data.

Tier

Source
   ftp://oceancolour.org/occci-v3.1/geographic/netcdf/monthly/chlor_a/
   user: oc-cci-data
   pass: ELaiWai8ae

Last access
   20190227

Download and processing instructions
   In case of issues with data download, check also the information provided at
       OceanColour webpage https://esa-oceancolour-cci.org/
   Put all files under a single directory (no subdirectories with years)
   in ${RAWOBS}/Tier2/ESACCI-OC

Modification history
   20190227-lovato_tomas: written.

"""

import logging
import os
import urllib

import iris
import xarray as xr

from .utilities import (constant_metadata, fix_coords, fix_var_metadata,
                        save_variable, set_global_atts)

logger = logging.getLogger(__name__)


def _fix_data(cube, var):
    """Specific data fixes for different variables."""
    logger.info("Fixing data ...")
    with constant_metadata(cube):
        if var == 'chl':
            cube *= 1.e-06
    return cube


def _add_depth_coord(cube):
    """Add depth auxiliary coordinate for CMIP5 standard."""
    if not cube.coords('depth'):
        depth = 1.
        depth_coord = iris.coords.AuxCoord(depth,
                                           standard_name='depth',
                                           long_name='depth',
                                           var_name='depth',
                                           units='m',
                                           attributes={'positive': 'down'})
        cube.add_aux_coord(depth_coord)
        cube.coordinates = 'depth'


def collect_files(in_dir, var, cfg):
    """Compose input file list and download if missing."""
    file_list = []
    var_dict = cfg['variables'][var]
    in_dir = os.path.join(in_dir, var_dict['raw'])
    year_start = cfg['custom']['year_start']
    year_end = cfg['custom']['year_end']

    # create list of monthly data
    for year in range(year_start, year_end + 1):
        for mon in range(1, 13):
            fname = var_dict['file'] + "-{:04d}".format(
                year) + "{:02d}".format(
                    mon) + '-' + cfg['attributes']['version'] + '.nc'
            in_file = os.path.join(in_dir, fname)

            # download if missing
            if not os.path.isfile(in_file):
                if not os.path.isdir(in_dir):
                    os.makedirs(in_dir)
                logger.info(
                    'Input file %s is missing. Start FTP download (~200Mb)...',
                    fname)
                url = os.path.join(cfg['attributes']['source'],
                                   var_dict['raw'], str(year), fname)
                urllib.request.urlretrieve(url, filename=in_file)

            file_list.append(in_file)

    return file_list


def extract_variable(var_info, raw_info, out_dir, attrs):
    """Extract to all vars."""
    var = var_info.short_name
    cubes = iris.load(raw_info['file'])
    rawvar = raw_info['name']

    for cube in cubes:
        if cube.var_name == rawvar:
            fix_var_metadata(cube, var_info)
            fix_coords(cube)
            _add_depth_coord(cube)
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


def merge_data(datafile, out_dir, raw_info, bins):
    """Merge all data into a single (regridded) file."""
    var = raw_info['name']
    do_bin = (bins != 0) and (bins % 2 == 0)
    for xval in datafile:
        ds = xr.open_dataset(xval)
        da = ds[var].sel(lat=slice(None, None, -1))
        # remove inconsistent attributes
        for thekeys in [
                'grid_mapping', 'ancillary_variables', 'parameter_vocab_uri'
        ]:
            del da.attrs[thekeys]

        if do_bin:
            da = da.coarsen(lat=bins, boundary='exact').mean()
            da = da.coarsen(lon=bins, boundary='exact').mean()

        if xval == datafile[0]:
            newda = da
            thekeys = [
                'creator_name', 'creator_url', 'license', 'sensor',
                'processing_level'
            ]
            dsmeta = dict((yval, ds.attrs[yval]) for yval in thekeys)
            if do_bin:
                dsmeta['BINNING'] = ' '.join([
                    'Data binned using ', "{}".format(bins), 'by',
                    "{}".format(bins), 'cells average'
                ])
            else:
                dsmeta['BINNING'] = ""
            continue

        newda = xr.concat((newda, da), dim='time')

    # create dataset
    ds = newda.to_dataset(name=var)
    for xval, yval in dsmeta.items():
        ds.attrs[xval] = yval
    ds['lon'].attrs = {'standard_name': 'longitude'}
    ds['lat'].attrs = {'standard_name': 'latitude'}

    # encoding
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
            '_FillValue': 1.e20
        }
    }

    # save to file
    datafile = os.path.join(out_dir, raw_info['file'] + '_merged.nc')
    ds.to_netcdf(datafile, encoding=thekeys, unlimited_dims='time')

    logger.info("Merged data written to: %s", datafile)

    return (datafile, dsmeta['BINNING'])


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

    # run the cmorization
    for var, vals in cfg['variables'].items():
        var_info = cmor_table.get_variable(vals['mip'], var)
        glob_attrs['mip'] = vals['mip']
        raw_info = {'name': vals['raw'], 'file': vals['file']}

        in_files = collect_files(in_dir, var, cfg)

        # merge yearly data and apply binning
        inpfile, addinfo = merge_data(in_files, out_dir, raw_info,
                                      cfg['custom']['bin_size'])

        logger.info("CMORizing var %s from file %s", var, inpfile)
        raw_info['file'] = inpfile
        glob_attrs['comment'] = addinfo + glob_attrs['comment']
        extract_variable(var_info, raw_info, out_dir, glob_attrs)

    # Remove temporary input file
    os.remove(inpfile)
