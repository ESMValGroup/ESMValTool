"""ESMValTool CMORizer for ESACCI-OC data.

Tier

Source
   ftp://oceancolour.org/occci-v5.0/geographic/netcdf/monthly/chlor_a/
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

import glob
import logging
import os

from datetime import datetime
import iris
import numpy as np
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
    def adjust_dim(dim):
        if dim == 0:
            return dim
        return dim + 1
    if not cube.coords('depth'):
        assert len(cube.shape) == 3
        depth_coord = iris.coords.DimCoord(
            [0.],
            standard_name='depth',
            long_name='depth',
            var_name='lev',
            units='m',
            bounds=[0., 2.5],
            attributes={'positive': 'down'})
        dim_coords = cube.coords(dim_coords=True)
        aux_coords = cube.coords(dim_coords=False)
        dim_coords_and_dims = [
            (coord, adjust_dim(cube.coord_dims(coord)[0]))
            for coord in dim_coords]
        dim_coords_and_dims.append((depth_coord, 1))
        aux_coords_and_dims = [
            (coord, (adjust_dim(d) for d in cube.coord_dims(coord)))
            for coord in aux_coords]
        old_cube = cube
        new_data = cube.core_data()[:, np.newaxis, :, :]
        cube = iris.cube.Cube(
            new_data,
            old_cube.standard_name,
            old_cube.long_name,
            old_cube.var_name,
            old_cube.units,
            old_cube.attributes,
            old_cube.cell_methods,
            dim_coords_and_dims,
            aux_coords_and_dims,
        )
    return cube


def _fix_time(cube, frequency):
    if frequency == "mon":
        time = cube.coord("time")
        units = time.units
        new_dates = units.date2num(np.array([
            [datetime(d.year, d.month, 1),
             datetime(d.year, d.month, 15),
             datetime(d.year+(d.month // 12), (d.month % 12)+1, 1)]
            for d in units.num2date(time.points)
        ]))
        np.savetxt("time.txt", new_dates)
        time.points = new_dates[:, 1]
        time.bounds = new_dates[:, (0, 2)]


def extract_variable(var_info, raw_info, out_dir, attrs):
    """Extract to all vars."""
    var = var_info.short_name
    cubes = iris.load(raw_info['file'])
    rawvar = raw_info['name']

    for cube in cubes:
        if cube.var_name == rawvar:
            fix_var_metadata(cube, var_info)
            _fix_time(cube, var_info.frequency)
            fix_coords(cube, overwrite_time_bounds=False)
            cube = _add_depth_coord(cube)
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


def merge_data(in_dir, out_dir, raw_info, bins):
    """Merge all data into a single (regridded) file."""
    var = raw_info['name']
    do_bin = (bins != 0) and (bins % 2 == 0)
    datafile = sorted(glob.glob(in_dir + '/' + raw_info['file'] + '*.nc'))
    for x in datafile:
        ds = xr.open_dataset(x)
        da = ds[var].sel(lat=slice(None, None, -1))
        # remove inconsistent attributes
        for thekeys in [
                'grid_mapping', 'ancillary_variables', 'parameter_vocab_uri'
        ]:
            da.attrs.pop(thekeys, None)

        if do_bin:
            da = da.coarsen(lat=bins, boundary='exact').mean()
            da = da.coarsen(lon=bins, boundary='exact').mean()

        if x == datafile[0]:
            newda = da
            thekeys = [
                'creator_name', 'creator_url', 'license', 'sensor',
                'processing_level'
            ]
            dsmeta = dict((y, ds.attrs[y]) for y in thekeys)
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
    for x, y in dsmeta.items():
        ds.attrs[x] = y
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

        # merge yearly data and apply binning
        inpfile, addinfo = merge_data(in_dir, out_dir, raw_info,
                                      cfg['custom']['bin_size'])

        logger.info("CMORizing var %s from file %s", var, inpfile)
        raw_info['file'] = inpfile
        glob_attrs['comment'] = addinfo + glob_attrs['comment']
        extract_variable(var_info, raw_info, out_dir, glob_attrs)

    # Remove temporary input file
    os.remove(inpfile)
