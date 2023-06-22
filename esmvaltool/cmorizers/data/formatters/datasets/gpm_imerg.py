"""
ESMValTool CMORizer for GPM-IMERG data from NASA GEODISC DATA ARCHIVE
(30 min.). Data are regridded to a regular 0.5x0.5 degrees grid
and aggregated to hourly values.

Tier
    Tier 2: other freely-available dataset.

Source
    https://www.ncei.noaa.gov/data/international-satellite-cloud-
        climate-project-isccp-h-series-data/access/isccp-basic/hgg

Last access
    20230524

Download and processing instructions
    see download script cmorizers/data/downloaders/datasets/gpm-imerg.py
"""

import copy
import glob
import logging
import os

from datetime import datetime

import iris
import h5py

import dask.array as da

from esmvaltool.cmorizers.data import utilities as utils

from dateutil import relativedelta

from iris.coords import DimCoord

from esmvalcore.preprocessor import regrid

logger = logging.getLogger(__name__)


def _extract_variable(short_name, var, in_files, cfg, out_dir):
    """Extract variable."""
    # load data
    raw_var = var.get('raw', short_name)
    cubes = iris.cube.CubeList([])
    for infile in in_files:
        print(infile)
        tim0 = datetime.now()
        hdf = h5py.File(infile, 'r')

        precipdata = hdf[f'Grid/{raw_var}']
        lat = hdf.get('Grid/lat')
        lon = hdf.get('Grid/lon')
        time = hdf.get('Grid/time')

        latitude = DimCoord(lat[()], var_name='lat',
                            standard_name='latitude', units='degrees')
        longitude = DimCoord(lon[()], var_name='lon',
                             standard_name='longitude', units='degrees')
        time = DimCoord(time[()], var_name='time', standard_name='time',
                        units='seconds since 1970-01-01 00:00:00 UTC')

        # define fill value
        xfilled = da.ma.masked_equal(precipdata, -9999.9)

        # convert units from 'mm hr-1' to 'kg m-2 s-1'
        xfilled /= 3600.0

        tmpcube = iris.cube.Cube(xfilled, dim_coords_and_dims=[
                                 (time, 0), (longitude, 1), (latitude, 2)])

        # flip longitudes and latitudes
        tmpcube.transpose([0, 2, 1])

        # regridding from 0.1x0.1 to 0.5x0.5
        # ----------------------------------
        # Linear regridding is computationally cheap enough to do this for each
        # field instead of doing it for the concatenated array in the very end.
        # Doing so reduces the memory usage as storing the full resolution
        # fields for a whole month (~48*30=1440 fields) requires ~150 GB of
        # memory.
        utils.fix_bounds(tmpcube, tmpcube.coord('longitude'))
        utils.fix_bounds(tmpcube, tmpcube.coord('latitude'))
        # 'area_weighted' regridding results in a vertical line of missing
        # values at lon=180 --> use 'linear' scheme instead, which is good
        # enough and much faster
        cube05 = regrid(tmpcube, target_grid='0.5x0.5', scheme='linear')

        # realize data to be able to close hdf file
        # -----------------------------------------
        # This is needed as keeping ~48*30=1440 files per month open quickly
        # exceeds the maximum number of open files (check with ulimit -n).
        cube05.data
        hdf.close()

        cubes.append(cube05)

        tim1 = datetime.now()
        print(tim1 - tim0)

    cube = cubes.concatenate_cube()

    # aggregate 30-minute data to hourly values
    iris.coord_categorisation.add_hour(cube, cube.coord('time'), name='hour')
    iris.coord_categorisation.add_day_of_year(cube, cube.coord('time'),
                                              name='day_of_year')
    outcube = cube.aggregated_by(['day_of_year', 'hour'], iris.analysis.MEAN)

    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)

    # Fix metadata
    attrs = copy.deepcopy(cfg['attributes'])
    attrs['mip'] = var['mip']
    utils.fix_var_metadata(outcube, cmor_info)
    utils.set_global_atts(outcube, attrs)

    # fix coordinates
    utils.fix_coords(outcube)

    # Save variable
    utils.save_variable(outcube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    if start_date is None:
        start_date = datetime(2006, 1, 1)
    if end_date is None:
        end_date = datetime(2006, 12, 31)
    loop_date = start_date

    while loop_date <= end_date:
        year = loop_date.year
        month = f'{loop_date.month:0>2}'

        for short_name, var in cfg['variables'].items():
            if 'short_name' not in var:
                var['short_name'] = short_name

            # Now get list of files
            filepattern = os.path.join(in_dir + '/' + str(year),
                                       var['file'].format(year=year,
                                                          month=month))
            print(filepattern)
            in_files = glob.glob(filepattern)
            if not in_files:
                logger.warning('Warning: no data found for %d-%s', year, month)
                continue
            _extract_variable(short_name, var, in_files, cfg, out_dir)

        loop_date += relativedelta.relativedelta(months=1)
