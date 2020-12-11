"""ESMValTool CMORizer for GPCC data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://opendata.dwd.de/climate_environment/GPCC/html/fulldata-monthly_v2018_doi_download.html
    https://opendata.dwd.de/climate_environment/GPCC/full_data_2018/full_data_monthly_v2018_[025 05 10 25].nc.gz # noqa
Last access
    20200225

Download and processing instructions
    Download the following files:
        full_data_monthly_{version}.nc.gz

"""

import gzip
import logging
import os
import shutil
from warnings import catch_warnings, filterwarnings
from cf_units import Unit
import cftime
import numpy as np
import iris

from . import utilities as utils

logger = logging.getLogger(__name__)


def _clean(filepath):
    """Remove unzipped input file."""
    if os.path.isfile(filepath):
        os.remove(filepath)
        logger.info("Removed cached file %s", filepath)


def _get_centered_timecoord(cube):
    # time points start at the beginning of month at 00:00:00

    # # tried this, doesn't work
    # # the correct point should be the 2nd bound?
    # cube.coord('time').guess_bounds()
    # cube.coord('time').points = cube.coord('time').bounds[:, 1]
    # cube.coord('time').bounds = None

    time = cube.coord('time')
    times = time.units.num2date(time.points)

    # get the bounds right
    starts = [
        cftime.DatetimeNoLeap(c.year, c.month, 1)
        for c in times
    ]
    ends = [
        cftime.DatetimeNoLeap(c.year, c.month +
                                1, 1) if c.month < 12 else
        cftime.DatetimeNoLeap(c.year + 1, 1, 1) for c in times
    ]
    time.bounds = time.units.date2num(np.stack([starts, ends], -1))

    # get the center right
    time.points = [np.mean((t1, t2)) for t1, t2 in time.bounds]


def _extract_variable(short_name, var, version, cfg, filepath, out_dir):
    """Extract variable."""
    raw_var = var.get('raw', short_name)
    with catch_warnings():
        filterwarnings(
            action='ignore',
            message='Ignoring netCDF variable .* invalid units .*',
            category=UserWarning,
            module='iris',
        )
        cube = iris.load_cube(filepath, utils.var_name_constraint(raw_var))

    # Fix units
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    utils._set_units(cube, var.get('raw_units', short_name))
    # fix calendar type
    cal_time = var.get('calendar', short_name)
    origin_time = cube.coord('time').units.origin
    cube.coord('time').units = Unit(origin_time, calendar=cal_time)
    cube.convert_units(cmor_info.units)
    utils.convert_timeunits(cube, 1950)

    _get_centered_timecoord(cube)
    cube.coord('time').convert_units(
        Unit('days since 1950-1-1 00:00:00', calendar='gregorian'))

    # Fix coordinates
    # latitude is flipped
    utils.flip_dim_coord(cube, 'latitude')
    # utils.fix_coords(cube) <- why are bounds for time re-guessed that are not
    # correct??
    utils._fix_dim_coordnames(cube)
    # fix longitude
    cube_coord = cube.coord('longitude')
    if cube_coord.points[0] < 0. and \
            cube_coord.points[-1] < 181.:
        cube_coord.points = \
            cube_coord.points + 180.
        utils._fix_bounds(cube, cube_coord)
        cube.attributes['geospatial_lon_min'] = 0.
        cube.attributes['geospatial_lon_max'] = 360.
        nlon = len(cube_coord.points)
        utils._roll_cube_data(cube, nlon // 2, -1)
    cube_coord = cube.coord('latitude')
    utils._fix_bounds(cube, cube_coord)

    # Fix metadata
    attrs = cfg['attributes']
    attrs['mip'] = var['mip']
    attrs['version'] = version
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def _unzip(short_name, var, version, raw_filepath, out_dir):
    """Unzip `*.gz` file."""
    raw_var = var.get('raw', short_name)
    zip_path = raw_filepath.format(version=version, raw_name=raw_var)
    if not os.path.isfile(zip_path):
        logger.debug("Skipping '%s', file '%s' not found", short_name,
                     zip_path)
        return None
    logger.info("Found input file '%s'", zip_path)
    filename = os.path.basename(zip_path.replace('.gz', ''))
    new_path = os.path.join(out_dir, filename)
    with gzip.open(zip_path, 'rb') as zip_file:
        with open(new_path, 'wb') as new_file:
            shutil.copyfileobj(zip_file, new_file)
    logger.info("Succefully extracted file to %s", new_path)
    return new_path


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    raw_filepath = os.path.join(in_dir, cfg['filename'])

    # Run the cmorization
    for version in cfg['attributes']['version'].values():
        for (short_name, var) in cfg['variables'].items():
            logger.info("CMORizing variable '%s'", short_name)
            filepath = _unzip(short_name, var, version, raw_filepath, out_dir)
            if filepath is None:
                continue
            _extract_variable(short_name, var, version, cfg, filepath, out_dir)
            _clean(filepath)
