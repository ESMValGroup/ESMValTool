"""Utils module for Python cmorizers."""
import datetime
import logging
import os
from contextlib import contextmanager

import iris
import iris.exceptions
import yaml
from cf_units import Unit
from dask import array as da

from esmvaltool import __version__ as version
from esmvaltool.cmor.table import CMOR_TABLES
from esmvaltool._config import get_tag_value

logger = logging.getLogger(__name__)


@contextmanager
def constant_metadata(cube):
    """Do cube math without modifying units etc."""
    metadata = cube.metadata
    yield metadata
    cube.metadata = metadata


def _read_cmor_config(cmor_config):
    """Read the associated dataset-specific config file."""
    reg_path = os.path.join(
        os.path.dirname(__file__), 'cmor_config', cmor_config)
    with open(reg_path, 'r') as file:
        cfg = yaml.safe_load(file)
    cfg['cmor_table'] = \
        CMOR_TABLES[cfg['attributes']['project_id']]
    if 'comment' not in cfg.keys():
        cfg['attributes']['comment'] = ''
    return cfg


def _fix_var_metadata(cube, var_info):
    """Fix var metadata from CMOR table."""
    cube.var_name = var_info.short_name
    cube.standard_name = var_info.standard_name
    cube.long_name = var_info.long_name
    _set_units(cube, var_info.units)
    return cube


def _set_units(cube, units):
    """Set units in compliance with cf_unit."""
    special = {'psu': 1.e-3, 'Sv': '1e6 m3 s-1'}
    if units in list(special.keys()):
        cube.units = special[units]
    else:
        cube.units = Unit(units)
    return cube


def _convert_timeunits(cube, start_year):
    """Convert time axis from malformed Year 0."""
    # TODO any more weird cases?
    if cube.coord('time').units == 'months since 0000-01-01 00:00:00':
        real_unit = 'months since {}-01-01 00:00:00'.format(str(start_year))
    elif cube.coord('time').units == 'days since 0000-01-01 00:00:00':
        real_unit = 'days since {}-01-01 00:00:00'.format(str(start_year))
    elif cube.coord('time').units == 'days since 1950-1-1':
        real_unit = 'days since 1950-1-1 00:00:00'
    else:
        real_unit = cube.coord('time').units
    cube.coord('time').units = real_unit
    return cube


def _fix_dim_coordnames(cube):
    """Perform a check on dim coordinate names."""
    # first check for CMOR standard coord;
    for coord in cube.coords():
        # guess the CMOR-standard x, y, z and t axes if not there
        coord_type = iris.util.guess_coord_axis(coord)
        try:
            coord = cube.coord(axis=coord_type)
        except iris.exceptions.CoordinateNotFoundError:
            logger.warning(
                'Multiple coordinates for axis %s. '
                'This may be an error, specially for regular grids',
                coord_type
            )
            continue

        if coord_type == 'T':
            coord.var_name = 'time'
            coord.attributes = {}

        if coord_type == 'X':
            coord.var_name = 'lon'
            coord.standard_name = 'longitude'
            coord.long_name = 'longitude coordinate'
            coord.units = Unit('degrees')
            coord.attributes = {}

        if coord_type == 'Y':
            coord.var_name = 'lat'
            coord.standard_name = 'latitude'
            coord.long_name = 'latitude coordinate'
            coord.units = Unit('degrees')
            coord.attributes = {}

        if coord_type == 'Z':
            if coord.var_name == 'depth':
                coord.standard_name = 'depth'
                coord.long_name = \
                    'ocean depth coordinate'
                coord.var_name = 'lev'
                coord.attributes['positive'] = 'down'
            if coord.var_name == 'pressure':
                coord.standard_name = 'air_pressure'
                coord.long_name = 'pressure'
                coord.var_name = 'air_pressure'
                coord.attributes['positive'] = 'up'

    return cube


def _fix_bounds(cube, dim_coord):
    """Reset and fix all bounds."""
    if len(cube.coord(dim_coord).points) > 1:
        if not cube.coord(dim_coord).has_bounds():
            cube.coord(dim_coord).bounds = None
            try:
                cube.coord(dim_coord).guess_bounds()
            except iris.exceptions.CoordinateMultiDimError:
                logger.warning(
                    'Bounds could not be guessed for multidimensional %s '
                    'coordinate',
                    dim_coord.standard_name
                )

    if cube.coord(dim_coord).has_bounds():
        cube.coord(dim_coord).bounds = da.array(
            cube.coord(dim_coord).core_bounds(), dtype='float64')
    return cube


def _fix_coords(cube):
    """Fix the time units and values to CMOR standards."""
    # first fix any completely missing coord var names
    _fix_dim_coordnames(cube)
    # fix individual coords
    lon_coord = cube.coord('longitude')
    for cube_coord in cube.coords():
        # fix time
        if cube_coord.var_name == 'time':
            logger.info("Fixing time...")
            cube.coord('time').convert_units(
                Unit('days since 1950-1-1 00:00:00', calendar='gregorian'))
            _fix_bounds(cube, cube.coord('time'))

        # fix longitude
        if cube_coord.var_name == 'lon':
            logger.info("Fixing longitude...")
            if lon_coord.ndim == 1:
                if lon_coord.points[0] < 0. and \
                        lon_coord.points[-1] < 181.:
                    lon_coord.points = \
                        lon_coord.points + 180.
                    _fix_bounds(cube, lon_coord)
                    cube.attributes['geospatial_lon_min'] = 0.
                    cube.attributes['geospatial_lon_max'] = 360.
                    nlon = len(lon_coord.points)
                    _roll_cube_data(cube, int(nlon / 2), -1)

        # fix latitude
        if cube_coord.var_name == 'lat':
            logger.info("Fixing latitude...")
            _fix_bounds(cube, cube.coord('latitude'))

        # fix depth
        if cube_coord.var_name == 'lev':
            logger.info("Fixing depth...")
            _fix_bounds(cube, cube.coord('depth'))

        # fix air_pressure
        if cube_coord.var_name == 'air_pressure':
            logger.info("Fixing air pressure...")
            _fix_bounds(cube, cube.coord('air_pressure'))

    # remove CS
    cube.coord('latitude').coord_system = None
    lon_coord.coord_system = None

    return cube


def _set_global_atts(cube, attrs):
    """Complete the cmorized file with global metadata."""
    logger.info("Set global metadata...")

    if bool(cube.metadata.attributes):
        cube.metadata.attributes.clear()

    timestamp = datetime.datetime.utcnow()
    timestamp_format = "%Y-%m-%d %H:%M:%S"
    now_time = timestamp.strftime(timestamp_format)
    glob_dict = {
        'title':
        attrs['dataset_id'] + ' data reformatted for ESMValTool v' + version,
        'version': attrs['version'],
        'tier': str(attrs['tier']),
        'source': attrs['source'],
        'reference': get_tag_value('references', attrs['reference']),
        'comment': attrs['comment'],
        'user': os.environ["USER"],
        'host': os.environ["HOSTNAME"],
        'history': 'Created on ' + now_time
    }

    for att, value in glob_dict.items():
        cube.metadata.attributes[att] = value


def _roll_cube_data(cube, shift, axis):
    """Roll a cube data on specified axis."""
    cube.data = da.roll(cube.core_data(), shift, axis=axis)
    return cube


def _save_variable(cube, var, outdir, attrs, **kwargs):
    """Saver function."""
    # CMOR standard
    cube_time = cube.coord('time')
    reftime = Unit(cube_time.units.origin, cube_time.units.calendar)
    dates = reftime.num2date(cube_time.points[[0, -1]])
    if len(cube_time.points) == 1:
        year = str(dates[0].year)
        time_suffix = '-'.join([year + '01', year + '12'])
    else:
        date1 = str(dates[0].year) + '%02d' % dates[0].month
        date2 = str(dates[1].year) + '%02d' % dates[1].month
        time_suffix = '-'.join([date1, date2])

    file_name = '_'.join([
        'OBS',
        attrs['dataset_id'],
        attrs['modeling_realm'],
        attrs['version'],
        attrs['mip'],
        var,
        time_suffix,
    ]) + '.nc'
    file_path = os.path.join(outdir, file_name)
    logger.info('Saving: %s', file_path)
    status = 'lazy' if cube.has_lazy_data() else 'realized'
    logger.info('Cube has %s data [lazy is preferred]', status)
    iris.save(cube, file_path, fill_value=1e20, **kwargs)
