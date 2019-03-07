"""Utils module for Python cmorizers."""
import logging
import os

import iris
import numpy as np
from cf_units import Unit

import yaml

logger = logging.getLogger(__name__)


# read the associated dataset-specific config file
def _read_cmor_config(cmor_config):
    """Read cmor configuration in a dict."""
    reg_path = os.path.join(
        os.path.dirname(__file__), 'cmor_config', cmor_config)
    with open(reg_path, 'r') as file:
        cfg = yaml.safe_load(file)
    return cfg


def _convert_timeunits(cube, start_year):
    """Convert time axis from malformed Year 0."""
    # TODO any more weird cases?
    if cube.coord('time').units == 'months since 0000-01-01 00:00:00':
        real_unit = 'months since {}-01-01 00:00:00'.format(str(start_year))
    if cube.coord('time').units == 'days since 0000-01-01 00:00:00':
        real_unit = 'days since {}-01-01 00:00:00'.format(str(start_year))
    if cube.coord('time').units == 'days since 1950-1-1':
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

        if coord_type == 'T':
            cube.coord(axis=coord_type).var_name = 'time'
            cube.coord(axis=coord_type).attributes = {}

        if coord_type == 'X':
            cube.coord(axis=coord_type).var_name = 'lon'
            cube.coord(axis=coord_type).standard_name = 'longitude'
            cube.coord(axis=coord_type).long_name = 'longitude coordinate'
            cube.coord(axis=coord_type).units = Unit('degrees')
            cube.coord(axis=coord_type).attributes = {}

        if coord_type == 'Y':
            cube.coord(axis=coord_type).var_name = 'lat'
            cube.coord(axis=coord_type).standard_name = 'latitude'
            cube.coord(axis=coord_type).long_name = 'latitude coordinate'
            cube.coord(axis=coord_type).units = Unit('degrees')
            cube.coord(axis=coord_type).attributes = {}

        if coord_type == 'Z':
            if cube.coord(axis=coord_type).var_name == 'depth':
                cube.coord(axis=coord_type).standard_name = 'depth'
                cube.coord(axis=coord_type).long_name = \
                    'ocean depth coordinate'
                cube.coord(axis=coord_type).var_name = 'lev'
                cube.coord(axis=coord_type).attributes['positive'] = 'down'
            if cube.coord(axis=coord_type).var_name == 'pressure':
                cube.coord(axis=coord_type).standard_name = 'air_pressure'
                cube.coord(axis=coord_type).long_name = 'pressure'
                cube.coord(axis=coord_type).var_name = 'air_pressure'
                cube.coord(axis=coord_type).attributes['positive'] = 'up'

    return cube


def _fix_bounds(cube, dim_coord):
    """Reset and fix all bounds."""
    if len(cube.coord(dim_coord).points) > 1:
        if not cube.coord(dim_coord).has_bounds():
            cube.coord(dim_coord).guess_bounds()
        else:
            cube.coord(dim_coord).bounds = None
            cube.coord(dim_coord).guess_bounds()

    if cube.coord(dim_coord).has_bounds():
        cube.coord(dim_coord).bounds = np.array(
            cube.coord(dim_coord).bounds, dtype='float64')
    return cube


def _fix_coords(cube):
    """Fix the time units and values to CMOR standards."""
    # first fix any completely missing coord var names
    _fix_dim_coordnames(cube)
    # fix individual coords
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
            if cube.coord('longitude').points[0] < 0. and \
                    cube.coord('longitude').points[-1] < 181.:
                cube.coord('longitude').points = \
                    cube.coord('longitude').points + 180.
                _fix_bounds(cube, cube.coord('longitude'))
                cube.attributes['geospatial_lon_min'] = 0.
                cube.attributes['geospatial_lon_max'] = 360.

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
    cube.coord('longitude').coord_system = None

    return cube


def _add_metadata(cube, proj):
    """Complete the cmorized file with useful metadata."""
    logger.info("Add Global metadata...")
    for att in proj['metadata_attributes']:
        if att not in cube.metadata.attributes:
            cube.metadata.attributes[att] = proj['metadata_attributes'][att]


def _roll_cube_data(cube, shift, axis):
    """Roll a cube data on specified axis."""
    cube.data = np.roll(cube.data, shift, axis=axis)
    return cube


def _save_variable(cube, var, outdir, year, proj, **kwargs):
    """Saver function."""
    # CMOR standard
    if not isinstance(year, list):
        time_suffix = '-'.join([str(year) + '01', str(year) + '12'])
    else:
        yr1, yr2 = year
        time_suffix = '-'.join([str(yr1) + '01', str(yr2) + '12'])
    cmor_prefix = '_'.join([
        'OBS', proj['dataset'], proj['realm'], proj['version'],
        proj['frequency'][var], var
    ])
    file_name = cmor_prefix + '_' + time_suffix + '.nc'
    file_path = os.path.join(outdir, file_name)
    logger.info('Saving: %s', file_path)
    iris.save(cube, file_path, **kwargs)
