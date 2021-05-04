"""Utils module for Python cmorizers."""
from pathlib import Path
import datetime
import logging
import os
import re
from contextlib import contextmanager

import iris
import numpy as np
import yaml
from cf_units import Unit
from dask import array as da

from esmvalcore.cmor.table import CMOR_TABLES
from esmvaltool import __version__ as version, __file__ as esmvaltool_file

logger = logging.getLogger(__name__)

REFERENCES_PATH = Path(esmvaltool_file).absolute().parent / 'references'


def add_height2m(cube):
    """Add scalar coordinate 'height' with value of 2m."""
    add_scalar_height_coord(cube, height=2.)


def add_scalar_height_coord(cube, height=2.):
    """Add scalar coordinate 'height' with value of `height`m."""
    logger.debug("Adding height coordinate (%sm)", height)
    height_coord = iris.coords.AuxCoord(
        height,
        var_name='height',
        standard_name='height',
        long_name='height',
        units=Unit('m'),
        attributes={'positive': 'up'})
    cube.add_aux_coord(height_coord, ())


@contextmanager
def constant_metadata(cube):
    """Do cube math without modifying units etc."""
    metadata = cube.metadata
    yield metadata
    cube.metadata = metadata


def convert_timeunits(cube, start_year):
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


def fix_coords(cube, overwrite_time_bounds=True, overwrite_lon_bounds=True,
               overwrite_lat_bounds=True, overwrite_lev_bounds=True,
               overwrite_airpres_bounds=True):
    """
    Fix coordinates to CMOR standards.

    Fixes coordinates eg time to have correct units, bounds etc;
    longitude to be CMOR-compliant 0-360deg; fixes some attributes
    and bounds - the user can avert bounds fixing by using supplied
    arguments; if bounds are None they will be fixed regardless.

    Parameters
    ----------
    cube: iris.cube.Cube
        data cube with coordinates to be fixed.

    overwrite_time_bounds: bool (optional)
        set to False not to overwrite time bounds.

    overwrite_lon_bounds: bool (optional)
        set to False not to overwrite longitude bounds.

    overwrite_lat_bounds: bool (optional)
        set to False not to overwrite latitude bounds.

    overwrite_lev_bounds: bool (optional)
        set to False not to overwrite depth bounds.

    overwrite_airpres_bounds: bool (optional)
        set to False not to overwrite air pressure bounds.

    Returns
    -------
    cube: iris.cube.Cube
        data cube with fixed coordinates.

    """
    # first fix any completely missing coord var names
    _fix_dim_coordnames(cube)
    # fix individual coords
    for cube_coord in cube.coords():
        # fix time
        if cube_coord.var_name == 'time':
            logger.info("Fixing time...")
            cube.coord('time').convert_units(
                Unit('days since 1950-1-1 00:00:00', calendar='gregorian'))
            if overwrite_time_bounds or not cube.coord('time').has_bounds():
                _fix_bounds(cube, cube.coord('time'))

        # fix longitude
        if cube_coord.var_name == 'lon':
            logger.info("Fixing longitude...")
            if cube_coord.ndim == 1:
                if cube_coord.points[0] < 0. and \
                        cube_coord.points[-1] < 181.:
                    cube_coord.points = \
                        cube_coord.points + 180.
                    if overwrite_lon_bounds or not cube_coord.has_bounds():
                        _fix_bounds(cube, cube_coord)
                    cube.attributes['geospatial_lon_min'] = 0.
                    cube.attributes['geospatial_lon_max'] = 360.
                    nlon = len(cube_coord.points)
                    _roll_cube_data(cube, nlon // 2, -1)

        # fix latitude
        if cube_coord.var_name == 'lat':
            logger.info("Fixing latitude...")
            if overwrite_lat_bounds or not cube.coord('latitude').has_bounds():
                _fix_bounds(cube, cube.coord('latitude'))

        # fix depth
        if cube_coord.var_name == 'lev':
            logger.info("Fixing depth...")
            if overwrite_lev_bounds or not cube.coord('depth').has_bounds():
                _fix_bounds(cube, cube.coord('depth'))

        # fix air_pressure
        if cube_coord.var_name == 'air_pressure':
            logger.info("Fixing air pressure...")
            if overwrite_airpres_bounds \
                    or not cube.coord('air_pressure').has_bounds():
                _fix_bounds(cube, cube.coord('air_pressure'))

    # remove CS
    cube.coord('latitude').coord_system = None
    cube.coord('longitude').coord_system = None

    return cube


def fix_var_metadata(cube, var_info):
    """Fix var metadata from CMOR table."""
    if var_info.standard_name == '':
        cube.standard_name = None
    else:
        cube.standard_name = var_info.standard_name
    cube.var_name = var_info.short_name
    cube.long_name = var_info.long_name
    _set_units(cube, var_info.units)
    return cube


def flip_dim_coord(cube, coord_name):
    """Flip (reverse) dimensional coordinate of cube."""
    logger.info("Flipping dimensional coordinate %s...", coord_name)
    coord = cube.coord(coord_name, dim_coords=True)
    coord_idx = cube.coord_dims(coord)[0]
    coord.points = np.flip(coord.points)
    if coord.bounds is not None:
        coord.bounds = np.flip(coord.bounds, axis=0)
    cube.data = da.flip(cube.core_data(), axis=coord_idx)


def read_cmor_config(dataset):
    """Read the associated dataset-specific config file."""
    reg_path = os.path.join(os.path.dirname(__file__), 'cmor_config',
                            dataset + '.yml')
    with open(reg_path, 'r') as file:
        cfg = yaml.safe_load(file)
    cfg['cmor_table'] = \
        CMOR_TABLES[cfg['attributes']['project_id']]
    if 'comment' not in cfg['attributes']:
        cfg['attributes']['comment'] = ''
    return cfg


def save_variable(cube, var, outdir, attrs, **kwargs):
    """Saver function."""
    _fix_dtype(cube)
    # CMOR standard
    try:
        time = cube.coord('time')
    except iris.exceptions.CoordinateNotFoundError:
        time_suffix = None
    else:
        if len(time.points) == 1 and "mon" not in cube.attributes.get('mip'):
            year = str(time.cell(0).point.year)
            time_suffix = '-'.join([year + '01', year + '12'])
        else:
            date1 = str(time.cell(0).point.year) + '%02d' % \
                time.cell(0).point.month
            date2 = str(time.cell(-1).point.year) + '%02d' % \
                time.cell(-1).point.month
            time_suffix = '-'.join([date1, date2])

    name_elements = [
        attrs['project_id'],
        attrs['dataset_id'],
        attrs['modeling_realm'],
        attrs['version'],
        attrs['mip'],
        var,
    ]
    if time_suffix:
        name_elements.append(time_suffix)
    file_name = '_'.join(name_elements) + '.nc'
    file_path = os.path.join(outdir, file_name)
    logger.info('Saving: %s', file_path)
    status = 'lazy' if cube.has_lazy_data() else 'realized'
    logger.info('Cube has %s data [lazy is preferred]', status)
    iris.save(cube, file_path, fill_value=1e20, **kwargs)


def extract_doi_value(tags):
    """Extract doi(s) from a bibtex entry."""
    reference_doi = []
    pattern = r'doi\ = {(.*?)\},'

    if not isinstance(tags, list):
        tags = [tags]

    for tag in tags:
        bibtex_file = REFERENCES_PATH / f'{tag}.bibtex'
        if bibtex_file.is_file():
            reference_entry = bibtex_file.read_text()
            if re.search("doi", reference_entry):
                reference_doi.append(
                    f'doi:{re.search(pattern, reference_entry).group(1)}'
                )
            else:
                reference_doi.append('doi not found')
                logger.warning(
                    'The reference file %s does not have a doi.', bibtex_file
                )
        else:
            reference_doi.append('doi not found')
            logger.warning(
                'The reference file %s does not exist.', bibtex_file
            )
    return ', '.join(reference_doi)


def set_global_atts(cube, attrs):
    """Complete the cmorized file with global metadata."""
    logger.debug("Setting global metadata...")
    attrs = dict(attrs)
    cube.attributes.clear()
    timestamp = datetime.datetime.utcnow()
    timestamp_format = "%Y-%m-%d %H:%M:%S"
    now_time = timestamp.strftime(timestamp_format)

    # Necessary attributes
    try:
        glob_dict = {
            'title': (f"{attrs.pop('dataset_id')} data reformatted for "
                      f"ESMValTool v{version}"),
            'version':
            attrs.pop('version'),
            'tier':
            str(attrs.pop('tier')),
            'source':
            attrs.pop('source'),
            'reference':
            extract_doi_value(attrs.pop('reference')),
            'comment':
            attrs.pop('comment'),
            'user':
            os.environ.get("USER", "unknown user"),
            'host':
            os.environ.get("HOSTNAME", "unknown host"),
            'history':
            f'Created on {now_time}',
            'project_id':
            attrs.pop('project_id'),
        }
    except KeyError:
        raise KeyError(
            "All CMORized datasets need the global attributes 'dataset_id', "
            "'version', 'tier', 'source', 'reference', 'comment' and "
            "'project_id' specified in the configuration file")

    # Additional attributes
    glob_dict.update(attrs)
    cube.attributes = glob_dict


def var_name_constraint(var_name):
    """:mod:`iris.Constraint` using `var_name` of an :mod:`iris.cube.Cube`."""
    return iris.Constraint(cube_func=lambda c: c.var_name == var_name)


def _fix_bounds(cube, dim_coord):
    """Reset and fix all bounds."""
    if len(cube.coord(dim_coord).points) > 1:
        if cube.coord(dim_coord).has_bounds():
            cube.coord(dim_coord).bounds = None
        cube.coord(dim_coord).guess_bounds()

    if cube.coord(dim_coord).has_bounds():
        cube.coord(dim_coord).bounds = da.array(
            cube.coord(dim_coord).core_bounds(), dtype='float64')
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


def _fix_dtype(cube):
    """Fix `dtype` of a cube and its coordinates."""
    if cube.dtype != np.float32:
        logger.info("Converting data type of data from '%s' to 'float32'",
                    cube.dtype)
        cube.data = cube.core_data().astype(np.float32, casting='same_kind')
    for coord in cube.coords():
        if coord.dtype != np.float64:
            logger.info(
                "Converting data type of coordinate points of '%s' from '%s' "
                "to 'float64'", coord.name(), coord.dtype)
            coord.points = coord.core_points().astype(np.float64,
                                                      casting='same_kind')
        if coord.has_bounds() and coord.bounds_dtype != np.float64:
            logger.info(
                "Converting data type of coordinate bounds of '%s' from '%s' "
                "to 'float64'", coord.name(), coord.bounds_dtype)
            coord.bounds = coord.core_bounds().astype(np.float64,
                                                      casting='same_kind')


def _roll_cube_data(cube, shift, axis):
    """Roll a cube data on specified axis."""
    cube.data = da.roll(cube.core_data(), shift, axis=axis)
    return cube


def _set_units(cube, units):
    """Set units in compliance with cf_unit."""
    special = {'psu': 1.e-3, 'Sv': '1e6 m3 s-1'}
    if units in list(special.keys()):
        cube.units = special[units]
    else:
        cube.units = Unit(units)
    return cube
