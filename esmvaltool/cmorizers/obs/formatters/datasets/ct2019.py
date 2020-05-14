"""ESMValTool CMORizer for CT2019 data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://www.esrl.noaa.gov/gmd/ccgg/carbontracker/index.php

Last access
    20200323

Download and processing instructions
    Create a new empty directory ``$RAWOBSPATH/Tier2/CT2019`` (where
    ``$RAWOBSPATH`` is given by your user configuration file) where the raw
    data will be stored. The download of the data is automatically handled by
    this script. If data is already present in this directory, the download is
    skipped (to force a new download delete your old files).

"""

import fnmatch
import glob
import logging
import os
import warnings
from ftplib import FTP
from pprint import pformat

import dask.array as da
import iris

from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)


def _add_aux_coords(cube, input_files, coords_to_add):
    """Add additional auxiliary coordinates to cube."""
    for (coord_name, coord_dims) in coords_to_add.items():
        logger.info("Adding auxiliary coordinate '%s' to '%s'", coord_name,
                    cube.var_name)
        coord_cube = _load_cube(input_files, coord_name)
        utils.fix_coords(coord_cube)
        dim_coords = [c.name() for c in coord_cube.coords(dim_coords=True)]
        if 'boundary' in dim_coords:
            (points, bounds) = _interpolate_center(coord_cube)
            attributes = {
                'comment': 'Coordinate points where estimated as arithmetic '
                           'mean from given coordinate bounds',
            }
        else:
            points = coord_cube.core_data()
            bounds = None
            attributes = {}
        if coord_cube.long_name == 'air_pressure':
            coord_cube.long_name = 'pressure'
            coord_cube.standard_name = 'air_pressure'
            coord_cube.var_name = 'plev'
        aux_coord = iris.coords.AuxCoord(
            points,
            bounds=bounds,
            var_name=coord_cube.var_name,
            standard_name=coord_cube.standard_name,
            long_name=coord_cube.long_name,
            units=coord_cube.units,
            attributes=attributes,
        )
        cube.add_aux_coord(aux_coord, coord_dims)


def _download_files(in_dir, cfg):
    """Download input files using FTP."""
    logger.info("Downloading data from FTP server %s", cfg['ftp_host'])
    logger.info("Looking for files matching %s", os.path.join(
        cfg['data_dir'], cfg['input_file_pattern']))
    input_files = []
    with FTP(cfg['ftp_host']) as ftp_client:
        logger.info(ftp_client.getwelcome())
        ftp_client.login()
        ftp_client.cwd(cfg['data_dir'])
        files_to_download = fnmatch.filter(ftp_client.nlst(),
                                           cfg['input_file_pattern'])
        for filename in files_to_download:
            logger.info("Downloading %s", filename)
            new_path = os.path.join(in_dir, filename)
            with open(new_path, mode='wb') as outfile:
                ftp_client.retrbinary(f'RETR {filename}', outfile.write)
            input_files.append(new_path)
    return input_files


def _get_input_files(in_dir, cfg):
    """Get input files."""
    pattern = os.path.join(in_dir, cfg['input_file_pattern'])
    input_files = glob.glob(pattern)
    if not input_files:
        input_files = _download_files(in_dir, cfg)
    logger.debug("Found input files:\n%s", pformat(input_files))
    return input_files


def _interpolate_center(cube, axis=1):
    """Interpolate center value for grid cells when only boundary is given."""
    indices = [slice(None)] * cube.ndim
    idx_all_but_first = indices.copy()
    idx_all_but_first[axis] = slice(1, None, None)
    idx_all_but_last = indices.copy()
    idx_all_but_last[axis] = slice(None, -1, None)
    data_all_but_first = cube.core_data()[tuple(idx_all_but_first)]
    data_all_but_last = cube.core_data()[tuple(idx_all_but_last)]
    points = (data_all_but_first + data_all_but_last) / 2.0
    bounds = da.stack((data_all_but_last, data_all_but_first), axis=-1)
    return (points, bounds)


def _remove_attributes(cubes):
    """Remove attributes from cubes that prevent concatenation."""
    for cube in cubes:
        cube.attributes.pop('history', None)
        cube.attributes.pop('nco_input_file_list', None)
        cube.attributes.pop('nco_input_file_number', None)
        cube.attributes.pop('nco_openmp_thread_number', None)
        cube.attributes.pop('NCO', None)
        cube.attributes.pop('version', None)


def _load_cube(input_files, constraints):
    """Load single :class:`iris.cube.Cube`."""
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore',
            message='Ignoring netCDF variable',
            category=UserWarning,
            module='iris',
        )
        cubes = iris.load(input_files, constraints)
    _remove_attributes(cubes)
    try:
        cube = cubes.concatenate_cube()
    except iris.exceptions.ConcatenateError:
        if cubes[0].coords('time'):
            raise
        cube = cubes[0]
    return cube


def _extract_variable(short_name, var, cfg, input_files, out_dir):
    """Extract variable."""
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)

    # Extract data
    constraint = var.get('raw_long_name', cmor_info.standard_name)
    cube = _load_cube(input_files, constraint)
    cube.var_name = short_name

    # Add auxiliary variables
    _add_aux_coords(cube, input_files, var.get('add_aux_coords', {}))

    # Variable specific operations
    if short_name == 'co2s':
        cube = cube[:, 0, :, :]
        cube.remove_coord('level')

    # Fix units
    cube.convert_units(cmor_info.units)
    utils.convert_timeunits(cube, 1950)

    # Fix coordinates
    utils.fix_coords(cube)

    # Fix metadata
    attrs = cfg['attributes']
    attrs['mip'] = var['mip']
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    input_files = _get_input_files(in_dir, cfg)

    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", short_name)
        _extract_variable(short_name, var, cfg, input_files, out_dir)
