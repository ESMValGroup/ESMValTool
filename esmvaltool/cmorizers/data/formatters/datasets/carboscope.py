"""ESMValTool CMORizer for Jena CarboScope data.

Tier
    Tier 2: other freely-available dataset.

Source
    v2020:
        https://www.bgc-jena.mpg.de/CarboScope/?ID=sEXTocNEET_v2020
    v2024:
        https://www.bgc-jena.mpg.de/CarboScope/?file=nbetEXToc_v2024E.flux.nc

Last access
    20240909

Download and processing instructions
    Download the file corresponding to the version you require:
    v2020: 'sEXTocNEET_v2020_daily.nc.gz'.
    v2024: nbetEXToc_v2024E.flux.nc

"""

import gzip
import logging
import os
import shutil
import warnings
from pathlib import Path

import dask.array as da
import iris
from iris import NameConstraint
from cf_units import Unit

from esmvalcore.preprocessor import monthly_statistics

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _clean(filepath):
    """Remove unzipped input file."""
    if os.path.isfile(filepath):
        os.remove(filepath)
        logger.info("Removed cached file %s", filepath)


def _calculate_flux(cube, filename, area_type):
    """Calculate flux (dividing by land/sea area) and mask land/sea."""
    if area_type == 'land':
        area_idx = 0
    elif area_type == 'ocean':
        area_idx = 1
    else:
        raise ValueError(
            f"Expected 'land' or 'ocean' for 'area_type', got '{area_type}'")

    # Get land/sea area fraction
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore',
            message="Ignoring netCDF variable '.*?' invalid units '.*?'",
            category=UserWarning,
            module='iris',
        )
        area_cube = iris.load_cube(str(filename),
                                   NameConstraint(var_name='area'))
    area = area_cube[area_idx].core_data()
    area = da.broadcast_to(area, cube.shape)

    # Mask
    mask = (area == 0.0)
    area = da.ma.masked_array(area, mask=mask)
    cube.data = da.ma.masked_array(cube.core_data(), mask=mask)

    # Calculate flux (sign change since input data and CMOR use different
    # conventions)
    cube.data = -cube.core_data() / area

    # Adapt metadata
    cube.units /= area_cube.units
    cube.attributes['positive'] = 'down'

    return cube


def _load_cube(filename, raw_name):
    """Load single cube."""
    logger.info("Loading '%s' from file %s", raw_name, filename)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore',
            message="Ignoring netCDF variable '.*?' invalid units '.*?'",
            category=UserWarning,
            module='iris',
        )
        cube = iris.load_cube(str(filename),
                              NameConstraint(var_name=raw_name))
    return cube


def _time_operations(cube):
    """Temporal operations on cube."""
    valid_start_year = cube.attributes['yrfi_valid']
    valid_end_year = cube.attributes['yrfe_valid']
    logger.info("Extracting valid years: %d-%d", valid_start_year,
                valid_end_year)
    cube = cube.extract(iris.Constraint(
        time=lambda cell: valid_start_year <= cell.point.year <= valid_end_year
    ))
    logger.info("Calculating monthly means")
    cube = monthly_statistics(cube)
    return cube


def _fix_metadata(cube, short_name, var, cfg):
    """Fix metadata of cube."""
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)

    # Fix cell measures
    if 'fx' not in var['mip']:
        cell_measure = cube.cell_measure('area per grid cell')
        cell_measure.attributes = {}
        cell_measure.standard_name = 'cell_area'
        if short_name == 'nbp':
            cell_measure.var_name = 'areacella'
            cell_measure.long_name = ('Grid-Cell Area for Atmospheric Grid '
                                      'Variables')
        elif short_name == 'fgco2':
            cell_measure.var_name = 'areacello'
            cell_measure.long_name = 'Grid-Cell Area for Ocean Variables'

    # Fix cell methods
    if short_name == 'nbp':
        cube.add_cell_method(iris.coords.CellMethod('mean where land',
                                                    coords='area'))
    elif short_name == 'fgco2':
        cube.add_cell_method(iris.coords.CellMethod('mean where sea',
                                                    coords='area'))
    elif short_name in ('areacella', 'areacello'):
        cube.add_cell_method(iris.coords.CellMethod('sum', coords='area'))
    elif short_name in ('sftlf', 'sftof'):
        cube.add_cell_method(iris.coords.CellMethod('mean', coords='area'))

    # Fix coordinates
    if 'fx' not in var['mip']:
        cube.remove_coord('month_number')
        cube.remove_coord('year')
        cube.coord('time').var_name = 'time'
        cube.coord('time').long_name = 'time'
        cube.coord('time').points = da.around(
            cube.coord('time').core_points(), 3)
        cube.coord('time').bounds = da.around(
            cube.coord('time').core_bounds(), 3)
    utils.fix_coords(cube)
    if 'depth0m' in cmor_info.dimensions:
        depth_coord = iris.coords.AuxCoord(
            0.0,
            var_name='depth',
            standard_name='depth',
            long_name='depth',
            units=Unit('m'),
            attributes={'positive': 'down'},
        )
        cube.add_aux_coord(depth_coord, ())

    # Fix variable metadata
    cube.convert_units(cmor_info.units)
    utils.fix_var_metadata(cube, cmor_info)

    # Fix global attributes
    attrs = cfg['attributes']
    attrs['mip'] = var['mip']
    utils.set_global_atts(cube, attrs)

    return (cube, attrs)


def _extract_variable(short_name, var, cfg, filename, out_dir):
    """Extract variable."""
    raw_name = var.get('raw_name', short_name)
    cube = _load_cube(filename, raw_name)

    # Fix metadata of raw input
    cube.cell_methods = ()
    if 'raw_units' in var:
        cube.units = var['raw_units']
        cube.attributes.pop('invalid_units', None)

    # Temporal operations
    if 'fx' not in var['mip']:
        cube = _time_operations(cube)

    # Variable-specific calculations
    if short_name == 'nbp':
        cube = _calculate_flux(cube, filename, 'land')
    elif short_name == 'fgco2':
        cube = _calculate_flux(cube, filename, 'ocean')
    elif short_name in ('areacella', 'areacello'):
        cube = cube[2]
        cube.remove_coord('rt')
    elif short_name in ('sftlf', 'sftof'):
        total_area = cube[2]
        area_idx = 0 if short_name == 'sftlf' else 1
        cube = cube[area_idx]
        cube.data = cube.core_data() / total_area.core_data()
        cube.units = '1'
        cube.remove_coord('rt')

    # Fix metadata
    (cube, attrs) = _fix_metadata(cube, short_name, var, cfg)

    # Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def _unzip(zip_path, out_dir):
    """Unzip `*.gz` file."""
    logger.info("Found input file '%s'", zip_path)
    unzip_path = Path(out_dir) / zip_path.with_suffix('').name
    logger.info("Unzipping file")
    with gzip.open(zip_path, 'rb') as zip_file:
        with open(unzip_path, 'wb') as unzip_file:
            shutil.copyfileobj(zip_file, unzip_file)
    logger.info("Succefully extracted file to %s", unzip_path)
    return unzip_path


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    in_path = Path(in_dir) / cfg['filename']

    # Unzip file if necessary
    if '.gz' in cfg['filename']:
        data_path = _unzip(in_path, out_dir)
    else:
        data_path = in_path

    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", short_name)
        _extract_variable(short_name, var, cfg, data_path, out_dir)

    # Remove unzipped file if necessary
    if '.gz' in cfg['filename']:
        _clean(data_path)
