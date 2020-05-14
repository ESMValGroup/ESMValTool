"""ESMValTool CMORizer for PERSIANN-CDR data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://www.ncei.noaa.gov/data/precipitation-persiann/access/

Last access
    20200422

Download and processing instructions
    Files are available free for download on the indicated site.
    Files are stored as daily nc-files in individual year
        folders.
    Please copy all files in a single directory.

"""

import glob
import logging
import os
import warnings
from pprint import pformat

import dask.array as da
import iris
import iris.coord_categorisation
import numpy as np

from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)


def fix_coord_dimensions(cube):
    """Fix the order of lat and lon."""
    # Swap latitude and longitude coordinates
    flipped_data = np.moveaxis(cube.core_data(), 2, 1)
    coord_spec = [(cube.coord('time'), 0),
                  (cube.coord('latitude'), 1),
                  (cube.coord('longitude'), 2)]
    new_cube = iris.cube.Cube(flipped_data, dim_coords_and_dims=coord_spec)
    new_cube.metadata = cube.metadata

    # Reverse cube along latitude coordinate and fix latitude coordinate
    lat_coord = new_cube.coord('latitude')
    new_cube = iris.util.reverse(new_cube, lat_coord)
    if not lat_coord.is_contiguous():
        lat_coord.bounds = None
        lat_coord.guess_bounds()
    return new_cube


def _get_input_files(in_dir, cfg):
    """Get input files."""
    pattern = os.path.join(in_dir, cfg['input_file_pattern'])
    input_files = glob.glob(pattern)
    logger.debug("Found input files:\n%s", pformat(input_files))
    return input_files


def _preprocess_cubes(cubes):
    """Remove attributes from cubes that prevent concatenation."""
    new_cubes = iris.cube.CubeList()
    for cube in cubes:
        cube.attributes.pop('datetime', None)
        cube.attributes.pop('date_created', None)
        cube.attributes.pop('id', None)
        cube.attributes.pop('time_coverage_start', None)
        cube.attributes.pop('time_coverage_end', None)
        cube.attributes.pop('metadata_link', None)
        cube.attributes.pop('source', None)
        new_cube = fix_coord_dimensions(cube)
        new_cubes.append(new_cube)
    return new_cubes


def _load_cube(input_files):
    """Load single :class:`iris.cube.Cube`."""
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore',
            message='Ignoring netCDF variable',
            category=UserWarning,
            module='iris',
        )
        cubes = iris.load(input_files)
    cubes = _preprocess_cubes(cubes)
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
    cube = _load_cube(input_files)
    cube.var_name = short_name

    # Fix fill values
    cube.data = da.ma.masked_equal(cube.core_data(), -9999.0)

    # Convert data from precipitation_amount to precipitation_flux
    # divide 'mm' values by the number of seconds in one day
    cube.data = cube.core_data() / 86400.0

    # Fix units
    cube.units = 'kg m-2 s-1'

    # Fix coordinates
    utils.fix_coords(cube)

    # Fix metadata
    attrs = cfg['attributes']
    attrs['mip'] = var['mip']
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    # Save variable in daily files
    iris.coord_categorisation.add_year(cube, 'time', name='year')
    years = list(set(cube.coord('year').points))
    for year in years:
        cube_slice = cube.extract(iris.Constraint(year=year))
        utils.save_variable(cube_slice,
                            short_name,
                            out_dir,
                            attrs,
                            unlimited_dimensions=['time'])

    # Save variable in monthly files
    iris.coord_categorisation.add_month_number(cube, 'time', name='month')
    cube = cube.aggregated_by(['month', 'year'], iris.analysis.MEAN)
    attrs['mip'] = "Amon"
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
