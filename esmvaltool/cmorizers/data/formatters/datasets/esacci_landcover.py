"""ESMValTool CMORizer for ESACCI-LANDCOVER pft data.

Tier
    Tier 2: other freely-available dataset.

Source
    ftp://anon-ftp.ceda.ac.uk/neodc/esacci/land_cover/data/pft/

Last access
    20240626

Download and processing instructions
    Download the data from:
      pft/v2.0.8/
    Put all files under a single directory (no subdirectories with years).
    in ${RAWOBS}/Tier2/ESACCI-LANDCOVER

"""

import os
import glob
import gc
import logging
from datetime import datetime
import iris
import numpy as np

from ...utilities import (
    fix_coords_esacci,
    fix_var_metadata,
    fix_dtype,
    set_global_atts,
    set_units,
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Enable the new split-attributes handling mode
iris.FUTURE.save_split_attrs = True


def extract_variable(raw_info):
    """Extract the variable from the raw data file."""
    cube_list = iris.load(raw_info['file'], raw_info['name'])
    if not cube_list:
        logger.warning("No cubes found for {raw_info['name']} in file "
                       "{raw_info['file']}")
    else:
        logger.info("Extracted cubes: %s", cube_list)
    return cube_list


def average_block(data, block_size):
    """Average the data within each block of size block_size."""
    shape = data.shape
    reshaped_data = data.reshape(shape[0], shape[1] // block_size,
                                 block_size, shape[2] // block_size,
                                 block_size)
    averaged_data = reshaped_data.mean(axis=(2, 4))
    return averaged_data


def regrid_iris(cube):
    """Regrid the cubes using block averaging."""
    logger.info("Regridding using block averaging")

    block_size = 100

    combined_data = average_block(cube.data, block_size)

    # Define target latitude and longitude ranges
    target_lats = np.linspace(90 - 0.5 * (180 / combined_data.shape[1]),
                              -90 + 0.5 * (180 / combined_data.shape[1]),
                              combined_data.shape[1])
    target_lons = np.linspace(-180 + 0.5 * (360 / combined_data.shape[2]),
                              180 - 0.5 * (360 / combined_data.shape[2]),
                              combined_data.shape[2])

    # Flip the latitude points and data
    combined_data = combined_data[:, ::-1, :]

    lat_bounds = calculate_bounds(target_lats[::-1])
    lon_bounds = calculate_bounds(target_lons)

    combined_cube = iris.cube.Cube(combined_data,
                                   dim_coords_and_dims=[
                                       (cube.coord('time'), 0),
                                       (iris.coords.DimCoord(
                                           target_lats[::-1],
                                           standard_name='latitude',
                                           units='degrees',
                                           bounds=lat_bounds), 1),
                                       (iris.coords.DimCoord(
                                           target_lons,
                                           standard_name='longitude',
                                           units='degrees',
                                           bounds=lon_bounds), 2)])

    return combined_cube


def calculate_bounds(points):
    """Calculate bounds for a set of points."""
    bounds = np.zeros((len(points), 2))
    bounds[1:, 0] = (points[:-1] + points[1:]) / 2
    bounds[:-1, 1] = bounds[1:, 0]
    bounds[0, 0] = points[0] - (bounds[1, 0] - points[0])
    bounds[-1, 1] = points[-1] + (points[-1] - bounds[-2, 1])
    return bounds


def save_variable(cube, var, outdir, attrs, **kwargs):
    """Saver function.

    Saves iris cubes (data variables) in CMOR-standard named files.

    Parameters
    ----------
    cube: iris.cube.Cube
        data cube to be saved.

    var: str
        Variable short_name e.g. ts or tas.

    outdir: str
        root directory where the file will be saved.

    attrs: dict
        dictionary holding cube metadata attributes like
        project_id, version etc.

    **kwargs: kwargs
        Keyword arguments to be passed to `iris.save`
    """
    fix_dtype(cube)

    # Ensure the variable metadata is correctly set
    cube.var_name = var

    # Only set standard_name if it is a valid CF standard name
    if hasattr(cube, 'standard_name') and cube.standard_name in \
            iris.std_names.STD_NAMES:
        cube.standard_name = cube.standard_name
    else:
        cube.standard_name = None  # or keep it unset

    # Ensure long_name is set
    cube.long_name = cube.long_name or var

    # Debugging metadata before saving
    logger.debug("Saving cube with var_name=%s, standard_name=%s, "
                 "long_name=%s",
                 cube.var_name, cube.standard_name, cube.long_name)

    # CMOR standard
    try:
        time = cube.coord('time')
    except iris.exceptions.CoordinateNotFoundError:
        time_suffix = None
    else:
        year = f"{time.cell(0).point.year:d}"
        time_suffix = '-'.join([year])

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


def calculate_shrubfrac(shrub_cubes, glob_attrs, out_dir):
    if shrub_cubes:
        logger.info("Summing shrub cubes to create shrubFrac")
        shrub_fraction_cube = sum(shrub_cubes)
        regridded_shrub_fraction_cube = regrid_iris(shrub_fraction_cube)
        set_units(regridded_shrub_fraction_cube, '%')
        regridded_shrub_fraction_cube = fix_coords_esacci(
            regridded_shrub_fraction_cube)
        set_global_atts(regridded_shrub_fraction_cube, glob_attrs)
        save_variable(regridded_shrub_fraction_cube, "shrubFrac",
                      out_dir, glob_attrs, unlimited_dimensions=['time'])


def calculate_treefrac(tree_cubes, glob_attrs, out_dir):
    if tree_cubes:
        logger.info("Summing tree cubes to create treeFrac")
        tree_fraction_cube = sum(tree_cubes)
        regridded_tree_fraction_cube = regrid_iris(tree_fraction_cube)
        set_units(regridded_tree_fraction_cube, '%')
        regridded_tree_fraction_cube = fix_coords_esacci(
            regridded_tree_fraction_cube)
        set_global_atts(regridded_tree_fraction_cube, glob_attrs)
        save_variable(regridded_tree_fraction_cube, "treeFrac",
                      out_dir, glob_attrs, unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize data."""
    glob_attrs = cfg['attributes']
    if not start_date:
        start_date = datetime(1992, 1, 1)
    if not end_date:
        end_date = datetime(1992, 12, 31)

    shrub_vars = {'shrubs-bd', 'shrubs-be', 'shrubs-nd', 'shrubs-ne'}
    shrub_cubes = []
    tree_vars = {'trees-bd', 'trees-be', 'trees-nd', 'trees-ne'}
    tree_cubes = []

    for var_name, vals in cfg['variables'].items():
        if not isinstance(vals, dict):
            raise ValueError(f"Invalid format for variable {var_name}: "
                             f"{type(vals)}")
        var_info = cfg['cmor_table'].get_variable(vals['mip'], var_name)
        glob_attrs['mip'] = vals['mip']
        raw_info = {'name': vals['long_name']}
        inpfile_pattern = os.path.join(in_dir, vals['filename'])
        logger.info("CMORizing var %s from file type %s", var_name,
                    inpfile_pattern)

        for year in range(start_date.year, end_date.year + 1):
            year_inpfile_pattern = inpfile_pattern.format(year=year)
            inpfiles = sorted(glob.glob(year_inpfile_pattern))
            for inpfile in inpfiles:
                raw_info['file'] = inpfile
                cube_list = extract_variable(raw_info)
                for cube in cube_list:
                    if var_name in shrub_vars:
                        shrub_cubes.append(cube)
                    elif var_name in tree_vars:
                        tree_cubes.append(cube)
                    else:
                        logger.info("Regridding cube for %s", var_name)
                        regridded_cube = regrid_iris(cube)
                        fix_var_metadata(regridded_cube, var_info)
                        regridded_cube = fix_coords_esacci(regridded_cube)
                        set_global_atts(regridded_cube, glob_attrs)
                        save_variable(regridded_cube, var_name, out_dir,
                                      glob_attrs,
                                      unlimited_dimensions=['time'])
                del cube_list  # Free memory
                gc.collect()  # Explicitly call garbage collection

    calculate_shrubfrac(shrub_cubes, glob_attrs, out_dir)
    calculate_treefrac(tree_cubes, glob_attrs, out_dir)
