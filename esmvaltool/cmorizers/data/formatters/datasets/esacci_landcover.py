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
import logging
from datetime import datetime
import iris
import numpy as np

from ...utilities import (
    fix_coords,
    fix_var_metadata,
    set_global_atts,
    add_typebare,
    save_variable,
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

    combined_cube = iris.cube.Cube(combined_data,
                                   dim_coords_and_dims=[
                                       (cube.coord('time'), 0),
                                       (iris.coords.DimCoord(
                                           target_lats[::-1],
                                           standard_name='latitude',
                                           units='degrees'), 1),
                                       (iris.coords.DimCoord(
                                           target_lons,
                                           standard_name='longitude',
                                           units='degrees'), 2)])

    combined_cube.coord('latitude').guess_bounds()
    combined_cube.coord('longitude').guess_bounds()

    return combined_cube


def regrid_fix(cube, vals, glob_attrs, var_name, var_info):
    """Regrid cube and fixes.

    Regrids the cube, fixes metadata, coordinates and glob_attrs.

    Parameters
    ----------
    cube: iris.cube.Cube
          Data cube to be regridded.

    vals: dict
          Variable long_name.

    glob_attrs: dict
          Dictionary holding cube metadata attributes.

    var_name: str
          Variable name.

    var_info: dict
          Dictionary holding cube metadata attributes.

    Returns
    -------
    cube: iris.cube.Cube
        data cube regreidded and with fixed coordinates.
    """
    logger.info("Regridding cube for %s", var_name)
    regridded_cube = regrid_iris(cube)
    fix_var_metadata(regridded_cube, var_info)
    regridded_cube = fix_coords(regridded_cube)
    set_global_atts(regridded_cube, glob_attrs)

    return regridded_cube


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize data."""
    glob_attrs = cfg['attributes']
    if not start_date:
        start_date = datetime(1992, 1, 1)
    if not end_date:
        end_date = datetime(1992, 12, 31)

    for year in range(start_date.year, end_date.year + 1):
        inpfile_pattern = os.path.join(in_dir, cfg['filename'])
        year_inpfile_pattern = inpfile_pattern.format(year=year)
        inpfiles = sorted(glob.glob(year_inpfile_pattern))
        for inpfile in inpfiles:
            cubes = iris.load(inpfile)
            for var_name, vals in cfg['variables'].items():
                var_info = cfg['cmor_table'].get_variable(vals['mip'],
                                                          var_name)
                glob_attrs['mip'] = vals['mip']
                if var_name == 'shrubFrac':
                    cube = cubes.extract_cube('SHRUBS-BD') + \
                        cubes.extract_cube('SHRUBS-BE') + \
                        cubes.extract_cube('SHRUBS-ND') + \
                        cubes.extract_cube('SHRUBS-NE')
                elif var_name == 'treeFrac':
                    cube = cubes.extract_cube('TREES-BD') + \
                        cubes.extract_cube('TREES-BE') + \
                        cubes.extract_cube('TREES-ND') + \
                        cubes.extract_cube('TREES-NE')
                else:
                    cube = cubes.extract_cube(vals['long_name'])

                    regridded_cube = regrid_fix(cube, vals, glob_attrs,
                                                var_name, var_info)
                if var_name == 'baresoilFrac':
                    add_typebare(regridded_cube)
                save_variable(regridded_cube, var_name, out_dir, glob_attrs,
                              unlimited_dimensions=['time'])
