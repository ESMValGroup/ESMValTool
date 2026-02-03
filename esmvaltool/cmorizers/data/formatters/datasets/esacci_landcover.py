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
    in Tier2/ESACCI-LANDCOVER

"""

import glob
import logging
import os
from datetime import datetime

import iris
import numpy as np

from esmvaltool.cmorizers.data.utilities import (
    add_typebare,
    fix_coords,
    fix_var_metadata,
    save_variable,
    set_global_atts,
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Enable the new split-attributes handling mode
iris.FUTURE.save_split_attrs = True


def average_block(data, block_size):
    """Average the data within each block of size block_size.

    Parameters
    ----------
    data : numpy.ndarray
        The input data array to be block averaged.
    block_size : int
        The size of the block used for averaging. The data is averaged
        within non-overlapping blocks of this size along the spatial dimensions
        (latitude and longitude).

    Returns
    -------
    numpy.ndarray
        The block-averaged data array.
    """
    shape = data.shape
    reshaped_data = data.reshape(
        shape[0],
        shape[1] // block_size,
        block_size,
        shape[2] // block_size,
        block_size,
    )
    averaged_data = reshaped_data.mean(axis=(2, 4))
    return averaged_data


def regrid_iris(cube):
    """Regrid the cubes using block averaging.

    Parameters
    ----------
    cube : iris.cube.Cube
        The input data cube to be regridded.

    Returns
    -------
    iris.cube.Cube
        The regridded data cube.

    Notes
    -----
    The block size is set to 100, which means the data will be averaged within
    non-overlapping blocks of 100x100 grid cells along the spatial dimensions.
    """
    logger.info("Regridding using block averaging")

    block_size = 100  # Number of grid cells to average in each block

    combined_data = average_block(cube.data, block_size)

    # Define target latitude and longitude ranges
    target_lats = np.linspace(
        90 - 0.5 * (180 / combined_data.shape[1]),
        -90 + 0.5 * (180 / combined_data.shape[1]),
        combined_data.shape[1],
    )
    target_lons = np.linspace(
        -180 + 0.5 * (360 / combined_data.shape[2]),
        180 - 0.5 * (360 / combined_data.shape[2]),
        combined_data.shape[2],
    )

    combined_cube = iris.cube.Cube(
        combined_data,
        dim_coords_and_dims=[
            (cube.coord("time"), 0),
            (
                iris.coords.DimCoord(
                    target_lats,
                    standard_name="latitude",
                    units="degrees",
                ),
                1,
            ),
            (
                iris.coords.DimCoord(
                    target_lons,
                    standard_name="longitude",
                    units="degrees",
                ),
                2,
            ),
        ],
    )

    combined_cube.coord("latitude").guess_bounds()
    combined_cube.coord("longitude").guess_bounds()

    return combined_cube


def regrid_fix(cube, glob_attrs, var_name, var_info):
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
        data cube regridded and with fixed coordinates.
    """
    logger.info("Regridding cube for %s", var_name)
    regridded_cube = regrid_iris(cube)
    fix_var_metadata(regridded_cube, var_info)
    regridded_cube = fix_coords(regridded_cube)
    set_global_atts(regridded_cube, glob_attrs)

    return regridded_cube


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize data."""
    glob_attrs = cfg["attributes"]
    if not start_date:
        start_date = datetime(1992, 1, 1)
    if not end_date:
        end_date = datetime(2020, 12, 31)

    for year in range(start_date.year, end_date.year + 1):
        inpfile_pattern = os.path.join(in_dir, cfg["filename"])
        year_inpfile_pattern = inpfile_pattern.format(year=year)
        inpfiles = sorted(glob.glob(year_inpfile_pattern))
        for inpfile in inpfiles:
            cubes = iris.load(inpfile)
            for var_name, vals in cfg["variables"].items():
                var_info = cfg["cmor_table"].get_variable(
                    vals["mip"],
                    var_name,
                )
                glob_attrs["mip"] = vals["mip"]
                glob_attrs["frequency"] = vals["frequency"]
                if var_name == "shrubFrac":
                    cube = (
                        cubes.extract_cube("SHRUBS-BD")
                        + cubes.extract_cube("SHRUBS-BE")
                        + cubes.extract_cube("SHRUBS-ND")
                        + cubes.extract_cube("SHRUBS-NE")
                    )
                elif var_name == "treeFrac":
                    cube = (
                        cubes.extract_cube("TREES-BD")
                        + cubes.extract_cube("TREES-BE")
                        + cubes.extract_cube("TREES-ND")
                        + cubes.extract_cube("TREES-NE")
                    )
                else:
                    cube = cubes.extract_cube(vals["long_name"])
                regridded_cube = regrid_fix(
                    cube,
                    glob_attrs,
                    var_name,
                    var_info,
                )
                if var_name == "baresoilFrac":
                    add_typebare(regridded_cube)
                save_variable(
                    regridded_cube,
                    var_name,
                    out_dir,
                    glob_attrs,
                    unlimited_dimensions=["time"],
                )
