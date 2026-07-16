# pylint: disable=unused-argument
# pylint: disable=too-many-arguments
# pylint: disable=too-many-function-args
# pylint: disable=R0917
# pylint: disable=E1121
"""ESMValTool CMORizer for IAP data.

Tier
   Tier 2: other freely-available dataset.

Source
   IAPv4.2: "http://www.ocean.iap.ac.cn/ftp/cheng/"
            "IAPv4.2_IAP_Temperature_gridded_1month_netcdf/Monthly/"

Last access: 20250220

Download and processing instructions
   All handled by the script (download only if local data are missing)

   Alternatively, download and unzip the following files:
     Temperature_IAPv4.2_gridded_data_1940_1949.zip
     Temperature_IAPv4.2_gridded_data_1950_1959.zip
     Temperature_IAPv4.2_gridded_data_1960_1969.zip
     Temperature_IAPv4.2_gridded_data_1970_1979.zip
     Temperature_IAPv4.2_gridded_data_1980_1989.zip
     Temperature_IAPv4.2_gridded_data_1990_1999.zip
     Temperature_IAPv4.2_gridded_data_2000_2009.zip
     Temperature_IAPv4.2_gridded_data_2010_2019.zip
     Temperature_IAPv4.2_gridded_data_2020_2023.zip
"""

import logging
import os
import warnings
from datetime import datetime
from warnings import catch_warnings

import iris
import numpy as np
from dateutil import relativedelta

from esmvaltool.cmorizers.data.utilities import (
    fix_coords,
    fix_var_metadata,
    save_variable,
    set_global_atts,
)

logger = logging.getLogger(__name__)

try:
    iris.FUTURE.date_microseconds = True
    iris.FUTURE.save_split_attrs = True
except AttributeError as e:
    # Handle cases where FUTURE or the attributes don't exist
    logger.warning("AttributeError: %s", e)
except (TypeError, ValueError) as e:
    # Handle specific errors if these might occur
    logger.warning("TypeError or ValueError: %s", e)


def collect_files(in_dir, cfg, start_date, end_date):
    """Create list of files path to be processed."""
    file_list = []

    if start_date is None:
        start_date = datetime(year=1940, month=1, day=1)
    if end_date is None:
        end_date = datetime(year=2024, month=12, day=31)

    loop_date = start_date

    while loop_date <= end_date:
        fname = (
            f"IAPv4_Temp_monthly_1_6000m_year_{loop_date.year}"
            f"_month_{loop_date.month:02d}.nc"
        )
        in_file = os.path.join(in_dir, fname)
        file_list.append(in_file)
        loop_date += relativedelta.relativedelta(months=1)

    return file_list


def process_data(cube):
    """Process raw data: concatenate the cubes and return the new cube."""
    # Add time dimension
    temperature_data = np.expand_dims(cube.data, axis=0)
    temperature_data = np.moveaxis(
        temperature_data,
        (0, 1, 2, 3),
        (0, 2, 3, 1),
    )  # Reorder axes

    # Create time coordinate
    start_date = datetime(
        int(cube.attributes["StartYear"]),
        int(cube.attributes["StartMonth"]),
        int(cube.attributes["StartDay"]),
    )
    reference_date = datetime(2000, 1, 1)
    time_points = [(start_date - reference_date).days]

    time_coord = iris.coords.DimCoord(
        time_points,
        standard_name="time",
        units=(
            f"days since {reference_date.year}-"
            f"{reference_date.month}-{reference_date.day}"
        ),
    )

    # Remove old date attributes
    for key in [
        "StartDay",
        "StartMonth",
        "StartYear",
        "EndDay",
        "EndMonth",
        "EndYear",
    ]:
        del cube.attributes[key]

    # Get existing coordinates and rename 'standard depth' to 'depth'
    latitude_coord = cube.coord("latitude")
    longitude_coord = cube.coord("longitude")
    depth_coord = cube.coord("standard depth")
    depth_coord.rename("depth")
    depth_coord.var_name = "lev"
    depth_coord.attributes["positive"] = "down"

    # Create and return the new cube
    return iris.cube.Cube(
        temperature_data,
        var_name="Temperature",
        dim_coords_and_dims=[
            (time_coord, 0),
            (depth_coord, 1),
            (latitude_coord, 2),
            (longitude_coord, 3),
        ],
        attributes=cube.attributes,
    )


def extract_variable(in_files, out_dir, attrs, raw_info, cmor_table):
    """Extract variables and create OBS dataset."""
    var = raw_info["var"]
    var_info = cmor_table.get_variable(raw_info["mip"], var)
    rawvar = raw_info["raw_var"]
    with catch_warnings():
        warnings.simplefilter("ignore")  # Ignore all warnings
        cubes = iris.load(in_files, rawvar)
        cubes = iris.cube.CubeList(
            [process_data(cube) for cube in cubes],
        )

    iris.util.equalise_attributes(cubes)
    cube = cubes.concatenate_cube()
    fix_var_metadata(cube, var_info)
    fix_coords(cube)
    set_global_atts(cube, attrs)
    save_variable(cube, var, out_dir, attrs, unlimited_dimensions=["time"])

    # derive ocean surface
    if "srf_var" in raw_info:
        var_info = cmor_table.get_variable(
            raw_info["mip"],
            raw_info["srf_var"],
        )
        logger.info("Extract surface OBS for %s", raw_info["srf_var"])
        level_constraint = iris.Constraint(cube.var_name, depth=1)
        cube_os = cube.extract(level_constraint)
        fix_var_metadata(cube_os, var_info)
        save_variable(
            cube_os,
            raw_info["srf_var"],
            out_dir,
            attrs,
            unlimited_dimensions=["time"],
        )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    cmor_table = cfg["cmor_table"]
    glob_attrs = cfg["attributes"]

    # run the cmorization
    for var, vals in cfg["variables"].items():
        in_files = collect_files(in_dir, cfg, start_date, end_date)
        logger.info("CMORizing var %s from input set %s", var, vals["name"])
        raw_info = cfg["variables"][var]
        raw_info.update(
            {
                "var": var,
                "reference_year": cfg["custom"]["reference_year"],
            },
        )
        glob_attrs["mip"] = vals["mip"]
        extract_variable(in_files, out_dir, glob_attrs, raw_info, cmor_table)
