"""ESMValTool CMORizer for NOAA-MBL-CH4 data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://gml.noaa.gov/ccgg/trends_ch4/

Last access
    20230717

Download and processing instructions
    Download the file:
    wget https://gml.noaa.gov/webdata/ccgg/trends/ch4/ch4_mm_gl.csv
"""

import logging
import warnings
from datetime import datetime
from pathlib import Path

import iris
import numpy as np
import pandas as pd
from cf_units import Unit

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)

LAT_COORD = iris.coords.DimCoord(
    [0.0],
    bounds=[[-90.0, 90.0]],
    var_name="lat",
    standard_name="latitude",
    long_name="latitude",
    units="degrees",
)
LON_COORD = iris.coords.DimCoord(
    [180.0],
    bounds=[[0.0, 360.0]],
    var_name="lon",
    standard_name="longitude",
    long_name="longitude",
    units="degrees",
)


def _fix_var_metadata(var_info, cmor_info, cube):
    """Fix variable metadata."""
    if "raw_units" in var_info:
        cube.units = var_info["raw_units"]

    cube.convert_units(cmor_info.units)

    utils.fix_var_metadata(cube, cmor_info)
    return cube


def _get_time_coord(year, month):
    """Get time coordinate."""
    point = datetime(year=year, month=month, day=15)
    bound_low = datetime(year=year, month=month, day=1)
    if month == 12:
        month_bound_up = 1
        year_bound_up = year + 1
    else:
        month_bound_up = month + 1
        year_bound_up = year
    bound_up = datetime(year=year_bound_up, month=month_bound_up, day=1)
    time_units = Unit("days since 1950-01-01 00:00:00", calendar="standard")
    time_coord = iris.coords.DimCoord(
        time_units.date2num(point),
        bounds=time_units.date2num([bound_low, bound_up]),
        var_name="time",
        standard_name="time",
        long_name="time",
        units=time_units,
    )
    return time_coord


def _get_cube(row, column_name):
    """Create :class:`iris.cube.Cube` from :class:`pandas.Series`."""
    time_coord = _get_time_coord(int(row["year"]), int(row["month"]))
    lat_coord = LAT_COORD.copy()
    lon_coord = LON_COORD.copy()
    data = np.ma.masked_invalid(row[column_name])
    cube = iris.cube.Cube(
        data.reshape((1, 1, 1)),
        dim_coords_and_dims=[(time_coord, 0), (lat_coord, 1), (lon_coord, 2)],
        units="ppb",
    )
    return cube


def _fix_coords(cube):
    """Fix coordinates."""
    utils.fix_dim_coordnames(cube)

    return cube


def _extract_variable(var_info, cmor_info, attrs, filepath, out_dir):
    """Extract variable."""
    var = cmor_info.short_name

    # Load data
    with warnings.catch_warnings():
        warnings.filterwarnings(
            action="ignore",
            message="Skipping global attribute 'units': 'units' is not a "
            "permitted attribute",
            category=UserWarning,
            module="iris",
        )
        skiprows = 0
        with open(filepath, encoding="utf-8") as csv:
            for line in csv:
                if line.startswith("#"):
                    skiprows = skiprows + 1

        data_frame = pd.read_csv(filepath, header=skiprows)

        # Extract cube
        cubes = iris.cube.CubeList()
        for _, row in data_frame.iterrows():
            cube = _get_cube(row, "average")
            cubes.append(cube)
        cube = cubes.concatenate_cube()
        cube.var_name = var

    # Fix coordinates
    cube = _fix_coords(cube)

    # Fix variable metadata
    cube = _fix_var_metadata(var_info, cmor_info, cube)

    # Fix global metadata
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(
        cube,
        var,
        out_dir,
        attrs,
        unlimited_dimensions=["time"],
    )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    cmor_table = cfg["cmor_table"]
    glob_attrs = cfg["attributes"]

    # Run the cmorization
    for var, var_info in cfg["variables"].items():
        filepath = Path(in_dir) / var_info["filename"]
        logger.info("CMORizing variable '%s' from file %s", var, filepath)
        glob_attrs["mip"] = var_info["mip"]
        cmor_info = cmor_table.get_variable(var_info["mip"], var)
        _extract_variable(var_info, cmor_info, glob_attrs, filepath, out_dir)
