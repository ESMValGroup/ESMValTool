"""ESMValTool CMORizer for GPCP-SG data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://psl.noaa.gov/data/gridded/data.gpcp.html
    https://downloads.psl.noaa.gov/Datasets/gpcp/precip.mon.mean.nc

Last access
    20230215

Download and processing instructions
    Download the file precip.mon.mean.nc
    wget https://downloads.psl.noaa.gov/Datasets/gpcp/precip.mon.mean.nc
"""

import logging
import warnings
from pathlib import Path

import iris
from iris import NameConstraint

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _fix_var_metadata(var_info, cmor_info, cube):
    """Fix variable metadata."""
    if "raw_units" in var_info:
        cube.units = var_info["raw_units"]

    if cube.units == "mm/day":
        cube.units = "kg m-2 day-1"

    cube.convert_units(cmor_info.units)

    utils.fix_var_metadata(cube, cmor_info)
    return cube


def _fix_coords(cube, filepath):
    """Fix coordinates."""
    utils.fix_dim_coordnames(cube)

    # Bounds

    # Time
    time_bnds = iris.load_cube(filepath, NameConstraint(var_name="time_bnds"))
    cube.coord("time").bounds = time_bnds.core_data()
    # Latitude
    lat_bnds = iris.load_cube(filepath, NameConstraint(var_name="lat_bnds"))
    cube.coord("latitude").bounds = lat_bnds.core_data()
    # Longitude
    lon_bnds = iris.load_cube(filepath, NameConstraint(var_name="lon_bnds"))
    cube.coord("longitude").bounds = lon_bnds.core_data()


def _extract_variable(var_info, cmor_info, attrs, filepath, out_dir):
    """Extract variable."""
    var = cmor_info.short_name
    raw_var = var_info.get("raw_name", var)

    # Load data
    with warnings.catch_warnings():
        warnings.filterwarnings(
            action="ignore",
            message="Skipping global attribute 'units': 'units' is not a "
            "permitted attribute",
            category=UserWarning,
            module="iris",
        )
        cube = iris.load_cube(filepath, NameConstraint(var_name=raw_var))

    # Fix variable metadata
    cube = _fix_var_metadata(var_info, cmor_info, cube)

    # Fix coordinates
    _fix_coords(cube, filepath)

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
