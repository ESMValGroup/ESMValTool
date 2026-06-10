"""CMORize Yang2020 data."""

import logging
from datetime import datetime
from pathlib import Path

import iris
from cf_units import Unit
from iris import NameConstraint
from iris.coords import CellMethod, DimCoord

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


TIME_UNITS = Unit("days since 1950-01-01 00:00:00", calendar="standard")


def _fix_climatological_time(cube):
    """Fix climatology coordinate."""
    # Time bounds of the climatology are set to 1971-2018 since measurements
    # that serve as input for this data product have been conducted between
    # 1971 and 2018 (note that the wind products from ERA5 and CCMP, which are
    # necessary to derive the N2O flux, are only available from 1988-2017
    # though). This has been recommended by a co-author of the original
    # publication (https://doi.org/10.1073/pnas.1921914117) via email, but is
    # NOT explicitly stated in the original publication. Thus, the period given
    # here needs to be treated with care.
    time_points = TIME_UNITS.date2num(
        [datetime(1994, m, 15) for m in range(1, 13)],
    )
    time_bounds = [
        [datetime(1971, m, 1), datetime(2018, m + 1, 1)] for m in range(1, 12)
    ]
    time_bounds.append([datetime(1971, 12, 1), datetime(2019, 1, 1)])
    time_bounds = TIME_UNITS.date2num(time_bounds)

    # Add new time coordinate to cube
    time_coord = DimCoord(
        time_points,
        bounds=time_bounds,
        standard_name="time",
        long_name="time",
        var_name="time",
        units=TIME_UNITS,
        climatological=True,
    )
    cube.remove_coord("time")
    cube.add_dim_coord(time_coord, 0)

    # Fix cell methods
    cube.add_cell_method(CellMethod("mean within years", coords=time_coord))
    cube.add_cell_method(CellMethod("mean over years", coords=time_coord))


def _fix_var_metadata(var_info, cmor_info, cube):
    """Fix variable metadata."""
    if "raw_units" in var_info:
        cube.units = var_info["raw_units"]
    if (
        "g" in str(cube.units)
        and "mol" in cmor_info.units
        and "molar_mass" in var_info
    ):
        cube = cube / var_info["molar_mass"]
        cube.units = cube.units / "g mol-1"
    cube.convert_units(cmor_info.units)
    utils.fix_var_metadata(cube, cmor_info)
    return cube


def _extract_variable(var_info, cmor_info, attrs, filepath, out_dir):
    """Extract variable."""
    var = cmor_info.short_name
    raw_var = var_info.get("raw_name", var)

    cube = iris.load_cube(filepath, NameConstraint(var_name=raw_var))

    # Fix variable metadata
    cube = _fix_var_metadata(var_info, cmor_info, cube)

    # Fix coordinates
    if "time" in cmor_info.coordinates:
        _fix_climatological_time(cube)
    cube = utils.fix_coords(cube, overwrite_time_bounds=False)

    # Fix global metadata
    utils.set_global_atts(cube, attrs)
    if cmor_info.positive:
        cube.attributes.locals["positive"] = cmor_info.positive

    # Save variable
    utils.save_variable(
        cube,
        var,
        out_dir,
        attrs,
        local_keys=["comment", "positive"],
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
        if "comment" in var_info:
            glob_attrs["comment"] = var_info["comment"]
        cmor_info = cmor_table.get_variable(var_info["mip"], var)
        _extract_variable(var_info, cmor_info, glob_attrs, filepath, out_dir)
