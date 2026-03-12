"""ESMValTool CMORizer for CAMS data.

Tier
   Tier 3

Source
   https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-global-greenhouse-gas-inversion?tab=form

Last access
   20240909

Download and processing instructions
   Select carbon dioxide, surface flux, surface air sample, monthly mean
    and download the year you require

"""

import logging
import os
import warnings
from datetime import datetime

import dask.array as da
import iris
from cf_units import Unit

from esmvaltool.cmorizers.data.utilities import (
    fix_coords,
    fix_var_metadata,
    save_variable,
    set_global_atts,
    set_units,
)

logger = logging.getLogger(__name__)


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


def add_timeunits(cube, filename):
    """Add timestamp to cube."""
    tmp = str.split(filename, "_")
    time_coord = _get_time_coord(int(tmp[-1][:4]), int(tmp[-1][4:6]))
    cube.add_aux_coord(time_coord)


def _calculate_flux(cube, filename, area_type):
    """Calculate flux (dividing by land/sea area) and mask land/sea."""
    # Get land/sea area fraction
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Ignoring netCDF variable '.*?' invalid units '.*?'",
            category=UserWarning,
            module="iris",
        )
        lsf_cube = iris.load_cube(filename, "lsf")
    lsf = lsf_cube.core_data()

    # Mask
    if area_type == "land":
        mask = lsf == 0.0
    elif area_type == "ocean":
        mask = lsf > 0
    cube.data = da.ma.masked_array(cube.core_data(), mask=mask)

    # Calculate flux (sign change since input data and CMOR use different
    # conventions)
    cube.data = -cube.core_data()

    cube.attributes["positive"] = "down"

    return cube


def fix_units(cube):
    """Fix units from invalid units through import."""
    set_units(cube, "kg m-2 month-1")
    cube.convert_units("kg m-2 s-1")
    del cube.attributes["invalid_units"]


def extract_variable(short_name, var, filename):
    """Extract variable."""
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Ignoring netCDF variable '.*?' invalid units '.*?'",
            category=UserWarning,
            module="iris",
        )
        cube = iris.load_cube(filename, var["varname"])
    if short_name == "sftof":
        cube.data = 100.0 * (1.0 - cube.core_data())
        cube.units = "%"
    elif short_name == "sftlf":
        cube.data = 100.0 * cube.core_data()
        cube.units = "%"
    elif short_name == "nbp":
        _calculate_flux(cube, filename, "land")
        add_timeunits(cube, filename)
        fix_units(cube)
    elif short_name == "fgco2":
        _calculate_flux(cube, filename, "ocean")
        add_timeunits(cube, filename)
        fix_units(cube)
    return cube


def _fix_depth(cube, short_name, var, cfg):
    """Fix metadata of cube."""
    cmor_info = cfg["cmor_table"].get_variable(var["mip"], short_name)
    if "depth0m" in cmor_info.dimensions:
        depth_coord = iris.coords.AuxCoord(
            0.0,
            var_name="depth",
            standard_name="depth",
            long_name="depth",
            units=Unit("m"),
            attributes={"positive": "down"},
        )
        cube.add_aux_coord(depth_coord, ())


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize data."""
    months = [f"{mo:02d}" for mo in range(1, 13)]
    fpattern = os.path.join(in_dir, cfg["filename"])

    # run the cmorization
    for short_name, var in cfg["variables"].items():
        var_info = cfg["cmor_table"].get_variable(var["mip"], short_name)
        var_cubes = iris.cube.CubeList()
        logger.info("CMORizing var %s from file type %s", short_name, fpattern)
        # fx files are time invariant
        if short_name in ["areacella", "areacello", "sftlf", "sftof"]:
            filename = fpattern.format(year=cfg["start_year"], month="01")
            cube = extract_variable(short_name, var, filename)
        else:
            for year in range(cfg["start_year"], cfg["end_year"] + 1):
                for month in months:
                    filename = fpattern.format(year=year, month=month)
                    var_cubes.append(
                        extract_variable(short_name, var, filename),
                    )

            cube = var_cubes.merge_cube()

        cube.var_name = short_name
        fix_coords(cube)
        _fix_depth(cube, short_name, var, cfg)
        fix_var_metadata(cube, var_info)
        attrs = cfg["attributes"]
        attrs["mip"] = var["mip"]
        set_global_atts(cube, attrs)

        save_variable(
            cube,
            short_name,
            out_dir,
            attrs,
            unlimited_dimensions=["time"],
        )
