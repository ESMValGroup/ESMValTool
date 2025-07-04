"""ESMValTool CMORizer for ESACCI-CLOUD data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://public.satproj.klima.dwd.de/data/ESA_Cloud_CCI/CLD_PRODUCTS/v3.0/L3U/AVHRR-PM/

Last access
    20230619

Download and processing instructions
    see downloading script
"""

import copy
import glob
import logging
import os
from calendar import monthrange
from datetime import datetime

import cf_units
import iris
from dask import array as da
from esmvalcore.preprocessor import (
    daily_statistics,
    monthly_statistics,
    regrid,
)
from iris import NameConstraint

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _create_nan_cube(cube, year, month, day):
    """Create cube containing only nan from existing cube."""
    nan_cube = cube.copy()
    nan_cube.data = da.ma.masked_greater(cube.core_data(), -1e20)

    # Read dataset time unit and calendar from file
    dataset_time_unit = str(nan_cube.coord("time").units)
    dataset_time_calender = nan_cube.coord("time").units.calendar
    newtime = datetime(
        year=year,
        month=month,
        day=day,
        hour=12,
        minute=0,
        second=0,
        microsecond=0,
    )
    newtime_num = cf_units.date2num(
        newtime, dataset_time_unit, dataset_time_calender
    )
    nan_cube.coord("time").points = float(newtime_num)

    return nan_cube


def _handle_missing_day(year, month, iday, short_name, cubes, cubes_day):
    """Fill missing day."""
    daily_cube = _create_nan_cube(cubes_day[0], year, month, iday)
    if short_name in ["clt", "ctp"]:
        cubes.append(daily_cube)
    cubes_day.append(daily_cube)


def _concatenate_and_save_daily_cubes(
    short_name, var, cubes, cubes_day, out_dir, cfg, cmor_info
):
    """Concatinate and save yearly cubes."""
    if short_name in ["clt", "ctp"]:
        cube = cubes.concatenate_cube()
    cube_day = cubes_day.concatenate_cube()

    # Calc daily
    if short_name in ["clt", "ctp"]:
        cube = daily_statistics(cube)
        cube.coord("time").points = [
            int(tpoint) + 0.5 for tpoint in cube.coord("time").points
        ]
    cube_day = daily_statistics(cube_day)
    cube_day.coord("time").points = [
        int(tpoint) + 0.5 for tpoint in cube_day.coord("time").points
    ]

    # Regridding from 0.05x0.05 to 0.5x0.5
    if short_name in ["clt", "ctp"]:
        cube = regrid(cube, target_grid="0.5x0.5", scheme="area_weighted")
    cube_day = regrid(cube_day, target_grid="0.5x0.5", scheme="area_weighted")

    # Fix units
    if short_name == "clt":
        cube.data = 100 * cube.core_data()
        cube_day.data = 100 * cube_day.core_data()
    else:
        if "raw_units" in var:
            if short_name == "ctp":
                cube.units = var["raw_units"]
            cube_day.units = var["raw_units"]
        if short_name == "ctp":
            cube.convert_units(cmor_info.units)
        cube_day.convert_units(cmor_info.units)

    if short_name in ["clt", "ctp"]:
        # Fix metadata and  update version information
        cube = utils.fix_coords(cube)
        utils.fix_var_metadata(cube, cmor_info)
        attrs = copy.deepcopy(cfg["attributes"])
        attrs["mip"] = var["mip"]
        attrs_day["version"] += "-AMPM"
        utils.set_global_atts(cube, attrs)
        # Save variable
        utils.save_variable(
            cube, short_name, out_dir, attrs, unlimited_dimensions=["time"]
        )

    # Fix metadata and  update version information
    cube_day = utils.fix_coords(cube_day)
    utils.fix_var_metadata(cube_day, cmor_info)
    attrs_day = copy.deepcopy(cfg["attributes"])
    attrs_day["mip"] = var["mip"]
    attrs_day["version"] += "-AMPM-daylight"
    utils.set_global_atts(cube_day, attrs_day)
    # Save variable
    utils.save_variable(
        cube_day,
        short_name,
        out_dir,
        attrs_day,
        unlimited_dimensions=["time"],
    )


def _concatenate_and_save_monthly_cubes(
    short_name, var, cubes, out_dir, attach, cfg, cmor_info
):
    """Concatinate monthly files and save."""
    # After gathering all cubes for all years, concatenate them
    cube = cubes.concatenate_cube()

    if attach == "-AMPM":
        cube = monthly_statistics(cube)

    # Regrid the cube to the target grid (e.g., 0.5x0.5)
    cube = regrid(cube, target_grid="0.5x0.5", scheme="area_weighted")

    # Fix units and handle any special cases like 'clt'
    if short_name == "clt":
        cube.data = 100 * cube.core_data()  # Example conversion
    else:
        if "raw_units" in var:
            cube.units = var["raw_units"]
        cube.convert_units(cmor_info.units)

    # Set global attributes and fix metadata
    utils.fix_var_metadata(cube, cmor_info)
    attrs = copy.deepcopy(cfg["attributes"])
    attrs["mip"] = var["mip"]
    attrs["version"] += attach
    utils.set_global_atts(cube, attrs)

    # Save the processed variable
    utils.save_variable(
        cube, short_name, out_dir, attrs, unlimited_dimensions=["time"]
    )


def _process_daily_file(
    ifile, inum, short_name, var, cmor_info, cubes, cubes_day
):
    """Extract variable from daily file."""
    logger.info("CMORizing file %s", ifile)
    # Extract raw names from the variable dictionary
    raw_var = var.get("raw", short_name)

    for ivar, raw_name in enumerate(raw_var):
        logger.info("Extracting raw variable %s", raw_name)

        # Define variable for daylight
        if "_asc" in raw_name:
            illum = "illum_asc"
        else:
            illum = "illum_desc"

        # Load cube using a constraint based on the raw_name
        daily_cube = iris.load_cube(ifile, NameConstraint(var_name=raw_name))
        daily_cube_ilum = iris.load_cube(ifile, NameConstraint(var_name=illum))

        # Set arbitrary time of day (in the end a daily mean is calculated)
        daily_cube.coord("time").points = (
            daily_cube.coord("time").points + (inum + 0.5 * ivar) * 0.1
        )
        daily_cube.attributes.clear()
        daily_cube.coord("time").long_name = "time"

        # Fix coordinates
        daily_cube = utils.fix_coords(daily_cube)
        # Fix dtype
        utils.fix_dtype(daily_cube)
        ## Fix metadata
        #utils.fix_var_metadata(daily_cube, cmor_info)

        # Check for daylight
        daily_cube_day = daily_cube.copy()
        daily_cube_day.data = da.ma.masked_where(
            daily_cube_ilum.core_data() > 1, daily_cube_day.core_data()
        )

        if short_name in ["clt", "ctp"]:
            cubes.append(daily_cube)
        cubes_day.append(daily_cube_day)


def _process_monthly_file(
    ifile, short_name, var, cmor_info, cubes_am, cubes_pm
):
    """Extract variable from monthly file."""
    logger.info("CMORizing file %s for variable %s", ifile, short_name)
    # Extract raw names from the variable dictionary
    raw_name = var.get("raw", short_name)
    # Try to load the cube using a constraint based on the raw_name
    monthly_cube = iris.load_cube(ifile, NameConstraint(var_name=raw_name))

    if short_name == "clwvi":
        logger.info("Adding lwp and clivi")
        cube_lwp = iris.load_cube(ifile, NameConstraint(var_name="lwp_allsky"))
        monthly_cube.data = monthly_cube.core_data() + cube_lwp.core_data()

    if monthly_cube is None:
        logger.warning("Cube could not be loaded for file '%s'", ifile)
        return  # Skip this file and move to the next

    monthly_cube.attributes.clear()
    monthly_cube.coord("time").long_name = "time"

    # Fix coordinates
    monthly_cube = utils.fix_coords(monthly_cube)
    # Fix data type
    utils.fix_dtype(monthly_cube)
    # Fix metadata
    utils.fix_var_metadata(monthly_cube, cmor_info)

    # Add the cube to the list
    try:
        if any(
            sat_am in ifile
            for sat_am in (
                "AVHRR_NOAA-12",
                "AVHRR_NOAA-15",
                "AVHRR_NOAA-17",
                "AVHRR_METOPA",
            )
        ):
            cubes_am.append(monthly_cube)
        elif any(
            sat_pm in ifile
            for sat_pm in (
                "AVHRR_NOAA-7",
                "AVHRR_NOAA-9",
                "AVHRR_NOAA-11",
                "AVHRR_NOAA-14",
                "AVHRR_NOAA-16",
                "AVHRR_NOAA-18",
                "AVHRR_NOAA-19",
            )
        ):
            cubes_pm.append(monthly_cube)
        else:
            raise ValueError(f"The file {ifile} is not assigned to AM or PM")
    except ValueError as e:
        logger.error(e)


def _extract_variable_daily(
    short_name, var, cfg, in_dir, out_dir, start_date, end_date
):
    """Extract daily variable."""
    cmor_info = cfg["cmor_table"].get_variable(var["mip"], short_name)

    if not start_date:
        start_date = datetime(cfg["start_year_daily"], 1, 1)
    if not end_date:
        end_date = datetime(cfg["end_year_daily"], 12, 31)

    for year in range(start_date.year, end_date.year + 1):
        cubes = iris.cube.CubeList()
        cubes_day = iris.cube.CubeList()

        for month in range(1, 13):
            num_days = monthrange(year, month)[1]
            for iday in range(1, num_days + 1):
                filelist = glob.glob(
                    os.path.join(
                        in_dir,
                        f"{year}{month:02}{iday:02}{var['file']}",
                    )
                )

                if filelist:
                    for inum, ifile in enumerate(filelist):
                        try:
                            _process_daily_file(
                                ifile,
                                inum,
                                short_name,
                                var,
                                cmor_info,
                                cubes,
                                cubes_day,
                            )
                        except Exception as e:
                            logger.error(
                                "Error processing file '%s': %s", ifile, e
                            )
                else:
                    _handle_missing_day(
                        year, month, iday, short_name, cubes, cubes_day
                    )

        _concatenate_and_save_daily_cubes(
            short_name, var, cubes, cubes_day, out_dir, cfg, cmor_info
        )


def _extract_variable_monthly(
    short_name, var, cfg, in_dir, out_dir, start_date, end_date
):
    """Extract monthly variable with improved handling for multiple cubes."""
    cmor_info = cfg["cmor_table"].get_variable(var["mip"], short_name)

    if not start_date:
        start_date = datetime(cfg["start_year_monthly"], 1, 1)
    if not end_date:
        end_date = datetime(cfg["end_year_monthly"], 12, 31)

    cubes_am = iris.cube.CubeList()
    cubes_pm = iris.cube.CubeList()

    for year in range(start_date.year, end_date.year + 1):
        for month in range(1, 13):  # Loop through all months (1-12)
            # Construct the file list for the current month
            filelist = glob.glob(
                os.path.join(in_dir, f"{year}{month:02}" + var["file"])
            )

            if not filelist:
                logger.warning(
                    "No monthly file found for %s-%02d", year, month
                )
                continue

            for ifile in filelist:
                try:
                    _process_monthly_file(
                        ifile, short_name, var, cmor_info, cubes_am, cubes_pm
                    )
                except Exception as e:
                    logger.error("Error processing file '%s': %s", ifile, e)

    if cubes_am:
        _concatenate_and_save_monthly_cubes(
            short_name, var, cubes_am, out_dir, "-AM", cfg, cmor_info
        )

    if cubes_pm:
        _concatenate_and_save_monthly_cubes(
            short_name, var, cubes_pm, out_dir, "-PM", cfg, cmor_info
        )

    if cubes_am and cubes_pm:
        # change day value in cubes_pm for concatinating
        for cube in cubes_am:
            time_coord = cube.coord("time")
            new_time_points = [tpoint + 1 for tpoint in time_coord.points]
            time_coord.points = new_time_points

        cubes_combined = cubes_am + cubes_pm
        _concatenate_and_save_monthly_cubes(
            short_name, var, cubes_combined, out_dir, "-AMPM", cfg, cmor_info
        )

    if not (cubes_am or cubes_pm):
        logger.warning("No valid cubes processed")


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """CMORization function call."""

    # Run the cmorization
    for var_name, var in cfg["variables"].items():
        short_name = var["short_name"]
        logger.info("CMORizing variable '%s'", var_name)
        if "L3U" in var["file"]:
            _extract_variable_daily(
                short_name, var, cfg, in_dir, out_dir, start_date, end_date
            )
        elif "L3C" in var["file"]:
            _extract_variable_monthly(
                short_name, var, cfg, in_dir, out_dir, start_date, end_date
            )
