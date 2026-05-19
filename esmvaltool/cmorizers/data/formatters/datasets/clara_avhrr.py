"""ESMValTool CMORizer for CLARA-AVHRR data.

Tier
   Tier 3: restricted dataset (registration required).

Source
   Copernicus Climate Data Store (CDS):
       https://cds.climate.copernicus.eu/datasets/satellite-cloud-properties?tab=download

Last access
   20251126

Download and processing instructions
   Select the following from the CDS:
       Product family: CLARA-A3
       Origin: EUMETSAT
       Variable: Cloud fraction, Cloud physical properties of the ice/liquid phase
       Climate data record type: TCDR
       Time aggregation: Daily mean / Monthly mean
       Year: select all
       Month: select all
       Day: select all
       Geographical area: Whole available region
       Put all daily files for one month (mm) of one year (yyyy) under a single
       directory "daily/<yyyymm>", monthly files for one year (yyyy) under "monthly/<yyyy>".
   Note: you must accept the terms of use for CC-BY, EUMETSAT CM SAF, ESA CCI to be able
         to download the data from the CDS

   Alternatively, use the automatic downloader (recommended):
     esmvaltool data download CLARA-AVHRR


Modification history
   20251126-lauer_axel: written.
"""

import datetime
import glob
import logging
import os
from copy import deepcopy

import cf_units
import iris
import numpy as np
from dask import array as da
from dateutil import relativedelta
from esmvalcore.cmor.table import CMOR_TABLES

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _create_masked_cube(cube, year, month, day):
    """Create cube containing only nan from existing cube."""
    masked_cube = cube.copy()
    masked_cube.data = da.ma.masked_greater(cube.core_data(), -1e20)

    # Read dataset time unit and calendar from file
    dataset_time_unit = str(masked_cube.coord("time").units)
    dataset_time_calender = masked_cube.coord("time").units.calendar
    # Convert datetime
    newtime = datetime.datetime(
        year=year, month=month, day=day, tzinfo=datetime.UTC
    )
    newtime = cf_units.date2num(
        newtime,
        dataset_time_unit,
        dataset_time_calender,
    )
    masked_cube.coord("time").points = float(newtime)

    return masked_cube


def _fix_coordinates(cube, definition):
    """Fix coordinates."""
    axis2def = {"T": "time", "X": "longitude", "Y": "latitude"}
    axes = ["T", "X", "Y"]

    for axis in axes:
        coord_def = definition.coordinates.get(axis2def[axis])
        if coord_def:
            coord = cube.coord(axis=axis)
            if axis == "T":
                coord.convert_units("days since 1850-1-1 00:00:00.0")
            coord.standard_name = coord_def.standard_name
            coord.var_name = coord_def.out_name
            coord.long_name = coord_def.long_name
            coord.points = coord.core_points().astype("float64")
            if len(coord.points) > 1:
                if coord.bounds is not None:
                    coord.bounds = None
                coord.guess_bounds()

    return cube


def _extract_variable(cube_list, var, cfg, out_dir, is_daily):
    timefreq = "daily" if is_daily else "monthly"
    logger.info("CMORizing variable '%s' (%s)", var["short_name"], timefreq)

    attributes = deepcopy(cfg["attributes"])
    attributes["mip"] = var["mip"]
    attributes["raw"] = var["raw"]
    cmor_table = CMOR_TABLES[attributes["project_id"]]
    definition = cmor_table.get_variable(var["mip"], var["short_name"])

    # make sure there is one cube for every day (daily data) or
    # every month (monthly data) of the year
    # (print debug info about missing days/months and add cube with
    # nan to fill gaps

    full_list = iris.cube.CubeList()
    time_list = []

    # round latitude and longitude points to 3 digits to avoid rounding issues
    for cube in cube_list:
        loncoord = cube.coord("longitude")
        latcoord = cube.coord("latitude")
        loncoord.points = np.round(loncoord.core_points(), 3)
        latcoord.points = np.round(latcoord.core_points(), 3)

    # create list of available days/months ('time_list')
    year0 = 0
    for cube in cube_list:
        timecoord = cube.coord("time")
        if year0 == 0:
            year0 = timecoord.units.num2date(timecoord.points[0]).year
        cubetime = timecoord.units.num2date(timecoord.points)
        time_list.append(cubetime)

    # create cube list for every day/month of the year by adding
    # cubes containing only nan to fill possible gaps

    if is_daily:
        loop_date = datetime.datetime(year0, 1, 1, tzinfo=datetime.UTC)
        while loop_date <= datetime.datetime(
            year0, 12, 31, tzinfo=datetime.UTC
        ):
            date_available = False
            for idx, cubetime in enumerate(time_list):
                if (
                    loop_date.year == cubetime[0].year
                    and loop_date.month == cubetime[0].month
                    and loop_date.day == cubetime[0].day
                ):
                    date_available = True
                    full_list.append(cube_list[idx])
                    break
            if not date_available:
                logger.debug(
                    "No data available for %s",
                    loop_date.strftime("%Y-%m-%d"),
                )
                masked_cube = _create_masked_cube(
                    cube_list[0],
                    loop_date.year,
                    loop_date.month,
                    loop_date.day,
                )
                full_list.append(masked_cube)
            loop_date += relativedelta.relativedelta(days=1)
    else:
        loop_date = datetime.datetime(year0, 1, 1, tzinfo=datetime.UTC)
        while loop_date <= datetime.datetime(
            year0, 12, 31, tzinfo=datetime.UTC
        ):
            date_available = False
            for idx, cubetime in enumerate(time_list):
                if (
                    loop_date.year == cubetime[0].year
                    and loop_date.month == cubetime[0].month
                ):
                    date_available = True
                    full_list.append(cube_list[idx])
                    break
            if not date_available:
                logger.debug(
                    "No data available for %s",
                    loop_date.strftime("%Y-%m"),
                )
                masked_cube = _create_masked_cube(
                    cube_list[0],
                    loop_date.year,
                    loop_date.month,
                    loop_date.day,
                )
                full_list.append(masked_cube)
            loop_date += relativedelta.relativedelta(months=1)

    iris.util.unify_time_units(full_list)
    cube = full_list.concatenate_cube()
    cube.coord("time").points = (
        cube.coord("time").core_points().astype("float64")
    )

    # Set correct names
    cube.var_name = definition.short_name
    cube.standard_name = definition.standard_name
    cube.long_name = definition.long_name

    # Fix units
    cube.units = definition.units

    # Roll longitude
    cube.coord("longitude").points = cube.coord("longitude").points + 180.0
    nlon = len(cube.coord("longitude").points)
    cube.data = da.roll(cube.core_data(), int(nlon / 2), axis=-1)
    cube.attributes.update(
        {"geospatial_lon_min": "0", "geospatial_lon_max": "360"},
    )

    # Fix coordinates
    cube = _fix_coordinates(cube, definition)
    cube.coord("latitude").attributes = None
    cube.coord("longitude").attributes = None

    # Save results
    logger.debug("Saving cube\n%s", cube)
    logger.debug("Setting time dimension to UNLIMITED while saving!")
    version = attributes["version"]
    if is_daily:
        attributes["version"] = f"{version}-DAILY"
    else:
        attributes["version"] = f"{version}-MONTHLY"

    utils.save_variable(
        cube,
        cube.var_name,
        out_dir,
        attributes,
        unlimited_dimensions=["time"],
    )

    logger.info(
        "Finished CMORizing variable '%s' (%s) for current year",
        var["short_name"],
        timefreq,
    )


def _load_files(var, in_dir, year, daily):
    """Load all input files for one year. If requested, add different variables."""
    varlist = (
        var["raw"] if isinstance(var["raw"], list) else var["raw"].split()
    )
    filelist = (
        var["filename"]
        if isinstance(var["filename"], list)
        else var["filename"].split()
    )

    # create a list of filenames to be read

    in_files = []

    for filemask in filelist:
        if daily:
            for month in range(1, 13):
                filepattern = os.path.join(
                    in_dir,
                    f"daily/{year}{month:02d}",
                    filemask.format(year=year, month=f"{month:02d}"),
                )
                in_files.extend(glob.glob(filepattern))
        else:
            filepattern = os.path.join(
                in_dir,
                f"monthly/{year}",
                filemask.format(year=year),
            )
            in_files.extend(glob.glob(filepattern))

    if len(varlist) == 1:
        cube_list = iris.load(in_files, varlist[0])
    else:
        cube_list = []
        for raw_name in varlist:
            cube_list.extend(iris.load(in_files, raw_name))

    # (global) attributes to remove
    drop_attrs = [
        "date_created",
        "time_coverage_start",
        "time_coverage_end",
        "CMSAF_included_Daily_Means",
        "CMSAF_platform_and_orbits",
        "platform",
    ]

    # remove global attributes that might prevent concatenation

    for cube in cube_list:
        for attr in drop_attrs:
            if attr in cube.attributes:
                cube.attributes.pop(attr)

    # If "operator" is defined in the CMOR config file, then
    # do calculations now. So far, only the "sum" of 2 or more
    # variables is implemented

    cube_list_sum = []
    if var.get("operator", "") == "sum":
        for raw_name in varlist:
            sublist = [c for c in cube_list if c.var_name == raw_name]
            if not cube_list_sum:
                cube_list_sum = sublist
            else:
                logger.debug("Adding cubes (%s)...\n", raw_name)
                for cube in sublist:
                    # get time of cube to be added to the sum
                    timecoord = cube.coord("time")
                    # find sum cube with matching time
                    for sumcube in cube_list_sum:
                        sumtimecoord = sumcube.coord("time")
                        if timecoord == sumtimecoord:
                            result = sumcube
                            result += cube
                            logger.debug(
                                "cube added for time %s",
                                timecoord.units.num2date(timecoord.points),
                            )
                            break
        cube_list = cube_list_sum
    elif var.get("operator"):
        raise ValueError(
            "Multiple input files found, with operator '{}' configured: {}".format(
                var.get("operator"),
                ", ".join(in_files),
            ),
        )

    return cube_list


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize CLARA-AVHRR dataset."""
    glob_attrs = cfg["attributes"]
    glob_version = glob_attrs["version"] if "version" in glob_attrs else ""

    logger.info(
        "Starting cmorization for tier%s OBS files: %s",
        glob_attrs["tier"],
        glob_attrs["dataset_id"],
    )
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)
    logger.info("CMORizing CLARA-AVHRR version %s", glob_attrs["version"])

    if start_date is None:
        start_date_mm = datetime.datetime(1979, 1, 1, tzinfo=datetime.UTC)
        start_date_dd = datetime.datetime(2020, 1, 1, tzinfo=datetime.UTC)
    else:
        start_date_mm = start_date
        start_date_dd = start_date

    if end_date is None:
        end_date_mm = datetime.datetime(2020, 12, 31, tzinfo=datetime.UTC)
        end_date_dd = datetime.datetime(2020, 12, 31, tzinfo=datetime.UTC)
    else:
        end_date_mm = end_date
        end_date_dd = end_date

    for var_name, var in cfg["variables"].items():
        var["var_name"] = var_name

        glob_attrs["mip"] = var["mip"]
        if "version" in var:
            glob_attrs["version"] = var["version"]
        else:
            glob_attrs["version"] = glob_version

        if "day" in var_name:
            logger.info("Input data for %s is daily data", var_name)
            daily = True
            start_date = start_date_dd
            end_date = end_date_dd
        else:
            logger.info("Input data for %s is monthly data", var_name)
            daily = False
            start_date = start_date_mm
            end_date = end_date_mm

        for year in range(start_date.year, end_date.year + 1):
            logger.info("Processing year %s", year)
            cube_list = _load_files(var, in_dir, year, daily)

            if not cube_list:
                logger.info(
                    "%d: no data not found for variable %s",
                    year,
                    var_name,
                )
            else:
                _extract_variable(cube_list, var, cfg, out_dir, daily)
