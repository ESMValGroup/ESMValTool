"""ESMValTool CMORizer for ESACCI-AEROSOL data.

Tier
   Tier 2: other freely-available dataset.

Source
   CCI CEDA ftp: ftp://anon-ftp.ceda.ac.uk/neodc/esacci/aerosol/data/
     ATSR2_SU/L3/v4.3/MONTHLY/ (1997-2002)
     AATSR_SU/L3/v4.3/MONTHLY/ (2003-2011)
     ATSR2_SU/L3/v4.3/DAILY/ (1997-2002)
     AATSR_SU/L3/v4.3/DAILY/ (2003-2011)
     Other years are not considered since they are not complete.
   Copernicus Climate Data Store (CDS):
       https://cds.climate.copernicus.eu/datasets/satellite-aerosol-properties?tab=download
       SLSTR (see downloading instructions below)

Last access
   20250811

Download and processing instructions
   AATSR/ATSR2
   -----------
   Download the following files from CEDA:
       ftp: ftp://anon-ftp.ceda.ac.uk/neodc/esacci/aerosol/data/
           ATSR2_SU/L3/v4.3/MONTHLY/YYYY/*.nc
           AATSR_SU/L3/v4.3/MONTHLY/YYYY/*.nc
           ATSR2_SU/L3/v4.3/DAILY/YYYY/MM/*.nc
           AATSR_SU/L3/v4.3/DAILY/YYYY/MM/*.nc
   and put all monthly files into one directory named '{version}-monthly'
   all daily files into one directory named '{version}-daily'
   (the version string has to match the definition in the config file
   cmor_config/ESACCI-AEROSOL.yml, i.e. version = 'AATSR-SU-v4.3').

   SLSTR
   -----
   Select the following from the CDS:
     Time aggregation: "Daily average" and "Monthly average"
     Variable: "Aerosol optical depth" and "Fine-mode aerosol optical depth"
     Sensor on satellite: "SLSTR on SENTINEL 3A" and "SLSTR on SENTINEL 3B"
     Algorithm": "SWANSEA (Swansea University)"
     Year: select all
     Month: select all
     Day: select all
     Version: "SLSTR", "v1.12"
   and put all monthly files into one directory named '{version}-monthly'
   all daily files into one directory named '{version}-daily'.
   (the version string has to match the definition in the config file
   cmor_config/ESACCI-AEROSOL.yml, i.e. version = 'SLSTR-SU-v1.12').

   Alternatively, use the automatic downloader (recommended):
     esmvaltool data download ESACCI-AEROSOL


Modification history
   20250811-lauer_axel: written.
"""

import datetime
import logging
from copy import deepcopy
from pathlib import Path

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
        year=year,
        month=month,
        day=day,
        tzinfo=datetime.UTC,
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


def _extract_variable(in_files, var, cfg, out_dir, is_daily):
    """Extract, process and save variable to netCDF."""
    logger.info(
        "CMORizing variable '%s' from input files '%s'",
        var["short_name"],
        ", ".join(in_files),
    )
    attributes = deepcopy(cfg["attributes"])
    attributes["mip"] = var["mip"]
    attributes["raw"] = var["raw"]
    cmor_table = CMOR_TABLES[attributes["project_id"]]
    definition = cmor_table.get_variable(var["mip"], var["short_name"])

    # load all input files (1 year) into 1 cube
    # --> drop attributes that differ among input files
    cube_list = iris.load(in_files, var["raw"])

    # (global) attributes to remove
    drop_attrs = [
        "tracking_id",
        "id",
        "time_coverage_start",
        "time_coverage_end",
        "date_created",
        "inputfilelist",
        "platform",
        "source",
    ]

    time_unit = "days since 1850-01-01 00:00:00"
    time_calendar = "standard"

    new_list = iris.cube.CubeList()

    for cube in cube_list:
        # get time from attributes (no time coordinate)
        time0 = cube.attributes["time_coverage_start"]
        year0 = int(time0[0:4])
        month0 = int(time0[4:6])
        day0 = int(time0[6:8])
        if is_daily:
            timestamp = datetime.datetime(
                year0,
                month0,
                day0,
                tzinfo=datetime.UTC,
            )
        else:
            timestamp = datetime.datetime(
                year0,
                month0,
                15,
                tzinfo=datetime.UTC,
            )
        time_coord = iris.coords.DimCoord(
            cf_units.date2num(timestamp, time_unit, time_calendar),
            standard_name="time",
            var_name="time",
            units=cf_units.Unit(time_unit, calendar=time_calendar),
        )
        newcube = iris.util.new_axis(cube)
        newcube.add_dim_coord(time_coord, 0)

        for attr in drop_attrs:
            if attr in newcube.attributes:
                newcube.attributes.pop(attr)

        new_list.append(newcube)

    # make sure there is one cube for every day (daily data) or
    # every month (monthly data) of the year
    # (print debug info about missing days/months and add cube with
    # nan to fill gaps

    full_list = iris.cube.CubeList()
    time_list = []

    for cube in new_list:
        loncoord = cube.coord("longitude")
        latcoord = cube.coord("latitude")
        loncoord.points = np.round(loncoord.core_points(), 3)
        latcoord.points = np.round(latcoord.core_points(), 3)

    # create list of available days/months ('time_list')

    for cube in new_list:
        timecoord = cube.coord("time")
        cubetime = timecoord.units.num2date(timecoord.points)
        time_list.append(cubetime)

    # create cube list for every day/month of the year by adding
    # cubes containing only nan to fill possible gaps

    if is_daily:
        loop_date = datetime.datetime(year0, 1, 1, tzinfo=datetime.UTC)
        while loop_date <= datetime.datetime(
            year0,
            12,
            31,
            tzinfo=datetime.UTC,
        ):
            date_available = False
            for idx, cubetime in enumerate(time_list):
                if loop_date == cubetime:
                    date_available = True
                    full_list.append(new_list[idx])
                    break
            if not date_available:
                logger.debug(
                    "No data available for %s",
                    loop_date.strftime("%Y-%m-%d"),
                )
                masked_cube = _create_masked_cube(
                    new_list[0],
                    loop_date.year,
                    loop_date.month,
                    loop_date.day,
                )
                full_list.append(masked_cube)
            loop_date += relativedelta.relativedelta(days=1)
    else:
        loop_date = datetime.datetime(year0, 1, 15, tzinfo=datetime.UTC)
        while loop_date <= datetime.datetime(
            year0,
            12,
            31,
            tzinfo=datetime.UTC,
        ):
            date_available = False
            for idx, cubetime in enumerate(time_list):
                if loop_date == cubetime:
                    date_available = True
                    full_list.append(new_list[idx])
                    break
            if not date_available:
                logger.debug(
                    "No data available for %s",
                    loop_date.strftime("%Y-%m"),
                )
                masked_cube = _create_masked_cube(
                    new_list[0],
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
    attributes["version"] = (
        f"{version}-DAILY" if is_daily else f"{version}-MONTHLY"
    )
    utils.save_variable(
        cube,
        cube.var_name,
        out_dir,
        attributes,
        unlimited_dimensions=["time"],
    )
    logger.info("Finished CMORizing %s", ", ".join(in_files))


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize ESACCI-AEROSOL dataset."""
    glob_attrs = cfg["attributes"]
    glob_version = glob_attrs["version"] if "version" in glob_attrs else ""

    logger.info(
        "Starting cmorization for tier%s OBS files: %s",
        glob_attrs["tier"],
        glob_attrs["dataset_id"],
    )
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)
    logger.info("CMORizing ESACCI-AEROSOL version %s", glob_attrs["version"])

    for var_name, var in cfg["variables"].items():
        # Define dataset-specific default time ranges
        if var_name.find("aatsr") >= 0:
            dataset_start = datetime.datetime(1997, 1, 1, tzinfo=datetime.UTC)
            dataset_end = datetime.datetime(2011, 12, 31, tzinfo=datetime.UTC)
        elif var_name.find("slstr") >= 0:
            dataset_start = datetime.datetime(2017, 1, 1, tzinfo=datetime.UTC)
            dataset_end = datetime.datetime(2022, 12, 31, tzinfo=datetime.UTC)
        else:
            errmsg = f"Unknown dataset for variable {var_name}"
            raise ValueError(errmsg)

        var["var_name"] = var_name

        glob_attrs["mip"] = var["mip"]
        if "version" in var:
            glob_attrs["version"] = var["version"]
        else:
            glob_attrs["version"] = glob_version

        # Adjust start and end dates if not provided
        start_date_x = start_date or dataset_start
        end_date_x = end_date or dataset_end

        # Ensure the requested date range falls within the dataset limits
        start_date_x = max(start_date_x, dataset_start)
        end_date_x = min(end_date_x, dataset_end)

        loop_date = start_date_x
        if "day" in var_name:
            logger.info("Input data for %s is daily data", var_name)
            daily = True
        else:
            logger.info("Input data for %s is monthly data", var_name)
            daily = False
        while loop_date <= end_date_x:
            freqstr = "daily" if daily else "monthly"
            filepattern = var["file"].format(year=loop_date.year)
            srcdir = Path(in_dir) / f"{glob_attrs['version']}-{freqstr}"
            in_files = [str(p) for p in srcdir.glob(pattern=filepattern)]
            if not in_files:
                logger.info(
                    "%d: no data not found for variable %s",
                    loop_date.year,
                    var_name,
                )
            else:
                _extract_variable(in_files, var, cfg, out_dir, daily)

            loop_date += relativedelta.relativedelta(years=1)
