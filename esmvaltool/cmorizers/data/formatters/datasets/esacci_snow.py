"""ESMValTool CMORizer for ESACCI-SNOW data.

Tier
   Tier 2: other freely-available dataset.

Source
   ftp://anon-ftp.ceda.ac.uk/neodc/esacci/snow/data

Last access
   20240214

Download and processing instructions
   Download the data from:
        scfg/AVHRR-MERGED/v2.0/
        swe/MERGED/v2.0/
      Put all scfg files in a single directory called 'scfg'
      (no subdirectories with years or months).
      Put all swe files in a single directory called 'swe'
      (no subdirectories with years or months).
"""

import logging
from copy import deepcopy
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import cf_units
import iris
import numpy as np
from dask import array as da
from dateutil import relativedelta
from esmvalcore.preprocessor import regrid

from esmvaltool.cmorizers.data.utilities import (
    save_variable,
)

logger = logging.getLogger(__name__)


def _load_callback(raw_cube, field, _):
    """Use this callback to fix anything Iris tries to break."""
    # Remove ancilliary variables that cause issues with merging
    # and concatenation.
    # Using a callback function is faster than removing the ancillary
    # variables later.
    for ancillary_variable in raw_cube.ancillary_variables():
        raw_cube.remove_ancillary_variable(ancillary_variable.standard_name)


def _create_nan_cube(cube, year, month, day):
    """Create cube containing only nan from existing cube."""
    nan_cube = cube.copy()
    nan_cube.data = da.ma.masked_greater(cube.core_data(), -1e20)

    # Read dataset time unit and calendar from file
    dataset_time_unit = str(nan_cube.coord("time").units)
    dataset_time_calender = nan_cube.coord("time").units.calendar
    # Convert datetime
    newtime = datetime(year=year, month=month, day=day, tzinfo=ZoneInfo("UTC"))
    newtime = cf_units.date2num(
        newtime,
        dataset_time_unit,
        dataset_time_calender,
    )
    nan_cube.coord("time").points = float(newtime)

    return nan_cube


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


def _extract_variable(in_files, var, cfg, out_dir, year):
    logger.info(
        "CMORizing variable '%s' from input files '%s'",
        var["short_name"],
        ", ".join(in_files),
    )
    attributes = deepcopy(cfg["attributes"])
    attributes["mip"] = var["mip"]
    attributes["raw"] = var["raw"]
    definition = cfg["cmor_table"].get_variable(var["mip"], var["short_name"])

    # load all input files (1 year) into 1 cube
    # --> drop attributes that differ among input files
    cube_list = iris.load(in_files, var["raw"], callback=_load_callback)
    # global attributes to remove
    drop_attrs = [
        "source",
        "date_created",
        "history",
        "tracking_id",
        "id",
        "time_coverage_start",
        "time_coverage_end",
        "platform",
        "sensor",
        "keywords",
    ]
    # variable attributes to remove
    drop_var_attrs = [
        "flag_meanings",
        "flag_values",
        "grid_mapping",
        "actual_range",
        "ancillary_variables",
    ]
    for cube in cube_list:
        for attr in drop_attrs:
            if attr in cube.attributes:
                cube.attributes.pop(attr)
        for attr in drop_var_attrs:
            if attr in cube.attributes:
                cube.attributes.pop(attr)

    # make sure there is one cube for every day of the year
    # (print debug info about missing days and add cube with
    # nan to fill gaps

    full_list = iris.cube.CubeList()
    loop_date = datetime(year, 1, 1, tzinfo=ZoneInfo("UTC"))
    time_list = []

    for cube in cube_list:
        loncoord = cube.coord("longitude")
        latcoord = cube.coord("latitude")
        loncoord.points = np.round(loncoord.core_points(), 3)
        latcoord.points = np.round(latcoord.core_points(), 3)

    # create list of available days ('time_list')

    for cube in cube_list:
        timecoord = cube.coord("time")
        cubetime = timecoord.units.num2date(timecoord.points)
        time_list.append(cubetime)

    # create cube list for every day of the year by adding
    # cubes containing only nan to fill possible gaps

    while loop_date <= datetime(year, 12, 31, tzinfo=ZoneInfo("UTC")):
        date_available = False
        for _idx, cubetime in enumerate(time_list):
            if loop_date == cubetime:
                date_available = True
                break
        if date_available:
            full_list.append(cube_list[_idx])
        else:
            logger.debug("No data available for %s", str(loop_date))
            nan_cube = _create_nan_cube(
                cube_list[0],
                loop_date.year,
                loop_date.month,
                loop_date.day,
            )
            full_list.append(nan_cube)
        loop_date += relativedelta.relativedelta(days=1)

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
    # input variable for snc (sncf) reports 'percent' --> rename to '%'
    # input variable for snw (swe) reports 'mm' --> rename to 'kg m-2'
    cube.units = definition.units

    # Fix data type
    cube.data = cube.core_data().astype("float32")

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

    # remove flags and other invalid points
    if cube.var_name == "snc":
        cube.data = da.ma.masked_greater(cube.core_data(), 100)
    elif cube.var_name == "snw":
        cube.data = da.ma.masked_less(cube.core_data(), 0)

    # regridding from 0.05x0.05 to 0.5x0.5 to save space
    cube = regrid(cube, target_grid="0.5x0.5", scheme="area_weighted")
    cube.attributes.update(
        {
            "geospatial_lon_resolution": "0.5",
            "geospatial_lat_resolution": "0.5",
            "spatial_resolution": "0.5",
        },
    )

    # Save results
    logger.debug("Saving cube\n%s", cube)
    logger.debug("Setting time dimension to UNLIMITED while saving!")
    save_variable(
        cube,
        cube.var_name,
        out_dir,
        attributes,
        unlimited_dimensions=["time"],
    )
    logger.info("Finished CMORizing %s", ", ".join(in_files))


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize ESACCI-SNOW dataset."""
    glob_attrs = cfg["attributes"]

    logger.info(
        "Starting cmorization for tier%s OBS files: %s",
        glob_attrs["tier"],
        glob_attrs["dataset_id"],
    )
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)
    logger.info("CMORizing ESACCI-SNOW version %s", glob_attrs["version"])

    if start_date is None:
        start_date = datetime(1979, 1, 1, tzinfo=ZoneInfo("UTC"))
    if end_date is None:
        end_date = datetime(2019, 12, 31, tzinfo=ZoneInfo("UTC"))

    for short_name, var in cfg["variables"].items():
        if "short_name" not in var:
            var["short_name"] = short_name
        loop_date = start_date
        while loop_date <= end_date:
            in_files = [
                str(x.resolve())
                for x in Path(in_dir).glob(
                    var["file"].format(year=loop_date.year),
                )
            ]
            if not in_files:
                logger.info(
                    "%d: no data not found for variable %s",
                    loop_date.year,
                    short_name,
                )
            else:
                _extract_variable(in_files, var, cfg, out_dir, loop_date.year)

            loop_date += relativedelta.relativedelta(years=1)
