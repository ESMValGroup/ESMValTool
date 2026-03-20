"""ESMValTool CMORizer for ESACCI-SEAICE data.

Tier
   Tier 2: other freely-available dataset.

Source
   ftp://anon-ftp.ceda.ac.uk/neodc/esacci/sea_ice/data

Last access
   20241107

Download and processing instructions
   Download the data from:
       sea_ice_concentration/L4/ssmi_ssmis/12.5km/v3.0/{year}/{month}
   Put all files under a single directory (no subdirectories
   with years or months).
"""

import glob
import logging
import os
from copy import deepcopy
from datetime import datetime

import cf_units
import iris
import numpy as np
from dask import array as da
from dateutil import relativedelta
from esmvalcore.cmor._fixes.common import OceanFixGrid
from esmvalcore.preprocessor import monthly_statistics
from iris.coords import AuxCoord

from esmvaltool.cmorizers.data import utilities as utils

from ...utilities import save_variable

logger = logging.getLogger(__name__)


def _create_nan_cube(cube, year, month, day):
    """Create cube containing only nan from existing cube."""
    nan_cube = cube.copy()
    nan_cube.data = da.ma.masked_greater(cube.core_data(), -1e20)

    # Read dataset time unit and calendar from file
    dataset_time_unit = str(nan_cube.coord("time").units)
    dataset_time_calender = nan_cube.coord("time").units.calendar
    # Convert datetime
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
        newtime,
        dataset_time_unit,
        dataset_time_calender,
    )
    nan_cube.coord("time").points = float(newtime_num)

    # remove existing time bounds and create new bounds
    coord = nan_cube.coord("time")
    bnd1 = newtime + relativedelta.relativedelta(hours=-12)
    bnd2 = bnd1 + relativedelta.relativedelta(days=1)
    coord.bounds = [
        cf_units.date2num(bnd1, dataset_time_unit, dataset_time_calender),
        cf_units.date2num(bnd2, dataset_time_unit, dataset_time_calender),
    ]

    return nan_cube


def _create_areacello(cfg, cube, glob_attrs, out_dir):
    var_info = cfg["cmor_table"].get_variable("Ofx", "areacello")
    glob_attrs["mip"] = "Ofx"
    lat_coord = cube.coord("latitude")

    arcube = iris.cube.Cube(
        np.zeros(lat_coord.shape, np.float32),
        standard_name=var_info.standard_name,
        long_name=var_info.long_name,
        var_name=var_info.short_name,
        units="m2",
        # time is index 0, add cell index dim
        dim_coords_and_dims=[(cube.coords()[1], 0), (cube.coords()[2], 1)],
    )

    # each grid cell is 12.5 km x 12.5 km
    arcube.data = arcube.core_data() + 12500 * 12500

    arcube.add_aux_coord(lat_coord, (0, 1))
    arcube.add_aux_coord(cube.coord("longitude"), (0, 1))
    utils.fix_var_metadata(arcube, var_info)
    utils.set_global_atts(arcube, glob_attrs)
    utils.save_variable(
        arcube,
        var_info.short_name,
        out_dir,
        glob_attrs,
        zlib=True,
    )


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
                coord.points = coord.core_points().astype("float64")
                if len(coord.points) > 1:
                    if coord.bounds is not None:
                        coord.bounds = None
                    coord.guess_bounds()
            coord.standard_name = coord_def.standard_name
            coord.var_name = coord_def.out_name
            coord.long_name = coord_def.long_name

    return cube


def _extract_variable(in_files, var, cfg, out_dir, year0, region):
    logger.info(
        "CMORizing variable '%s' from input files '%s'",
        var["short_name"],
        ", ".join(in_files),
    )
    attributes = deepcopy(cfg["attributes"])
    attributes["mip"] = var["mip_day"]
    attributes["raw"] = var["raw"]
    definition = cfg["cmor_table"].get_variable(
        var["mip_day"],
        var["short_name"],
    )

    # load all input files (1 year) into 1 cube
    # --> drop attributes that differ among input files
    cube_list = iris.load(in_files, var["raw"])

    # remove ancillary variables
    for cube in cube_list:
        for ancillary_variable in cube.ancillary_variables():
            cube.remove_ancillary_variable(ancillary_variable.standard_name)

    # (global) attributes to remove
    drop_attrs = [
        "tracking_id",
        "id",
        "time_coverage_start",
        "time_coverage_end",
        "date_created",
        "inputfilelist",
        "history",
        "valid_min",
        "valid_max",
    ]

    new_list = iris.cube.CubeList()

    for cube in cube_list:
        for attr in drop_attrs:
            if attr in cube.attributes.keys():
                cube.attributes.pop(attr)

        new_list.append(cube)

    # make sure there is one cube for every day of the year
    # (print debug info about missing days and add cube with
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
        ctnew = cubetime[0].replace(hour=0, minute=0, second=0, microsecond=0)
        time_list.append(ctnew)

    # create cube list for every day/month of the year by adding
    # cubes containing only nan to fill possible gaps

    loop_date = datetime(year0, 1, 1)
    while loop_date <= datetime(year0, 12, 31):
        for idx, cubetime in enumerate(time_list):
            if loop_date == cubetime:
                full_list.append(new_list[idx])
                break
        else:
            logger.debug(
                "No data available for %d/%d/%d",
                loop_date.month,
                loop_date.day,
                loop_date.year,
            )
            nan_cube = _create_nan_cube(
                new_list[0],
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
    cube.units = definition.units

    # Fix ocean-type grid (2-dim lat + lon)
    fixcube = OceanFixGrid(definition)
    cube = fixcube.fix_metadata(cubes=[cube])[0]

    # Fix coordinates
    cube = _fix_coordinates(cube, definition)
    cube.coord("latitude").attributes = None
    cube.coord("longitude").attributes = None

    # add aux coord 'typesi'
    area_type = AuxCoord(
        [1.0],
        standard_name="area_type",
        var_name="type",
        long_name="Sea Ice area type",
    )
    cube.add_aux_coord(area_type)

    # add attribute cell_measures
    cube.attributes.locals["cell_measures"] = "area: areacello"

    # Fix data type
    cube.data = cube.core_data().astype("float32")

    # save daily results
    logger.debug("Saving cube\n%s", cube)
    logger.debug("Setting time dimension to UNLIMITED while saving!")
    version = attributes["version"]
    attributes["version"] = f"{version}-{region}"
    save_variable(
        cube,
        cube.var_name,
        out_dir,
        attributes,
        unlimited_dimensions=["time"],
    )

    # calculate monthly means
    cube = monthly_statistics(cube, operator="mean")
    # Remove monthly statistics aux coordinates
    cube.remove_coord(cube.coord("month_number"))
    cube.remove_coord(cube.coord("year"))
    # save monthly results
    logger.debug("Saving cube\n%s", cube)
    logger.debug("Setting time dimension to UNLIMITED while saving!")
    version = attributes["version"]
    attributes["mip"] = var["mip_mon"]
    definition = cfg["cmor_table"].get_variable(
        var["mip_mon"],
        var["short_name"],
    )
    save_variable(
        cube,
        cube.var_name,
        out_dir,
        attributes,
        unlimited_dimensions=["time"],
    )

    # create and save areacello
    # (code adadapted from formatter 'nsidc_g02202_sh.py')
    _create_areacello(cfg, cube, attributes, out_dir)

    logger.info("Finished CMORizing %s", ", ".join(in_files))


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorize ESACCI-SEAICE dataset."""
    glob_attrs = cfg["attributes"]

    logger.info(
        "Starting cmorization for tier%s OBS files: %s",
        glob_attrs["tier"],
        glob_attrs["dataset_id"],
    )
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)
    logger.info("CMORizing ESACCI-SEAICE version %s", glob_attrs["version"])

    if start_date is None:
        start_date = datetime(1991, 1, 1)
    if end_date is None:
        end_date = datetime(2020, 12, 31)

    for short_name, var in cfg["variables"].items():
        if "short_name" not in var:
            var["short_name"] = short_name
        if "regions" not in var:
            regions = ("NH", "SH")
        else:
            regions = var["regions"]
        for region in regions:
            loop_date = start_date
            while loop_date <= end_date:
                filepattern = os.path.join(
                    in_dir,
                    region,
                    var["file"].format(year=loop_date.year, region=region),
                )
                in_files = glob.glob(filepattern)
                if not in_files:
                    logger.info(
                        "%d: no data not found for variable %s",
                        loop_date.year,
                        short_name,
                    )
                else:
                    _extract_variable(
                        in_files,
                        var,
                        cfg,
                        out_dir,
                        loop_date.year,
                        region,
                    )

                loop_date += relativedelta.relativedelta(years=1)
