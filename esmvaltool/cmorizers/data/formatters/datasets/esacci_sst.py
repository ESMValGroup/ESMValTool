"""ESMValTool CMORizer for ESACCI-SST data.

Tier
   Tier 2: other freely-available dataset.

Source
   https://catalogue.ceda.ac.uk/uuid/4a9654136a7148e39b7feb56f8bb02d2/

Last access
   20241211

#Download and processing instructions
#   A donwnloader is provided by ESMValTool.
#   (esmvaltool/cmorizers/data/downloaders/esacci_sst.py)

"""

import copy
import glob
import logging
import os
from datetime import datetime

import iris
from esmvalcore.cmor.fixes import get_time_bounds
from esmvalcore.preprocessor import concatenate, regrid

from esmvaltool.cmorizers.data import utilities as utils

from ...utilities import (
    fix_coords,
    fix_var_metadata,
    save_variable,
    set_global_atts,
)

logger = logging.getLogger(__name__)


def extract_variable(raw_info):
    """Extract to all vars."""
    rawvar = raw_info["name"]
    constraint = iris.NameConstraint(var_name=rawvar)
    if rawvar == "analysed_sst_uncertainty":
        tmp_cube = iris.load_cube(
            raw_info["file"], iris.NameConstraint(var_name="analysed_sst")
        )
        ancillary_var = tmp_cube.ancillary_variable(
            "sea_water_temperature standard_error"
        )
        cube = tmp_cube.copy(ancillary_var.core_data())
    else:
        try:
            cube = iris.load_cube(raw_info["file"], constraint)
        except iris.exceptions.ConstraintMismatchError as constraint_error:
            raise ValueError(
                f"No data available for variable {rawvar} in file"
                f" {raw_info['file']}"
            ) from constraint_error

    # Remove ancillary data
    for ancillary_variable in cube.ancillary_variables():
        cube.remove_ancillary_variable(ancillary_variable)
    return cube


def get_monthly_cube(
    cfg, var, vals, raw_info, attrs, inpfile_pattern, year, month
):
    data_cubes = []
    month_inpfile_pattern = inpfile_pattern.format(
        year=str(year) + f"{month:02}"
    )
    logger.info("Pattern: %s", month_inpfile_pattern)
    inpfiles = sorted(glob.glob(month_inpfile_pattern))
    if inpfiles == []:
        logger.error(
            "Could not find any files with this pattern %s",
            month_inpfile_pattern,
        )
        raise ValueError
    logger.info("Found input files: %s", inpfiles)

    for inpfile in inpfiles:
        raw_info["file"] = inpfile
        logger.info(
            "CMORizing var %s from file type %s", var, raw_info["file"]
        )
        data_cubes.append(extract_variable(raw_info))

    cube = concatenate(data_cubes)

    # regridding from 0.05x0.05 to 0.5x0.5 (not for uncertainty field
    if "Stderr" not in var:
        cube = regrid(cube, target_grid="0.5x0.5", scheme="area_weighted")

    # Fix dtype
    utils.fix_dtype(cube)
    # Fix units
    cmor_info = cfg["cmor_table"].get_variable(vals["mip"][0], var)
    cube.convert_units(cmor_info.units)
    # Fix metadata
    fix_var_metadata(cube, cmor_info)
    # Fix coordinates
    fix_coords(cube)
    cube.coord("time").long_name = "time"
    cube.coord("latitude").long_name = "latitude"
    cube.coord("longitude").long_name = "longitude"
    # Fix monthly time bounds
    time = cube.coord("time")
    time.bounds = get_time_bounds(time, vals["frequency"])

    # set global attributes
    set_global_atts(cube, attrs)
    # add comment to tosStderr
    if var == "tosStderr":
        cube.attributes["comment"] = (
            "Note that the variable tsStderr is an "
            "uncertainty not a standard error."
        )

    return cube


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    glob_attrs = copy.deepcopy(cfg["attributes"])

    # run the cmorization
    for var, vals in cfg["variables"].items():
        if not start_date:
            start_date = datetime(vals["start_year"], 1, 1)
        if not end_date:
            end_date = datetime(vals["end_year"], 12, 31)
        raw_info = {"name": vals["raw"]}
        inpfile_pattern = os.path.join(in_dir, "{year}*" + vals["filename"])
        logger.info("CMORizing var %s from file type %s", var, inpfile_pattern)
        mon_cubes = []
        for year in range(start_date.year, end_date.year + 1):
            logger.info("Processing year %s", year)
            glob_attrs["mip"] = vals["mip"][0]
            for month in range(start_date.month, end_date.month + 1):
                monthly_cube = get_monthly_cube(
                    cfg,
                    var,
                    vals,
                    raw_info,
                    glob_attrs,
                    inpfile_pattern,
                    year,
                    month,
                )
                # Save daily data
                save_variable(
                    monthly_cube,
                    var,
                    out_dir,
                    glob_attrs,
                    unlimited_dimensions=["time"],
                )
                # Calculate monthly mean
                if "Stderr" not in var:
                    logger.info("Calculating monthly mean")
                    iris.coord_categorisation.add_month_number(
                        monthly_cube, "time"
                    )
                    iris.coord_categorisation.add_year(monthly_cube, "time")
                    monthly_cube = monthly_cube.aggregated_by(
                        ["month_number", "year"], iris.analysis.MEAN
                    )
                    monthly_cube.remove_coord("month_number")
                    monthly_cube.remove_coord("year")
                    mon_cubes.append(monthly_cube)
            # Save monthly data
            if "Stderr" not in var:
                yearly_cube = concatenate(mon_cubes)
                glob_attrs["mip"] = vals["mip"][1]
                save_variable(
                    yearly_cube,
                    var,
                    out_dir,
                    glob_attrs,
                    unlimited_dimensions=["time"],
                )
                mon_cubes.clear()
