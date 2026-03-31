"""ESMValTool cmorizer for Midas data.

Tier
    Tier 2

Source
    https://zenodo.org/records/4244106#.Y5ekg3bMKMo

Modification history:
    20260226-Fruttarol_Noah: Written
"""

import logging
import os
from pathlib import Path

import iris
import numpy as np
from esmvalcore.cmor.table import CMOR_TABLES

logger = logging.getLogger(__name__)


def group_by_time(filepaths):
    year_groups = {}
    for path in filepaths:
        str_path = str(path)
        year = int(str_path[-11:-7])
        if year not in year_groups:
            year_groups[year] = []
        year_groups[year].append(path)

    groups = {}

    for year, paths in year_groups.items():
        month_groups = {}
        for path in paths:
            str_path = str(path)
            month = int(str_path[-7:-5])
            if month not in month_groups:
                month_groups[month] = []
            month_groups[month].append(path)
        groups[year] = month_groups

    return groups


def exstract_month(year, month, month_group, raw_var):
    # extract month from each file and combine into one cube
    const = iris.Constraint(raw_var)
    month_cubes = iris.cube.CubeList()
    for day, path in month_group.values():
        cubes = iris.load_cube(path)
        cube = cubes.extract_cube(const)
        date = dt.date(year, month, day)
        lb_date = dt.date(year, month, day)
        ub_date = dt.date(year, month, day) + dt.timedelta(days=1)
        time_coord = iris.coords.DimCoord(
            points=[date],
            standard_name="time",
            units="date",
            bounds=[[lb_date, ub_date]],
        )
        cube.add_dim_coord(time_coord)
        month_cubes.append(cube)
    return month_cubes


def _extract_variable(short_name, var, cfg, in_dir, out_dir):
    attrs = np.copy.deepcopy(cfg["attributes"])
    attrs["mip"] = var["mip"]
    ver = attrs["version"]
    files = attrs["files"]
    raw_var = var.get("raw_name", short_name)

    cmor_table = CMOR_TABLES[attrs["project_id"]]
    cmor_info = cmor_table.get_variable(var["mip"], short_name)

    logger.info("CMORizing variable '%s' from file(s) '%s'", short_name, files)
    filepaths = list(Path(os.path.join(in_dir)).glob(files))
    year_groups = group_by_time(filepaths)
    for year, month_groups in year_groups.items():
        cubes = iris.cube.CubeList()
        for month, month_group in month_groups.items():
            month_cubes = exstract_month(month_group, raw_var)
            mean_month_cube = month_cubes.concatenate_cube().collapsed(
                "time", iris.analysis.MEAN
            )
            cubes.append(
                mean_month_cube
            )  # TODO TODO need to check if this is correct!!!!!!!!!!!!!!!!!


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    # get variable information from config file
    for short_name, var in cfg["variables"].items():
        _extract_variable(short_name, var, cfg, in_dir, out_dir)
