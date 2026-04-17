"""ESMValTool cmorizer for Midas data.

Tier
    Tier 2

Source
    https://zenodo.org/records/4244106#.Y5ekg3bMKMo

Modification history:
    20260226-Fruttarol_Noah: Written
"""

import copy
import datetime
import logging
from pathlib import Path

import iris
import numpy as np
from esmvalcore.cmor.table import CMOR_TABLES
from iris import cube as iris_cube

from esmvaltool.cmorizers.data import utilities as utils

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


def exstract_month(year, month, month_group, raw_var, attrs, cmor_info):
    # extract month from each file and combine into one cube
    starting_date = datetime.date(1950, 1, 1)
    const = iris.Constraint(raw_var)
    month_cubes = iris_cube.CubeList()
    for day in range(1, len(month_group) + 1):
        path = month_group[day - 1]
        cubes = iris.load(path)
        cube = cubes.extract_cube(const)
        date = int((datetime.date(year, month, day) - starting_date).days)
        lb_date = int((datetime.date(year, month, day) - starting_date).days)
        ub_date = int(
            (
                datetime.date(year, month, day)
                + datetime.timedelta(days=1)
                - starting_date
            ).days
        )
        time_coord = iris.coords.AuxCoord(
            date,
            standard_name="time",
            units="days since 1950-01-01 00:00:00",
            bounds=[lb_date, ub_date],
        )
        cube.add_aux_coord(time_coord)

        # Promote scalar 'time' coord to a new dimension of length 1.
        lat_points = cube.coord("latitude").points
        lon_points = cube.coord("longitude").points

        if lat_points.ndim == 2:
            # MIDAS grid is regular: latitude varies by row, longitude by column
            if not np.allclose(lat_points, lat_points[:, [0]], equal_nan=True):
                raise ValueError(
                    "Latitude is truly 2D/curvilinear; cannot convert to DimCoord."
                )
            lats = lat_points[:, 0]
        else:
            lats = lat_points

        if lon_points.ndim == 2:
            if not np.allclose(lon_points, lon_points[[0], :], equal_nan=True):
                raise ValueError(
                    "Longitude is truly 2D/curvilinear; cannot convert to DimCoord."
                )
            lons = lon_points[0, :]
        else:
            lons = lon_points

        lat_bounds = np.column_stack((lats - 0.05, lats + 0.05))
        lon_bounds = np.column_stack((lons - 0.05, lons + 0.05))

        lat_1d = iris.coords.DimCoord(
            lats, standard_name="latitude", units="degrees", bounds=lat_bounds
        )
        lon_1d = iris.coords.DimCoord(
            lons, standard_name="longitude", units="degrees", bounds=lon_bounds
        )

        cube.remove_coord("latitude")
        cube.remove_coord("longitude")
        cube.add_dim_coord(lat_1d, 0)
        cube.add_dim_coord(lon_1d, 1)
        cube = iris.util.new_axis(cube, "time")
        cube.units = "1"

        utils.set_global_atts(cube, attrs)
        utils.fix_var_metadata(cube, cmor_info)
        utils.fix_dim_coordnames(cube)

        utils.fix_coords(cube)
        month_cubes.append(cube)
    return month_cubes


def _extract_variable(short_name, var, cfg, in_dir, out_dir):
    attrs = copy.deepcopy(cfg["attributes"])
    attrs["mip"] = var["mip"]
    raw_var = var.get("raw_name", short_name)

    cmor_table = CMOR_TABLES[attrs["project_id"]]
    cmor_info = cmor_table.get_variable(var["mip"], short_name)

    file_template = attrs["filename_fmt"]
    logger.info(
        f"CMORizing variable {short_name} from file(s) {file_template}"
    )
    filepaths = Path(in_dir).glob(file_template)
    filepaths = sorted(filepaths)
    year_groups = group_by_time(filepaths)
    for year, month_groups in year_groups.items():
        cubes = iris_cube.CubeList()
        for month, month_group in month_groups.items():
            month_cubes = exstract_month(
                year, month, month_group, raw_var, attrs, cmor_info
            )
            iris.util.equalise_attributes(month_cubes)
            month_cube = month_cubes.concatenate_cube()
            mean_month_cube = month_cube.collapsed("time", iris.analysis.MEAN)
            starting_date = datetime.date(1950, 1, 1)
            date = int((datetime.date(year, month, 15) - starting_date).days)
            mean_month_cube = iris.util.new_axis(mean_month_cube, "time")
            mean_month_cube.coord("time").points = [date]

            cubes.append(mean_month_cube)
        cube = cubes.concatenate_cube()
        if short_name in ["od550dust"]:
            wavelength_coord = iris.coords.AuxCoord(
                [550],
                standard_name="radiation_wavelength",
                var_name="wavelength",
                units="nm",
            )
            cube.add_aux_coord(wavelength_coord)
        logger.debug(f"cube: {cube}")
        utils.save_variable(
            cube, short_name, out_dir, attrs, unlimited_dimensions=["time"]
        )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    # get variable information from config file
    for short_name, var in cfg["variables"].items():
        _extract_variable(short_name, var, cfg, in_dir, out_dir)
