"""ESMValTool CMORizer for GPCC data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://opendata.dwd.de/climate_environment/GPCC/html/fulldata-monthly_v2018_doi_download.html
    https://opendata.dwd.de/climate_environment/GPCC/
        full_data_2018/full_data_monthly_v2018_[025 05 10 25].nc.gz
Last access
    20200225

Download and processing instructions
    Download the following files:
        full_data_monthly_{version}.nc.gz

Two files are generated per version, one with version_grid (i.e. v2018-25),
one with version-grid-numgauge1 (i.e. v2018-25-numgauge1), which is constrained
on holding gridpoint values relying on data from at least one station (i.e.
removing gridpoints solely relying on climatological infilling).
"""

import copy
import logging
import os
from warnings import catch_warnings, filterwarnings

import cftime
import iris
import numpy as np
from cf_units import Unit
from iris import NameConstraint

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _get_centered_timecoord(cube):
    """Fix time coordinate.

    Time points start at the beginning of month at 00:00:00.
    """
    time = cube.coord("time")
    times = time.units.num2date(time.points)

    # get bounds
    starts = [cftime.DatetimeNoLeap(c.year, c.month, 1) for c in times]
    ends = [
        cftime.DatetimeNoLeap(c.year, c.month + 1, 1)
        if c.month < 12
        else cftime.DatetimeNoLeap(c.year + 1, 1, 1)
        for c in times
    ]
    time.bounds = time.units.date2num(np.stack([starts, ends], -1))

    # get points
    time.points = [np.mean((t1, t2)) for t1, t2 in time.bounds]


def _extract_variable(short_name, var, version, cfg, filepath, out_dir):
    """Extract variable."""
    raw_var = var.get("raw", short_name)
    with catch_warnings():
        filterwarnings(
            action="ignore",
            message="Ignoring netCDF variable .* invalid units .*",
            category=UserWarning,
            module="iris",
        )
        cube = iris.load_cube(filepath, NameConstraint(var_name=raw_var))

    # Fix units (mm/month) -> 'kg m-2 month-1' -> 'kg m-2 s-1'
    cmor_info = cfg["cmor_table"].get_variable(var["mip"], short_name)
    cube.units = Unit(var.get("raw_units", short_name))
    cube.convert_units(cmor_info.units)

    # fix calendar type
    cube.coord("time").units = Unit(
        cube.coord("time").units.origin,
        calendar=var.get("calendar", short_name),
    )
    cube.coord("time").convert_units(
        Unit("days since 1950-1-1 00:00:00", calendar="gregorian")
    )

    # Fix coordinates
    # fix time
    _get_centered_timecoord(cube)

    # fix flipped latitude
    utils.flip_dim_coord(cube, "latitude")
    utils.fix_dim_coordnames(cube)
    cube_coord = cube.coord("latitude")
    utils.fix_bounds(cube, cube_coord)

    # fix longitude
    cube_coord = cube.coord("longitude")
    if cube_coord.points[0] < 0.0 and cube_coord.points[-1] < 181.0:
        cube_coord.points = cube_coord.points + 180.0
        utils.fix_bounds(cube, cube_coord)
        cube.attributes["geospatial_lon_min"] = 0.0
        cube.attributes["geospatial_lon_max"] = 360.0
        nlon = len(cube_coord.points)
        utils.roll_cube_data(cube, nlon // 2, -1)

    # Fix metadata
    attrs = cfg["attributes"]
    attrs["mip"] = var["mip"]
    attrs["version"] = version
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(
        cube, short_name, out_dir, attrs, unlimited_dimensions=["time"]
    )

    # build contrainted cube on numgauge < 1
    constraint_var = var.get("constraint", short_name)
    with catch_warnings():
        filterwarnings(
            action="ignore",
            message="Ignoring netCDF variable .* invalid units .*",
            category=UserWarning,
            module="iris",
        )
        constr_cube = iris.load_cube(
            filepath, NameConstraint(var_name=constraint_var)
        )

    # fix flipped latitude
    utils.flip_dim_coord(constr_cube, "latitude")
    utils.fix_dim_coordnames(constr_cube)
    cube_coord = constr_cube.coord("latitude")
    utils.fix_bounds(constr_cube, cube_coord)

    # fix longitude
    cube_coord = constr_cube.coord("longitude")
    if cube_coord.points[0] < 0.0 and cube_coord.points[-1] < 181.0:
        cube_coord.points = cube_coord.points + 180.0
        utils.fix_bounds(constr_cube, cube_coord)
        constr_cube.attributes["geospatial_lon_min"] = 0.0
        constr_cube.attributes["geospatial_lon_max"] = 360.0
        nlon = len(cube_coord.points)
        utils.roll_cube_data(constr_cube, nlon // 2, -1)

    cube.data = np.ma.masked_where(constr_cube.data < 1.0, cube.data)

    # Save variable
    attrs = copy.deepcopy(cfg["attributes"])
    attrs.update(
        {
            "comment": "constrained on gridpoint values being based on"
            "at least 1 station",
            "version": attrs["version"] + "-numgauge1",
        }
    )
    attrs["mip"] = var["mip"]

    utils.save_variable(
        cube, short_name, out_dir, attrs, unlimited_dimensions=["time"]
    )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    raw_filepath = os.path.join(in_dir, cfg["filename"])

    # Run the cmorization
    for version in cfg["attributes"]["version"].values():
        for short_name, var in cfg["variables"].items():
            raw_var = var.get("raw", short_name)
            filepath = raw_filepath.format(version=version, raw_name=raw_var)
            logger.info("CMORizing variable '%s'", short_name)
            _extract_variable(short_name, var, version, cfg, filepath, out_dir)
