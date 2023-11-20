"""ESMValTool CMORizer for Yang2020 data.
"""

from datetime import datetime
import logging
import pathlib

import iris
from cf_units import Unit
import iris.coords
import iris.cube
import netCDF4

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def build_cube(var, var_info, dataset):
    var_data = dataset.variables[var_info["raw"]][:, :, :]
    var_data *= var_info["factor"]

    lonobs = dataset.variables["lon"][:]
    latobs = dataset.variables["lat"][:]

    time_units = Unit("days since 1950-01-01 00:00:00", calendar="standard")
    time_points = time_units.date2num(
        [datetime(1999, m, 15) for m in range(1, 13)],
    )
    time_bounds = [
        [datetime(1979, m, 1), datetime(2018, m + 1, 1)] for m in range(1, 12)
    ]
    time_bounds.append([datetime(1979, 12, 1), datetime(2019, 1, 1)])
    time_bounds = time_units.date2num(time_bounds)

    time_coord = iris.coords.DimCoord(
        time_points,
        bounds=time_bounds,
        var_name="time",
        standard_name="time",
        long_name="time",
        units=time_units,
        climatological=True,
    )

    lon_coord = iris.coords.DimCoord(
        lonobs,
        var_name="lon",
        standard_name="longitude",
        long_name="Longitude",
        units="degrees_east",
    )
    lat_coord = iris.coords.DimCoord(
        latobs,
        var_name="lat",
        standard_name="latitude",
        long_name="latitude",
        units="degrees_north",
    )

    coord_spec = [
        (time_coord, 0),
        (lat_coord, 1),
        (lon_coord, 2),
    ]
    cube = iris.cube.Cube(
        var_data,
        var_name=var,
        dim_coords_and_dims=coord_spec,
        units=var_info["units"],
    )
    return cube


def extract_variable(var, var_info, cmor_table, attrs, filepath, out_dir):
    """Extract variable."""
    cmor_info = cmor_table.get_variable(var, var)

    with netCDF4.Dataset(filepath, "r") as dataset:
        cube = build_cube(var, var_info, dataset)
        utils.fix_var_metadata(cube, cmor_info)
        utils.convert_timeunits(cube, 1950)
        utils.fix_coords(cube)
        utils.set_global_atts(cube, attrs)
        utils.save_variable(
            cube,
            var,
            out_dir,
            attrs,
            unlimited_dimensions=["time"],
        )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    glob_attrs = cfg["attributes"]
    cmor_table = cfg["cmor_table"]

    for var, var_info in cfg["variables"].items():
        filepath = pathlib.Path(in_dir) / var_info["filename"]
        logger.info(f"CMORizing variable {var} from {filepath}")
        glob_attrs["mip"] = var_info["mip"]
        extract_variable(
            var, var_info, cmor_table, glob_attrs, filepath, out_dir
        )
