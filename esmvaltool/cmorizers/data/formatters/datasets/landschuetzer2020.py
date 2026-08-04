"""ESMValTool CMORizer for Landschuetzer2020 data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0209633/

Last access
    20221102

Download and processing instructions
    Download the file MPI-ULB-SOM_FFN_clim.nc

"""

import logging
import warnings
from datetime import datetime
from pathlib import Path

import iris
from cf_units import Unit
from dask import array as da
from iris import NameConstraint
from iris.coords import CellMethod, DimCoord

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _callback_fix_fillvalue(cube, field, _):
    """Create masked array from FillValue."""
    if hasattr(field.cf_data, "FillValue"):
        fill_value = int(field.cf_data.FillValue)
        logger.info("Fixing fill value (%i)", fill_value)
        cube.data = da.ma.masked_equal(cube.core_data(), fill_value)


def _fix_climatological_time(cube):
    """Fix climatology coordinate."""
    time_units = Unit("days since 1950-01-01 00:00:00", calendar="standard")

    # Following the doc the covered time period of the climatology is
    # 1988-01-01 to 2020-01-01 (Use 2004 as the "mean" year). See
    # https://www.ncei.noaa.gov/access/metadata/landing-page/bin/
    # iso?id=gov.noaa.nodc%3A0209633
    time_points = time_units.date2num(
        [datetime(2004, m, 15) for m in range(1, 13)],
    )
    time_bounds = [
        [datetime(1988, m, 1), datetime(2019, m + 1, 1)] for m in range(1, 12)
    ]
    time_bounds.append([datetime(1988, 12, 1), datetime(2020, 1, 1)])
    time_bounds = time_units.date2num(time_bounds)

    # Add new time coordinate to cube
    time_coord = DimCoord(
        time_points,
        bounds=time_bounds,
        standard_name="time",
        long_name="time",
        var_name="time",
        units=time_units,
        climatological=True,
    )
    cube.remove_coord("time")
    cube.add_dim_coord(time_coord, 0)

    # Fix cell methods
    cube.add_cell_method(CellMethod("mean within years", coords=time_coord))
    cube.add_cell_method(CellMethod("mean over years", coords=time_coord))


def _fix_scalar_coords(cube):
    """Fix scalar coordinates."""
    if cube.var_name == "spco2":
        utils.add_scalar_depth_coord(cube)


def _extract_variable(var_info, cmor_info, attrs, filepath, out_dir):
    """Extract variable."""
    var = cmor_info.short_name
    raw_var = var_info.get("raw_name", var)

    # Load data
    with warnings.catch_warnings():
        warnings.filterwarnings(
            action="ignore",
            message="Ignoring netCDF variable .* invalid units .*",
            category=UserWarning,
            module="iris",
        )
        cube = iris.load_cube(
            filepath,
            NameConstraint(var_name=raw_var),
            callback=_callback_fix_fillvalue,
        )

    # Fix variable metadata
    if "raw_units" in var_info:
        cube.units = var_info["raw_units"]
    cube.convert_units(cmor_info.units)
    utils.fix_var_metadata(cube, cmor_info)

    # Fix coordinates
    _fix_climatological_time(cube)
    cube = utils.fix_coords(
        cube,
        overwrite_lat_bounds=False,
        overwrite_lon_bounds=False,
        overwrite_time_bounds=False,
    )
    _fix_scalar_coords(cube)

    # Fix global metadata
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(
        cube,
        var,
        out_dir,
        attrs,
        local_keys=["positive"],
        unlimited_dimensions=["time"],
    )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    cmor_table = cfg["cmor_table"]
    glob_attrs = cfg["attributes"]

    # Run the cmorization
    for var, var_info in cfg["variables"].items():
        filepath = Path(in_dir) / var_info["filename"]
        logger.info("CMORizing variable '%s' from file %s", var, filepath)
        glob_attrs["mip"] = var_info["mip"]
        cmor_info = cmor_table.get_variable(var_info["mip"], var)
        _extract_variable(var_info, cmor_info, glob_attrs, filepath, out_dir)
