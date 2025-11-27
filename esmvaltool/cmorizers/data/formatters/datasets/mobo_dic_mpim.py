"""ESMValTool CMORizer for MOBO-DIC-MPIM data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0221526/

Last access
    20221103

Download and processing instructions
    Download the file MOBO-DIC_MPIM_monthly_clim.nc

"""

import logging
import warnings
from datetime import datetime
from pathlib import Path

import iris
import numpy as np
from cf_units import Unit
from dask import array as da
from iris import NameConstraint
from iris.coords import CellMethod, DimCoord

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


TIME_UNITS = Unit("days since 1950-01-01 00:00:00", calendar="standard")


def _callback_fix_missing_value(cube, field, _):
    """Create masked array from missing_value."""
    if hasattr(field.cf_data, "missing_value"):
        missing_value = float(field.cf_data.missing_value)
        logger.info("Fixing missing value (%f)", missing_value)
        cube.data = da.ma.masked_equal(cube.core_data(), missing_value)


def _fix_climatological_time(cube):
    """Fix climatology coordinate."""
    # Following the doc the covered time period of the climatology is
    # January 2004 to December 2017 (Use 2011 as the "mean" year). See
    # https://www.ncei.noaa.gov/access/metadata/landing-page/bin/
    # iso?id=gov.noaa.nodc%3A0221526
    time_points = TIME_UNITS.date2num(
        [datetime(2011, m, 15) for m in range(1, 13)],
    )
    time_bounds = [
        [datetime(2004, m, 1), datetime(2017, m + 1, 1)] for m in range(1, 12)
    ]
    time_bounds.append([datetime(2004, 12, 1), datetime(2018, 1, 1)])
    time_bounds = TIME_UNITS.date2num(time_bounds)

    # Add new time coordinate to cube
    time_coord = DimCoord(
        time_points,
        bounds=time_bounds,
        standard_name="time",
        long_name="time",
        var_name="time",
        units=TIME_UNITS,
        climatological=True,
    )
    cube.remove_coord("month of the year")
    cube.add_dim_coord(time_coord, 0)

    # Fix cell methods
    cube.add_cell_method(CellMethod("mean within years", coords=time_coord))
    cube.add_cell_method(CellMethod("mean over years", coords=time_coord))


def _fix_time(cube):
    """Fix time coordinate."""
    julian_day_coord = cube.coord("Julian Day")

    # Calculate bounds of new time coordinate
    # print(str(julian_day_coord.units))
    datetime_base = datetime.strptime(
        str(julian_day_coord.units).partition(" since ")[2],
        "%Y-%m-%d %H:%M:%S",
    )
    base_year = datetime_base.year
    base_month = datetime_base.month
    all_months = list(julian_day_coord.points.astype(int)) + [
        julian_day_coord.points.astype(int).max() + 1,  # 1 more month for bnds
    ]
    bounds_datetimes = [
        datetime(base_year + (m - 1) // 12, base_month + (m - 1) % 12, 1)
        for m in all_months
    ]
    time_bounds = np.stack(
        (
            TIME_UNITS.date2num(bounds_datetimes[:-1]),
            TIME_UNITS.date2num(bounds_datetimes[1:]),
        ),
        axis=-1,
    )

    # Calculate time points as mean of bounds
    time_points = np.mean(time_bounds, axis=1)

    # Add new time coordinate to cube
    time_coord = DimCoord(
        time_points,
        bounds=time_bounds,
        standard_name="time",
        long_name="time",
        var_name="time",
        units=TIME_UNITS,
    )
    cube.remove_coord("Julian Day")
    cube.add_dim_coord(time_coord, 0)


def _fix_var_metadata(var_info, cmor_info, cube):
    """Fix variable metadata.

    Note
    ----
    The original units of 'dissic' are mumol/kg. To convert to the CMOR units
    mol/m3, we assume a constant sea water density of 1032 kg/m3, which is
    approximately the sea water density for T=4°C, salinity=35PSU, and p=100bar
    according to the UNESCO formula (UNESCO, 1981, Tenth report of the joint
    panel on oceanographic tables and standards, UNESCO Technical Papers in
    Marine Science, see
    https://www.wkcgroup.com/tools-room/seawater-density-calculator/ and
    https://link.springer.com/content/pdf/bbm:978-3-319-18908-6/1.pdf).

    """
    if "raw_units" in var_info:
        cube.units = var_info["raw_units"]

    # Special conversion for dissic (see Note above)
    if cmor_info.short_name == "dissic":
        cube.data = cube.core_data() * 1032.0
        cube.units *= "kg m-3"

    cube.convert_units(cmor_info.units)

    utils.fix_var_metadata(cube, cmor_info)


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
            callback=_callback_fix_missing_value,
        )

    # Fix variable metadata
    _fix_var_metadata(var_info, cmor_info, cube)

    # Fix coordinates
    if cube.coords("month of the year"):  # MOBO-DIC-MPIM
        _fix_climatological_time(cube)
    elif cube.coords("Julian Day"):  # MOBO-DIC2004-2019
        _fix_time(cube)
    cube.coord("depth").units = "m"
    cube = utils.fix_coords(cube, overwrite_time_bounds=False)

    # Fix global metadata
    utils.set_global_atts(cube, attrs)

    # Save variable
    with warnings.catch_warnings():
        warnings.filterwarnings(
            action="ignore",
            message="WARNING: missing_value not used",
            category=UserWarning,
            module="iris",
        )
        utils.save_variable(
            cube,
            var,
            out_dir,
            attrs,
            local_keys=["comment", "positive"],
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
        glob_attrs["comment"] = var_info["comment"]
        glob_attrs["mip"] = var_info["mip"]
        cmor_info = cmor_table.get_variable(var_info["mip"], var)
        _extract_variable(var_info, cmor_info, glob_attrs, filepath, out_dir)
