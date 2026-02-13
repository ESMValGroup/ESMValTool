"""ESMValTool CMORizer for OceanSODA-ETHZ data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0220059/

Last access
    20240215

Download and processing instructions
    Download the file OceanSODA_ETHZ-v2023.OCADS.01_1982-2022.nc

"""

import logging
import warnings
from datetime import datetime
from pathlib import Path

import iris
from dask import array as da
from iris import NameConstraint

from esmvaltool.cmorizers.data import utilities as utils
from esmvaltool.diag_scripts.shared.iris_helpers import unify_time_coord

logger = logging.getLogger(__name__)


def _fix_coords(cube, cmor_info):
    """Fix coordinates."""
    # Dimensional coordinates
    if "time" in cmor_info.dimensions:
        unify_time_coord(cube, "days since 1950-01-01 00:00:00")

        # Move time points to center of month
        time_coord = cube.coord("time")
        old_dates = time_coord.units.num2date(time_coord.points)
        new_dates = [datetime(t.year, t.month, 15) for t in old_dates]
        time_coord.points = time_coord.units.date2num(new_dates)
    cube.coord("lat").standard_name = "latitude"
    cube.coord("lon").standard_name = "longitude"
    cube = utils.fix_coords(cube)

    # Scalar coordinates
    if cmor_info.short_name in ("fgco2", "spco2"):
        utils.add_scalar_depth_coord(cube)

    return cube


def _fix_data(cube, var):
    """Fix data."""
    if var == "areacello":
        cube.data = da.ma.masked_equal(cube.core_data(), 0.0)


def _fix_var_metadata(var_info, cmor_info, attrs, cube):
    """Fix variable metadata."""
    if "raw_units" in var_info:
        cube.units = var_info["raw_units"]

    # fgco2:
    # Convert from mol(CO2) to kgC (note that one CO2 molecule contains one C
    # atom) and fix wrong sign (the dataset reports sea->air flux, while CMOR
    # expects "positive into ocean")
    if cmor_info.short_name == "fgco2":
        cube.data = -cube.core_data() * 12.01  # molar mass of C [g/mol]
        cube.units *= "g mol-1"
        attrs["positive"] = "down"

    # co3os, dissicos, talkos:
    # The original units of these variables are mumol/kg. To convert to the
    # CMOR units mol/m3, we assume a constant sea water density of 1028 kg/m3,
    # which is approximately the sea water density for T=4°C, salinity=35PSU,
    # and p=0bar according to the UNESCO formula (UNESCO, 1981, Tenth report of
    # the joint panel on oceanographic tables and standards, UNESCO Technical
    # Papers in Marine Science, see
    # https://www.wkcgroup.com/tools-room/seawater-density-calculator/ and
    # https://link.springer.com/content/pdf/bbm:978-3-319-18908-6/1.pdf).
    if cmor_info.short_name in ("co3os", "dissicos", "talkos"):
        cube.data = cube.core_data() * 1028.0
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
        cube = iris.load_cube(filepath, NameConstraint(var_name=raw_var))

    # Fix data
    _fix_data(cube, var)

    # Fix variable metadata
    _fix_var_metadata(var_info, cmor_info, attrs, cube)

    # Fix coordinates
    cube = _fix_coords(cube, cmor_info)

    # Fix global metadata
    utils.set_global_atts(cube, attrs)

    # Save variable
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
        glob_attrs["comment"] = var_info.get("comment", "")
        glob_attrs["mip"] = var_info["mip"]
        cmor_info = cmor_table.get_variable(var_info["mip"], var)
        _extract_variable(var_info, cmor_info, glob_attrs, filepath, out_dir)
