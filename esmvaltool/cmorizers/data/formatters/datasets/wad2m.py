"""ESMValTool CMORizer for WAD2M data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://zenodo.org/record/5553187

Last access
    20250701

Download and processing instructions
    Download the file WAD2M_wetlands_2000-2020_025deg_Ver2.0.nc.

"""

import logging
import warnings
import zipfile
from datetime import datetime
from pathlib import Path

import iris
from cf_units import Unit
from iris import NameConstraint
from iris.coords import AuxCoord

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

    # Add typewetla coordinate
    typewetla_coord = AuxCoord(
        "wetland",
        var_name="typewetla",
        standard_name="area_type",
        long_name="Wetland",
        units=Unit("no unit"),
    )
    cube.add_aux_coord(typewetla_coord, ())

    cube = utils.fix_coords(cube)

    return cube


def _fix_var_metadata(var_info, cmor_info, cube):
    """Fix variable metadata."""
    if "raw_units" in var_info:
        cube.units = var_info["raw_units"]

    # wetlandFrac:
    # Convert from fraction to percentage
    cube.convert_units("%")

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

    # Fix variable metadata
    _fix_var_metadata(var_info, cmor_info, cube)

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
        unlimited_dimensions=["time"],
    )


def _unzip(filepath, out_dir):
    """Unzip `*.zip` file."""
    logger.info("Starting extraction of %s to %s", filepath, out_dir)
    with zipfile.ZipFile(filepath, "r") as zip_ref:
        zip_ref.extract(Path(filepath).stem, out_dir)
    extracted_file = Path(out_dir) / Path(filepath).stem
    logger.info("Succefully extracted file to %s", extracted_file)
    return extracted_file


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    cmor_table = cfg["cmor_table"]
    glob_attrs = cfg["attributes"]

    # Run the cmorization
    for var, var_info in cfg["variables"].items():
        logger.info("CMORizing variable '%s'", var)
        glob_attrs["comment"] = var_info.get("comment", "")
        glob_attrs["mip"] = var_info["mip"]
        cmor_info = cmor_table.get_variable(var_info["mip"], var)

        # Extract file from ZIP archive
        zip_file = Path(in_dir) / cfg["filename"]
        if not Path(zip_file).is_file():
            logger.debug("Skipping '%s', file '%s' not found", var, zip_file)
            continue
        logger.info("Found input file '%s'", zip_file)
        filepath = _unzip(zip_file, out_dir)

        # Extract and save variable file
        _extract_variable(var_info, cmor_info, glob_attrs, filepath, out_dir)

        # Remove extracted file
        Path(filepath).unlink()
        logger.info("Removed cached input file %s", filepath)
