"""ESMValTool CMORizer for NOAA ERSST data, version 3b.

   This is the CMORizer script for the NOAA Extended Reconstructed
   Sea Surface Temperature (ERSST) in its version 3b.

Tier
    Tier 2: open dataset.

Source
    https://doi.org/10.1175/1520-0442-16.10.1495

Last access
    20200520

Download and processing instructions
    The data is provided by NOAA at:
        https://www1.ncdc.noaa.gov/pub/data/cmb/ersst/v3b/netcdf/

"""

import logging
import os
import re

import iris
from cf_units import Unit

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _get_filepaths(in_dir, basename):
    """Find correct name of file (extend basename with timestamp)."""
    regex = re.compile(basename)
    return_files = []
    for files in os.listdir(in_dir):
        if regex.match(files):
            return_files.append(os.path.join(in_dir, files))

    return return_files


def _fix_time_coord(cube, _field, _filename):
    """Set time points to central day of month."""
    time_coord = cube.coord("time")
    new_unit = Unit("days since 1850-01-01 00:00:00", calendar="standard")
    time_coord.convert_units(new_unit)
    old_time = new_unit.num2date(time_coord.points)
    new_time = [d.replace(day=15) for d in old_time]
    time_coord.points = new_unit.date2num(new_time)


def _extract_variable(raw_var, cmor_info, attrs, filepath, out_dir):
    """Extract variable from all files."""
    var = cmor_info.short_name
    cubes = iris.load(filepath, raw_var, _fix_time_coord)
    iris.util.equalise_attributes(cubes)
    cube = cubes.concatenate_cube()
    cube = iris.util.squeeze(cube)

    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)
    utils.save_variable(
        cube, var, out_dir, attrs, unlimited_dimensions=["time"]
    )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    glob_attrs = cfg["attributes"]
    cmor_table = cfg["cmor_table"]

    filepaths = _get_filepaths(in_dir, cfg["filename"])

    if len(filepaths) > 0:
        logger.info("Found %d input files in '%s'", len(filepaths), in_dir)
    else:
        logger.info("No files found, basename: %s", cfg["filename"])

    for var, var_info in cfg["variables"].items():
        logger.info("CMORizing variable '%s'", var)
        glob_attrs["mip"] = var_info["mip"]
        cmor_info = cmor_table.get_variable(var_info["mip"], var)
        raw_var = var_info.get("raw", var)
        _extract_variable(raw_var, cmor_info, glob_attrs, filepaths, out_dir)
