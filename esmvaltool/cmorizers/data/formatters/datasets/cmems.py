"""ESMValTool CMORizer for CMEMS - Sea Level Thematic Assembly Center.

   This is the CMORizer script for Sea Surface Height Above Geoid.
   Copernicus Marine Environment Monitoring Service (CMEMS)
   Product: SEALEVEL_GLO_PHY_L4_MY_008_047

Tier
    Tier 2: open dataset.

Source
    https://doi.org/10.48670/moi-00148

Last access
    20250220

Download and processing instructions
    Download daily files from:
    https://data.marine.copernicus.eu/product/SEALEVEL_GLO_PHY_L4_MY_008_047/services
    Using daily data (sea surface height above geoid) as the month product is sea level anomalies.

"""

import logging
import os
import re
from pathlib import Path

import iris
from esmvalcore.preprocessor import monthly_statistics

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _get_filepaths(in_dir, basename):
    """Group daily files."""
    regex = re.compile(basename)
    return_ls = []

    for root, _dir, files in os.walk(in_dir, followlinks=True):
        if len(files) > 1:
            return_files = [
                (Path(root) / filename)
                for filename in files
                if regex.match(filename)
            ]

            return_ls.append(return_files)

    return return_ls


def _extract_variable(cmor_info, attrs, file_ls, out_dir):
    """Extract variable and aggregate months."""
    var = cmor_info.short_name
    standard_name = cmor_info.standard_name

    cube_prepls = iris.cube.CubeList()
    for filepaths in file_ls:
        cubels = iris.load(filepaths, standard_name)
        iris.util.equalise_attributes(cubels)
        iris.util.unify_time_units(cubels)
        # if cmor.frequency is monthly
        cube_mon = cubels.concatenate_cube()
        cube_mon = iris.util.squeeze(cube_mon)
        cube_prepls.append(monthly_statistics(cube_mon, "mean"))

    iris.util.equalise_attributes(cube_prepls)
    cube = cube_prepls.concatenate_cube()
    utils.fix_var_metadata(cube, cmor_info)
    cube = utils.fix_coords(cube)

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

    file_ls = _get_filepaths(in_dir, cfg["filename"])

    if len(file_ls) > 0:
        logger.info("Found %d input files in '%s'", len(file_ls), in_dir)
    else:
        logger.info("No files found, basename: %s", cfg["filename"])

    # Run the cmorization
    for var, var_info in cfg["variables"].items():
        logger.info("CMORizing variable '%s'", var)
        glob_attrs["mip"] = var_info["mip"]
        cmor_info = cmor_table.get_variable(var_info["mip"], var)

        _extract_variable(cmor_info, glob_attrs, file_ls, out_dir)
