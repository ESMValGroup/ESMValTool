"""
ESMValTool CMORizer for JRA-3Q data.

Tier
    Tier 2: other freely-available dataset.

Source
    Research Data Archive (RDA):
    https://rda.ucar.edu/datasets/d640000/


Last access
    20251014



"""

import copy
import logging
import os

import iris
from esmvalcore.preprocessor import daily_statistics
from iris import NameConstraint

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _extract_variable(short_name, var, filename, cfg, in_dir, out_dir):
    cmor_info = cfg["cmor_table"].get_variable(var["mip"], short_name)
    attrs = copy.deepcopy(cfg["attributes"])
    attrs["mip"] = var["mip"]
    ver = attrs["version"]
    freq = var["freq"]
    raw_var = var.get("raw_short_name", short_name)

    # load data
    filepath = os.path.join(in_dir, ver, freq, raw_var, filename)
    cubes = iris.load(filepath, NameConstraint(var_name=raw_var))

    for cube in cubes:
        # JRA3Q has hourly max/min temperature data, but we want daily max/min
        if short_name == "tasmax":
            cube = daily_statistics(cube, "max")
        elif short_name == "tasmin":
            cube = daily_statistics(cube, "min")

        # Fix metadata

        utils.fix_var_metadata(cube, cmor_info)
        utils.set_global_atts(cube, attrs)

        utils.fix_dim_coordnames(cube)
        utils.fix_coords(cube)
        if "height2m" in cmor_info.dimensions:
            utils.add_height2m(cube)
        utils.set_global_atts(cube, attrs)

        # Save variable
        utils.save_variable(
            cube, short_name, out_dir, attrs, unlimited_dimensions=["time"]
        )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    # Run the cmorization
    for short_name, var in cfg["variables"].items():
        filename = var["file"]
        logger.info(
            "CMORizing variable '%s' from file '%s'", short_name, filename
        )
        _extract_variable(short_name, var, filename, cfg, in_dir, out_dir)
