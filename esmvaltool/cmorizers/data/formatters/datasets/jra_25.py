"""
ESMValTool CMORizer for JRA-25 data.

Tier
    Tier 2: other freely-available dataset.

Source
    ESGF:
    https://esgf.nccs.nasa.gov/thredds/fileServer/CREATE-IP/
        reanalysis/JMA/JRA-25/JRA-25/mon/atmos/

Last access
    20221122

Download and processing instructions
    see download script cmorizers/data/downloaders/datasets/jra_25.py
"""

import copy
import logging
import os

import iris
from iris import NameConstraint

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _extract_variable(short_name, var, filename, cfg, in_dir, out_dir):
    """Extract variable."""
    # load data
    filepath = os.path.join(in_dir, filename)
    raw_var = var.get("raw", short_name)
    cube = iris.load_cube(filepath, NameConstraint(var_name=raw_var))

    # Fix metadata
    cmor_info = cfg["cmor_table"].get_variable(var["mip"], short_name)
    attrs = copy.deepcopy(cfg["attributes"])
    attrs["mip"] = var["mip"]
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(
        cube, short_name, out_dir, attrs, unlimited_dimensions=["time"]
    )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    # Run the cmorization
    for short_name, var in cfg["variables"].items():
        short_name = var["short_name"]
        filename = var["file"]
        logger.info(
            "CMORizing variable '%s' from file '%s'", short_name, filename
        )
        _extract_variable(short_name, var, filename, cfg, in_dir, out_dir)
