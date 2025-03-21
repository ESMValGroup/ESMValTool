"""ESMValTool CMORizer for GISTEMP data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://data.giss.nasa.gov/gistemp/
    https://data.giss.nasa.gov/pub/gistemp/gistemp250_GHCNv4.nc.gz

Last access
    20200303
"""

import logging
import os

import iris
from iris import NameConstraint

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _extract_variable(short_name, var, cfg, filepath, out_dir):
    """Extract variable."""
    raw_var = var.get("raw", short_name)
    cube = iris.load_cube(filepath, NameConstraint(var_name=raw_var))

    # Fix units
    cmor_info = cfg["cmor_table"].get_variable(var["mip"], short_name)
    cube.convert_units(cmor_info.units)
    utils.convert_timeunits(cube, 1950)

    # Fix coordinates
    cube = utils.fix_coords(cube)
    if "height2m" in cmor_info.dimensions:
        utils.add_height2m(cube)

    # Fix metadata
    attrs = cfg["attributes"]
    attrs["mip"] = var["mip"]
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(
        cube, short_name, out_dir, attrs, unlimited_dimensions=["time"]
    )


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    filepath = os.path.join(in_dir, cfg["filename"])

    # Run the cmorization
    for short_name, var in cfg["variables"].items():
        logger.info("CMORizing variable '%s'", short_name)
        _extract_variable(short_name, var, cfg, filepath, out_dir)
