"""ESMValTool CMORizer for CRU data.

Tier
    Tier 2: other freely-available dataset.

Source
    TS4.02: https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.02/cruts.1811131722.v4.02/  # noqa: E501
    TS4.06: https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.06/cruts.2205201912.v4.06/  # noqa: E501
    TS4.07: https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.07/cruts.2304141047.v4.07/  # noqa: E501

Last access
    TS4.02: 20190516
    TS4.06: 20231012
    TS4.07: 20231012

Download and processing instructions
    Download the following files:
    ``{raw_name}/cru_ts4.{X}.1901.{end_year}.{raw_name}.dat.nc.gz``
    where ``{raw_name}`` is the name of the desired variable(s) or run
    ``esmvaltool data download CRU`` for the latest version
"""

import logging
import os

import cftime
import iris
import numpy as np
from cf_units import Unit
from iris import NameConstraint

from esmvaltool.cmorizers.data import utilities as utils

logger = logging.getLogger(__name__)


def _center_timecoord(cube):
    """Set time coordinates to exact center of each month.

    CRU timepoints are not in the center of the month and added bounds
    by utils.fix_coords are incorrect. #1981
    """
    time = cube.coord("time")
    times = time.units.num2date(time.points)

    # get bounds
    starts = [cftime.DatetimeNoLeap(c.year, c.month, 1) for c in times]
    ends = [
        cftime.DatetimeNoLeap(c.year, c.month + 1, 1)
        if c.month < 12 else cftime.DatetimeNoLeap(c.year + 1, 1, 1)
        for c in times
    ]
    time.bounds = time.units.date2num(np.stack([starts, ends], -1))
    time.points = [np.mean((t1, t2)) for t1, t2 in time.bounds]


def _extract_variable(short_name, var, cfg, filepath, out_dir):
    """Extract variable."""
    raw_var = var.get("raw", short_name)
    version = cfg["attributes"]["version"]
    cube = iris.load_cube(filepath, NameConstraint(var_name=raw_var))

    # Fix units
    if "raw_units" in var:
        cube.units = var["raw_units"]
    cmor_info = cfg["cmor_table"].get_variable(var["mip"], short_name)
    cube.convert_units(cmor_info.units)
    if version in ["TS4.02"]:
        utils.convert_timeunits(cube, 1950)
    else:
        cube.coord("time").convert_units(
            Unit("days since 1950-1-1 00:00:00", calendar="gregorian"))

    # Fix coordinates
    cube = utils.fix_coords(cube)
    if "height2m" in cmor_info.dimensions:
        utils.add_height2m(cube)
    if version not in ["TS4.02"]:
        _center_timecoord(cube)

    # Fix metadata
    attrs = cfg["attributes"]
    attrs["mip"] = var["mip"]
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=["time"])


def cmorization(in_dir, out_dir, cfg, cfg_user, start_date, end_date):
    """Cmorization func call."""
    raw_filepath = os.path.join(in_dir, cfg["filename"])

    # Run the cmorization
    for short_name, var in cfg["variables"].items():
        logger.info("CMORizing variable '%s'", short_name)
        raw_var = var.get("raw", short_name)
        filepath = raw_filepath.format(raw_name=raw_var)
        if filepath is None:
            continue
        _extract_variable(short_name, var, cfg, filepath, out_dir)
