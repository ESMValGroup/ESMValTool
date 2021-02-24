"""ESMValTool CMORizer for GPCP data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://psl.noaa.gov/data/gridded/data.gpcp.html

Last access
    20210224

Download and processing instructions
    Download the following file:
        ftp://ftp.cdc.noaa.gov/Datasets/gpcp/precip.mon.mean.nc
"""

import logging
import os

import cftime
import iris
import numpy as np

from . import utilities as utils

logger = logging.getLogger(__name__)


def _get_centered_timecoord(cube):
    """Fix time coordinate.

    Time points start at the beginning of month at 00:00:00.
    """
    time = cube.coord('time')
    times = time.units.num2date(time.points)

    # get bounds
    starts = [
        cftime.DatetimeNoLeap(c.year, c.month, 1)
        for c in times
    ]
    ends = [
        cftime.DatetimeNoLeap(c.year, c.month + 1, 1) if c.month < 12
        else cftime.DatetimeNoLeap(c.year + 1, 1, 1)
        for c in times
    ]
    time.bounds = time.units.date2num(np.stack([starts, ends], -1))

    # get points
    time.points = [np.mean((t1, t2)) for t1, t2 in time.bounds]


def _cmorize_variable(short_name, var, cfg, filepath, out_dir):
    """Cmorize variable."""
    raw_var = var.get('raw', short_name)
    cube = iris.load_cube(filepath, utils.var_name_constraint(raw_var))

    # Fix units
    if 'raw_units' in var:
        cube.units = var['raw_units']
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    cube.convert_units(cmor_info.units)
    utils.convert_timeunits(cube, 1950)

    # Fix coordinates
    # fix time
    _get_centered_timecoord(cube)

    # Fix coordinates
    utils.fix_coords(cube, overwrite_time_bounds=False)

    # Fix metadata
    attrs = cfg['attributes']
    attrs['mip'] = var['mip']
    utils.fix_var_metadata(cube, cmor_info)
    utils.set_global_atts(cube, attrs)

    # Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    raw_filepath = os.path.join(in_dir, cfg['filename'])

    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", short_name)
        _cmorize_variable(short_name, var, cfg, raw_filepath, out_dir)
