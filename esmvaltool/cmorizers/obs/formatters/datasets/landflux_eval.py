"""ESMValTool CMORizer for LandFlux-EVAL data.

Tier
    Tier 3: restricted dataset.

Source
    https://data.iac.ethz.ch/landflux/

Last access
    20190516

Download and processing instructions
    Download the following files:
        LandFluxEVAL.merged.89-05.monthly.all.nc
    A registration is required for downloading the data (see
    <http://www.iac.ethz.ch/group/land-climate-dynamics/research/
    landflux-eval.html>).

"""

import logging
import os
from datetime import datetime

import iris
import numpy as np
from cf_units import Unit

from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)


def _extract_variable(raw_var, cmor_info, attrs, filepath, out_dir):
    """Extract variable."""
    var = cmor_info.short_name
    cube = iris.load_cube(filepath, utils.var_name_constraint(raw_var))
    _fix_time_coord(cube)
    utils.fix_var_metadata(cube, cmor_info)
    utils.convert_timeunits(cube, 1950)
    utils.fix_coords(cube)
    utils.set_global_atts(cube, attrs)
    utils.save_variable(
        cube, var, out_dir, attrs, unlimited_dimensions=['time'])


def _fix_time_coord(cube):
    """Fix time coordinate (given as month as %Y%m.%f)."""
    time_coord = cube.coord('time')
    new_units = Unit('days since 1950-1-1 00:00:00', calendar='standard')

    # Function to convert given date to correct number
    def _date2num(date_str):
        """Convert data given as %Y%m.%f to number."""
        date_str = str(date_str)
        year = int(date_str[:4])
        month = int(date_str[4:6])
        day = 15
        date = datetime(year, month, day)
        return new_units.date2num(date)

    # Convert time coordinate array and set correct units
    time_coord.points = np.vectorize(_date2num)(time_coord.points)
    time_coord.units = new_units
    time_coord.attributes = {}


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']
    filepath = os.path.join(in_dir, cfg['filename'])
    logger.info("Found input file '%s'", filepath)

    # Run the cmorization
    for (var, var_info) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", var)
        glob_attrs['mip'] = var_info['mip']
        cmor_info = cmor_table.get_variable(var_info['mip'], var)
        raw_var = var_info.get('raw', var)
        _extract_variable(raw_var, cmor_info, glob_attrs, filepath, out_dir)
