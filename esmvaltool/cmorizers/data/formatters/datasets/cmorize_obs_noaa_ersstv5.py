"""ESMValTool CMORizer for NOAA Extended Reconstructed Sea Surface
   Temperature (ERSST), Version 5 data.

Tier
    Tier 2: open dataset.

Source
    https://doi.org/10.7289/V5T72FNM

Last access
    20200520

Download and processing instructions
    TODO

"""

import logging
import os

import iris
from cf_units import Unit

from . import utilities as utils

logger = logging.getLogger(__name__)


def _fix_time_coord(cube):
    """Set time points to central day of month."""
    time_coord = cube.coord('time')
    new_unit = Unit('days since 1850-01-01 00:00:00', calendar='standard')
    time_coord.convert_units(new_unit)
    old_time = new_unit.num2date(time_coord.points)
    new_time = [d.replace(day=15) for d in old_time]
    time_coord.points = new_unit.date2num(new_time)


def _extract_variable(raw_var, cmor_info, attrs, filepath, out_dir):
    """Extract variable."""
    var = cmor_info.short_name
    cube = iris.load_cube(filepath, utils.var_name_constraint(raw_var))
    cube = iris.util.squeeze(cube)
    #_fix_time_coord(cube)
    utils.fix_var_metadata(cube, cmor_info)
    #utils.fix_coords(cube)
    utils.set_global_atts(cube, attrs)
    # utils.flip_dim_coord(cube, 'latitude')
    utils.save_variable(cube,
                        var,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


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
