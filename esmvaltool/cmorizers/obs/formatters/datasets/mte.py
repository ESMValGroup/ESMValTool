"""ESMValTool CMORizer for MTE data.

Tier
    Tier 3: restricted dataset.

Source
    http://www.bgc-jena.mpg.de/geodb/BGI/Home

Last access
    20190507

Download and processing instructions
    Download the following files:
        EnsembleGPP_GL.nc
    A registration is required for downloading the data.

"""

import logging
import os

import iris
from cf_units import Unit

from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)


def _fix_time_coord(cube):
    """Set time points to central day of month."""
    time_coord = cube.coord('time')
    new_unit = Unit('days since 1850-01-01 00:00:00', calendar='standard')
    time_coord.convert_units(new_unit)
    old_time = new_unit.num2date(time_coord.points)
    new_time = [d.replace(day=15) for d in old_time]
    time_coord.points = new_unit.date2num(new_time)


def _get_filepath(in_dir, basename):
    """Find correct name of file (extend basename with timestamp)."""
    all_files = [
        f for f in os.listdir(in_dir)
        if os.path.isfile(os.path.join(in_dir, f))
    ]
    for filename in all_files:
        if filename.endswith(basename):
            return os.path.join(in_dir, filename)
    raise OSError(
        f"Cannot find input file ending with '{basename}' in '{in_dir}'")


def _extract_variable(raw_var, cmor_info, attrs, filepath, out_dir):
    """Extract variable."""
    var = cmor_info.short_name
    cube = iris.load_cube(filepath, utils.var_name_constraint(raw_var))
    _fix_time_coord(cube)
    utils.fix_var_metadata(cube, cmor_info)
    utils.convert_timeunits(cube, 1950)
    utils.fix_coords(cube)
    utils.set_global_atts(cube, attrs)
    utils.flip_dim_coord(cube, 'latitude')
    utils.save_variable(cube,
                        var,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    glob_attrs = cfg['attributes']
    cmor_table = cfg['cmor_table']
    filepath = _get_filepath(in_dir, cfg['filename'])
    logger.info("Found input file '%s'", filepath)

    # Run the cmorization
    for (var, var_info) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", var)
        glob_attrs['mip'] = var_info['mip']
        cmor_info = cmor_table.get_variable(var_info['mip'], var)
        raw_var = var_info.get('raw', var)
        _extract_variable(raw_var, cmor_info, glob_attrs, filepath, out_dir)
