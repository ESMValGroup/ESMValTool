"""ESMValTool CMORizer for HWSD data.

Tier
    Tier 3: restricted dataset.

Source
    https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1247

Last access
    20191015

Download and processing instructions
    Download the following file:
        HWSD_SOIL_CLM_RES.nc4
    A registration is required for downloading the data.

"""

import logging
import os

import iris
from cf_units import Unit

from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)


def _extract_variable(short_name, var, cfg, filepath, out_dir):
    """Extract variable."""
    raw_var = var.get('raw', short_name)
    cube = iris.load_cube(filepath, utils.var_name_constraint(raw_var))

    # Sum over levels
    if short_name in ('cSoil', ):
        level_coord = iris.coords.DimCoord([0, 1], long_name='level')
        cube.add_dim_coord(level_coord, 0)
        cube = cube.collapsed('level', iris.analysis.SUM)

    # Fix coordinates
    if var['mip'] != 'fx':
        cube = iris.util.new_axis(cube)
        time_dim = iris.coords.DimCoord(
            [183.0],
            bounds=[0.0, 366.0],
            units=Unit('days since 2000-01-01 00:00:00'),
            standard_name='time',
            var_name='time',
            long_name='time')
        cube.add_dim_coord(time_dim, 0)
        utils.convert_timeunits(cube, 1950)
    utils.fix_coords(cube)

    # Fix units
    if 'kg C' in cube.units.origin:
        cube.units = Unit(cube.units.origin.replace('C', ''))
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    cube.convert_units(cmor_info.units)

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
    filepath = os.path.join(in_dir, cfg['filename'])
    logger.info("Reading file '%s'", filepath)

    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", short_name)
        _extract_variable(short_name, var, cfg, filepath, out_dir)
