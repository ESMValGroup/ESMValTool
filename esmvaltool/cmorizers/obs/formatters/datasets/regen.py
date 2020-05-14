"""ESMValTool CMORizer for REGEN data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://researchdata.ands.org.au/rainfall-estimates-gridded-v1-2019/1408744
Last access
    20200226

Download and processing instructions
    Download the following files:
        REGEN_AllStns_{version}_[1950..2016].nc

"""

import logging
from pathlib import Path

import cf_units
import iris

from esmvalcore.preprocessor import monthly_statistics
from esmvaltool.cmorizers.obs import utilities as utils

logger = logging.getLogger(__name__)


def _extract_variable(short_name, var, cfg, file_path, out_dir):
    """Extract variable."""

    raw_var = var.get('raw', short_name)
    cube = iris.load_cube(file_path, utils.var_name_constraint(raw_var))

    # Fix units
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    if 'raw_units' in var:
        cube.units = var['raw_units']
    cube.convert_units(cmor_info.units)

    # Fix calendar type
    cal_time = var.get('calendar', short_name)
    origin_time = cube.coord('time').units.origin
    cube.coord('time').units = cf_units.Unit(origin_time, calendar=cal_time)
    utils.convert_timeunits(cube, 1950)

    # Fix coordinates
    utils.fix_coords(cube)

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

    if 'add_mon' in var.keys():
        if var['add_mon']:
            logger.info("Building monthly means")

            # Calc monthly
            cube = monthly_statistics(cube)
            cube.remove_coord('month_number')
            cube.remove_coord('year')

            # Fix metadata
            attrs['mip'] = 'Amon'

            # Fix coordinates
            utils.fix_coords(cube)

            # Save variable
            utils.save_variable(cube,
                                short_name,
                                out_dir,
                                attrs,
                                unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    raw_filename = cfg['filename']
    file_names = raw_filename.format(version=cfg['attributes']['version'])

    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", short_name)
        for file_path in sorted(Path(in_dir).glob(file_names)):
            logger.info("Loading '%s'", file_path)
            _extract_variable(short_name, var, cfg, str(file_path), out_dir)
