import logging
import os

import iris
from iris.coord_categorisation import add_month, add_year
from iris.analysis import MEAN

from .utilities import (_set_global_atts, _convert_timeunits, _fix_coords,
                        _fix_var_metadata, _save_variable)

logger = logging.getLogger(__name__)

def cmorize_osi(in_dir, out_dir, cfg,  hemisphere):
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

    logger.info("Starting cmorization for Tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)

    # run the cmorization
    for var, vals in cfg['variables'].items():
        for year in os.listdir(in_dir):
            logger.info("CMORizing var %s for year %s", var, year)
            var_info = cmor_table.get_variable(vals['mip'], var)
            raw_info = {
                'name': vals['raw'],
                'file': os.path.join(
                    in_dir, year, '??',
                    '{0}_{1}_*.nc'.format(vals['raw'], hemisphere)
                )
            }
            glob_attrs['mip'] = vals['mip']
            extract_variable(var_info, raw_info, out_dir, glob_attrs, year)

def extract_variable(var_info, raw_info, out_dir, attrs, year):
    """Extract to all vars."""
    var = var_info.short_name
    cubes = iris.load(
        raw_info['file'],
        iris.Constraint(cube_func=lambda c: c.var_name == raw_info['name'])
    )

    tracking_ids = []
    for cube in cubes:
        tracking_ids.append(cube.attributes['tracking_id'])
        for attr in ['time_coverage_start', 'time_coverage_end', 'history', 'tracking_id']:
            del cube.attributes[attr]
    cube = cubes.concatenate_cube()
    add_month(cube, 'time')
    add_year(cube, 'time')
    cube.aggregated_by(('month', 'year'), MEAN)
    logger.debug(cube)
    cube.coord('projection_x_coordinate').var_name = 'x'
    cube.coord('projection_y_coordinate').var_name = 'y'
    cube.remove_coord('month')
    cube.remove_coord('year')
    _fix_var_metadata(cube, var_info)
    _convert_timeunits(cube, year)
    _fix_coords(cube)
    _set_global_atts(cube, attrs)
    _save_variable(cube, var, out_dir, attrs)