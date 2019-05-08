import logging
import os

import iris
from iris.coord_categorisation import add_month, add_year
from iris.analysis import MEAN

from .utilities import (_set_global_atts, _convert_timeunits, _fix_coords,
                        _fix_var_metadata, _save_variable)

logger = logging.getLogger(__name__)


def cmorize_osi(in_dir, out_dir, cfg, hemisphere):
    cmor_table = cfg['cmor_table']
    glob_attrs = cfg['attributes']

    logger.info("Starting cmorization for Tier%s OBS files: %s",
                glob_attrs['tier'], glob_attrs['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)

    # run the cmorization
    for var, vals in cfg['variables'].items():
        for year in os.listdir(in_dir):
            logger.info(
                "CMORizing var %s:%s for year %s", vals['mip'], var, year
            )
            var_info = cmor_table.get_variable(vals['mip'], var)
            file_pattern = '{0}_{1}_{2}_*.nc'.format(
                vals['raw'], hemisphere, vals['grid']
            )
            raw_info = {
                'name': vals['raw'],
                'file': os.path.join(in_dir, year, '??', file_pattern)
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
        # OSI-409 and OSI-450 do not have the same attributes
        try:
            tracking_ids.append(cube.attributes['tracking_id'])
        except KeyError:
            pass

        to_remove = [
            'time_coverage_start', 'time_coverage_end',
            'history', 'tracking_id', 'start_date', 'stop_date',
        ]
        for attr in to_remove:
            try:
                del cube.attributes[attr]
            except KeyError:
                pass
    cube = cubes.concatenate_cube()
    del cubes
    if var_info.frequency == 'mon':
        add_month(cube, 'time')
        add_year(cube, 'time')
        cube = cube.aggregated_by(('month', 'year'), MEAN)
        cube.remove_coord('month')
        cube.remove_coord('year')
    if var_info.frequency == 'day':
        if cube.coord('time').shape[0] < 300:
            logger.warning(
                'Only %s days available. Can not generate daily files',
                cube.coord('time').shape[0]
            )
            return
    logger.debug(cube)
    if tracking_ids:
        cube.attributes['tracking_ids'] = tracking_ids
    cube.coord('projection_x_coordinate').var_name = 'x'
    cube.coord('projection_y_coordinate').var_name = 'y'

    _fix_var_metadata(cube, var_info)
    _convert_timeunits(cube, year)
    _fix_coords(cube)
    _set_global_atts(cube, attrs)
    _save_variable(cube, var, out_dir, attrs)
    del cube
