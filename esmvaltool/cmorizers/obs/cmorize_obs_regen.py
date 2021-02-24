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

Masking Information
    Download the following files:
        REGEN_AllStns_{version}_1950-2016_QualityMask.nc

Monthly data
    Monthly data is additionally derived.
"""

import copy
import logging
import os
from pathlib import Path

import cftime
import iris
import numpy as np
from cf_units import Unit
from esmvalcore.preprocessor import (monthly_statistics, annual_statistics)

from . import utilities as utils

logger = logging.getLogger(__name__)

def _fix_time_coord(cube, var):
    """Correct wrong time points."""
    # Fix calendar type
    cube.coord('time').units = Unit(cube.coord('time').units.origin,
                                    calendar=var.get('calendar'))
    cube.coord('time').convert_units(
        Unit('days since 1950-1-1 00:00:00', calendar='gregorian'))

    # time points are 00:00:00, should be 12:00:00
    time = cube.coord('time')
    points = cube.coord('time').units.num2date(cube.coord('time').points)
    points = [
        cftime.DatetimeGregorian(c.year, c.month, c.day, 12) for c in points
    ]
    cube.coord('time').points = time.units.date2num(points)


def _extract_variable(short_name, qs, cfg, file_path, mask_path, out_dir):
    """Extract variable."""
    var = cfg['variables'].get(short_name)

    raw_var = var.get('raw')
    cube = iris.load_cube(file_path, utils.var_name_constraint(raw_var))

    # Fix units
    cmor_info = cfg['cmor_table'].get_variable(var['mip'], short_name)
    if 'raw_units' in var:
        cube.units = var['raw_units']
    cube.convert_units(cmor_info.units)

    # Fix time coord
    _fix_time_coord(cube, var)

    # Fix coordinates
    utils.fix_coords(cube)

    # Fix metadata
    attrs = copy.deepcopy(cfg['attributes'])
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
            cube_amon = monthly_statistics(cube)
            for cub in [cube, cube_amon]:
                cub.remove_coord('month_number')
                cub.remove_coord('year')

            # Fix metadata
            cmor_info_amon = cfg['cmor_table'].get_variable('Amon', short_name)
            attrs_amon = copy.deepcopy(attrs)
            attrs_amon['mip'] = 'Amon'
            utils.fix_var_metadata(cube_amon, cmor_info_amon)
            utils.set_global_atts(cube_amon, attrs_amon)

            # Save variable
            utils.save_variable(cube_amon,
                                short_name,
                                out_dir,
                                attrs_amon,
                                unlimited_dimensions=['time'])

    # masking reduced
    logger.info("Deriving long-term quality mask")
    mask_cube = iris.load_cube(mask_path)
    utils.fix_coords(mask_cube)
    mask = mask_cube.data.mask
    mask[np.where(mask_cube.data.data == 0)] = True
    mask_day = np.array([mask] * cube.shape[0])
    cube.data.mask = np.ma.mask_or(cube.data.mask, mask_day)

    ## Fix metadata
    cube.attributes.update({'masking': 'based on {}'.format(
        os.path.basename(mask_path))})
    attrs['version'] += '-masked'

    ## Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])

    if 'add_mon' in var.keys():
        if var['add_mon']:
            logger.info("Applying mask to monthly data")

            mask_amon = np.array([mask] * cube_amon.shape[0])
            cube_amon.data.mask = np.ma.mask_or(cube_amon.data.mask, mask_amon)

            cube_amon.attributes.update({'masking': 'based on {}'.format(
                os.path.basename(mask_path))})
            attrs_amon['version'] += '-masked'

            # Save variable
            utils.save_variable(cube_amon,
                                short_name,
                                out_dir,
                                attrs_amon,
                                unlimited_dimensions=['time'])


def cmorization(in_dir, out_dir, cfg, _):
    """Cmorization func call."""
    raw_filename = cfg['filename']
    file_names = raw_filename.format(version=cfg['attributes']['version'])
    mask_path = os.path.join(in_dir,
                             cfg['filename_mask'].format(
                                 version=cfg['attributes']['version'])
                            )

    # Run the cmorization
    for (short_name, var) in cfg['variables'].items():
        logger.info("CMORizing variable '%s'", short_name)
        file_paths = [str(fp) for fp in sorted(Path(in_dir).glob(file_names))]
        file_paths.remove(os.path.join(in_dir, mask_path))

        for file_path in file_paths:
            logger.info("Loading '%s'", file_path)
            _extract_variable(short_name, var, cfg, file_path, mask_path,
                              out_dir)
