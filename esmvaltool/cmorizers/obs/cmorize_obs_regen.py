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
    reducedmasked:
    - kriging error within the 95th percentile of all the grids on the day
    - coefficient of variation within the 95th percentile of all the grids on
        the day

    masked:
    - kriging error within the 95th percentile of all the grids on the day
    - coefficient of variation within the 95th percentile of all the grids on
        the day
    - at least one station per grid point

Monthly data
    Monthly data is additionally derived. For both masking cases, grid point
    data is set to missing if for any given month <80% is valid.
"""

import copy
import logging
from pathlib import Path

import cftime
import iris
import numpy as np
from cf_units import Unit
from esmvalcore.preprocessor import monthly_statistics

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


def _get_mask_qbased(cube, cfg, var_str):
    """Get percentile based mask."""
    q_target = cfg['masking_thresholds'].get(var_str)
    qs = [np.percentile(data.compressed(), q_target) for data in cube.data]
    # for data, q, m in zip(cube.data, qs, mask):
    #     m[np.where(data >= q)] = 1
    mask = np.zeros(cube.shape)
    for data, q, m in zip(cube.data, qs, mask):
        m[np.where(data > q)] = 1
    return mask.astype('bool')


def _get_mask_sbased(cube, cfg, var_str):
    """Get station based mask."""
    s_target = cfg['masking_thresholds'].get(var_str)

    mask = np.zeros(cube.shape)
    for data, m in zip(cube.data, mask):
        m[np.where(data < s_target)] = 1
    return mask.astype('bool')


def _build_monthy_masked_cube(cube, cfg):
    """Build mask based on valid data per month (80%)."""
    # Calc monthly
    cube_amon = copy.deepcopy(cube)
    iris.coord_categorisation.add_month_number(cube_amon,
                                               'time',
                                               name='month_number')
    monthly_masks = []
    months = list(set(cube_amon.coord('month_number').points))
    m_threshhold = cfg['masking_thresholds'].get('monthly_statisitcs')
    for mon_ind, mon in enumerate(months):
        const = iris.Constraint(month_number=mon)
        mon_cube = cube_amon.extract(const)
        mon_mask = mon_cube.data.mask.sum(axis=0)
        mon_mask = np.ma.masked_where((mon_mask / mon_cube.shape[0]) >= \
            (1 - m_threshhold), mon_mask)
        monthly_masks.append(mon_mask.mask)
    monthly_masks = np.array(monthly_masks)

    # Calc monthly
    cube_amon = monthly_statistics(cube_amon)
    cube_amon.remove_coord('month_number')
    cube_amon.remove_coord('year')
    cube_amon.data.mask = monthly_masks

    return cube_amon


def _extract_variable(short_name, var, cfg, file_path, out_dir):
    """Extract variable."""

    raw_var = var.get('raw', short_name)
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
    logger.info("Deriving mask based on krigging error and coefficient of "\
                "variation")

    ## krigging error
    raw_var = cfg['masking_variables_raw']['ek'].get('raw', short_name)
    ek_cube = iris.load_cube(file_path, utils.var_name_constraint(raw_var))
    _fix_time_coord(ek_cube, cfg['masking_variables_raw']['ek'])
    utils.fix_coords(ek_cube)

    ek_mask = _get_mask_qbased(ek_cube, cfg, 'ek')

    ## coefficient of variation
    raw_var = cfg['masking_variables_raw']['sd'].get('raw', short_name)
    sd_cube = iris.load_cube(file_path, utils.var_name_constraint(raw_var))
    _fix_time_coord(sd_cube, cfg['masking_variables_raw']['sd'])
    utils.fix_coords(sd_cube)

    sd_cube.units = cfg['masking_variables_raw']['sd']['raw_units']
    sd_cube.convert_units(cmor_info.units)
    cov_mask = _get_mask_qbased(sd_cube / cube, cfg, 'cov')

    ## resulting mask
    res_mask = np.ma.mask_or(ek_mask, cov_mask)
    cube.data.mask = np.ma.mask_or(cube.data.mask, res_mask)

    ## Fix metadata
    attrs = copy.deepcopy(cfg['attributes_reducedmasked'])
    attrs['mip'] = var['mip']
    utils.set_global_atts(cube, attrs)

    ## Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])

    if 'add_mon' in var.keys():
        if var['add_mon']:
            logger.info("Building monthly means")

            # get monthly masked data based on monthly threshold
            cube_amon = _build_monthy_masked_cube(cube, cfg)

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

    # masking
    logger.info("Deriving mask based on krigging error, coefficient of "\
                "variation and station number")

    ## station number
    raw_var = cfg['masking_variables_raw']['s'].get('raw', short_name)
    s_cube = iris.load_cube(file_path, utils.var_name_constraint(raw_var))
    _fix_time_coord(s_cube, cfg['masking_variables_raw']['s'])
    utils.fix_coords(s_cube)

    s_mask = _get_mask_sbased(s_cube, cfg, 's')
    cube.data.mask = np.ma.mask_or(cube.data.mask, s_mask)

    ## Fix metadata
    attrs = copy.deepcopy(cfg['attributes_masked'])
    attrs['mip'] = var['mip']
    utils.set_global_atts(cube, attrs)

    ## Save variable
    utils.save_variable(cube,
                        short_name,
                        out_dir,
                        attrs,
                        unlimited_dimensions=['time'])

    if 'add_mon' in var.keys():
        if var['add_mon']:
            logger.info("Building monthly means")

            # get monthly masked data based on monthly threshold
            cube_amon = _build_monthy_masked_cube(cube, cfg)

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
