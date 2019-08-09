"""Common functionalities for OSI-450 and OSI-407 dataset cmorization."""

import logging
import os
import glob
from datetime import datetime, timedelta
from calendar import monthrange, isleap

import numpy as np
from numba import vectorize
import iris
from iris.cube import Cube, CubeList
from iris.coords import AuxCoord
from iris.coord_categorisation import (
    add_month_number, add_year, add_day_of_year)
from iris.analysis import MEAN

from .utilities import (set_global_atts, convert_timeunits, fix_var_metadata,
                        save_variable)

logger = logging.getLogger(__name__)


def cmorize_osi(in_dir, out_dir, cfg, hemisphere):
    """Cmorize OSI-450 or OSI-407 dataset"""
    logger.info("Starting cmorization for Tier%s OBS files: %s",
                cfg['attributes']['tier'], cfg['attributes']['dataset_id'])
    logger.info("Input data from: %s", in_dir)
    logger.info("Output will be written to: %s", out_dir)

    # run the cmorization
    first_run = True
    for var, vals in cfg['variables'].items():
        var_info = {}
        for mip in vals['mip']:
            var_info[mip] = cfg['cmor_table'].get_variable(mip, var)
        file_pattern = '{0}_{1}_{2}_*.nc'.format(
            vals['raw'], hemisphere, vals['grid']
        )
        for year in os.listdir(in_dir):
            year = int(year)
            logger.info(
                "CMORizing var %s for year %s", var, year
            )
            raw_info = {
                'name': vals['raw'],
                'file': os.path.join(in_dir, str(year), '??', file_pattern)
            }
            extract_variable(var_info, raw_info, out_dir, cfg['attributes'],
                             year, vals['mip'])
            if first_run:
                sample_file = glob.glob(os.path.join(
                    in_dir, str(year), '01', file_pattern))[0]
                cube = iris.load_cube(
                    sample_file,
                    iris.Constraint(
                        cube_func=lambda c: c.var_name == raw_info['name'])
                )
                _create_areacello(cfg, cube, cfg['attributes'], out_dir)
                first_run = False


def extract_variable(var_infos, raw_info, out_dir, attrs, year, mips):
    """Extract to all vars."""
    cubes = iris.load(
        raw_info['file'],
        iris.Constraint(cube_func=lambda c: c.var_name == raw_info['name'])
    )
    tracking_ids = _unify_attributes(cubes)
    cube = cubes.concatenate_cube()
    del cubes
    if tracking_ids:
        cube.attributes['tracking_ids'] = tracking_ids
    cube.coord('projection_x_coordinate').var_name = 'x'
    cube.coord('projection_y_coordinate').var_name = 'y'
    lon_coord = cube.coord('longitude')
    lon_coord.points = _correct_lons(lon_coord.points)

    for mip in mips:
        var_info = var_infos[mip]
        attrs['mip'] = mip
        if var_info.frequency == 'mon':
            cube = _monthly_mean(cube)
            if cube.coord('time').shape[0] < 12:
                cubes = CubeList(cube.slices_over('time'))
                model_cube = cubes[0].copy()
                for month in range(1, 13):
                    month_constraint = iris.Constraint(
                        time=lambda cell: cell.point.month == month
                    )
                    if cubes.extract(month_constraint):
                        continue
                    cubes.append(_create_nan_cube(
                        model_cube, month, month=True))
                cube = cubes.merge_cube()
                del cubes
                del model_cube
        elif var_info.frequency == 'day':
            if cube.coord('time').shape[0] < 300:
                logger.warning(
                    'Only %s days available. Can not generate daily files',
                    cube.coord('time').shape[0]
                )
                continue
            total_days = 366 if isleap(year) else 365
            if cube.coord('time').shape[0] < total_days:
                cubes = add_nan_timesteps(cube, total_days, cubes)
                cube = cubes.merge_cube()
                cube.remove_coord('day_of_year')
                del cubes

        logger.debug(cube)
        fix_var_metadata(cube, var_info)
        convert_timeunits(cube, year)
        set_global_atts(cube, attrs)
        save_variable(cube, var_info.short_name, out_dir, attrs)
    return cube


def add_nan_timesteps(cube, total_days, cubes):
    add_day_of_year(cube, 'time')
    cubes = CubeList(cube.slices_over('time'))
    model_cube = cubes[0].copy()
    model_cube.remove_coord('day_of_year')
    for day_of_year in range(total_days):
        day_constraint = iris.Constraint(
            day_of_year=day_of_year + 1
        )
        if cubes.extract(day_constraint):
            continue
        nan_cube = _create_nan_cube(
            model_cube, day_of_year, month=False
        )
        add_day_of_year(nan_cube, 'time')
        cubes.append(nan_cube)
    del model_cube
    return cubes


def _create_nan_cube(model_cube, num, month):
    nan_cube = model_cube.copy(
        np.ma.masked_all(model_cube.shape)
    )
    time_coord = nan_cube.coord('time')
    nan_cube.remove_coord(time_coord)
    date = time_coord.cell(0).point
    if month:
        date = datetime(date.year, num, date.day)
        bounds = (
            datetime(date.year, num, 1),
            datetime(date.year, num, monthrange(date.year, num)[1])
        )
    else:
        date = datetime(date.year, 1, 1, 12) + timedelta(days=num)
        bounds = (
            datetime(date.year, 1, 1) + timedelta(days=num),
            datetime(date.year, 1, 1, 23, 59) + timedelta(days=num)
        )

    date = time_coord.units.date2num(date)
    bounds = (
        time_coord.units.date2num(bounds[0]),
        time_coord.units.date2num(bounds[1]),
    )
    nan_cube.add_aux_coord(AuxCoord(
        [date],
        standard_name=time_coord.standard_name,
        var_name=time_coord.var_name,
        long_name=time_coord.long_name,
        units=time_coord.units,
        attributes=time_coord.attributes,
        bounds=[bounds],
    ))
    return nan_cube


@vectorize(['float32(float32)'])
def _correct_lons(lon):
    if lon < 0:
        return lon + 360
    return lon


def _unify_attributes(cubes):
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
    return tracking_ids


def _monthly_mean(cube):
    add_month_number(cube, 'time')
    add_year(cube, 'time')
    cube = cube.aggregated_by(('month_number', 'year'), MEAN)
    cube.remove_coord('month_number')
    cube.remove_coord('year')
    return cube


def _create_areacello(cfg, sample_cube, glob_attrs, out_dir):
    if not cfg['custom'].get('create_areacello', False):
        return
    var_info = cfg['cmor_table'].get_variable('fx', 'areacello')
    lat_coord = sample_cube.coord('latitude')
    glob_attrs['mip'] = 'fx'
    cube = Cube(
        np.full(lat_coord.shape, cfg['custom']['grid_cell_size'], np.float32),
        standard_name=var_info.standard_name,
        long_name=var_info.long_name,
        var_name=var_info.short_name,
        units='m2',
    )
    cube.add_aux_coord(lat_coord, (0, 1))
    cube.add_aux_coord(sample_cube.coord('longitude'), (0, 1))
    cube.add_dim_coord(sample_cube.coord('projection_y_coordinate'), 0)
    cube.add_dim_coord(sample_cube.coord('projection_x_coordinate'), 1)
    cube.coord('projection_x_coordinate').var_name = 'x'
    cube.coord('projection_y_coordinate').var_name = 'y'
    fix_var_metadata(cube, var_info)
    set_global_atts(cube, glob_attrs)
    save_variable(
        cube, var_info.short_name, out_dir, cfg['attributes'], zlib=True
    )
