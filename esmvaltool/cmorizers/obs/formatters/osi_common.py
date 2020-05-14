"""Common functionalities for OSI-450 dataset cmorization."""

import logging
import os
import glob
from datetime import datetime, timedelta
from calendar import monthrange, isleap

import numpy as np
import iris
import iris.exceptions
from iris.cube import Cube, CubeList
from iris.coords import AuxCoord
from iris.coord_categorisation import add_day_of_year
from esmvalcore.preprocessor import monthly_statistics

from esmvaltool.cmorizers.obs.utilities import (
    set_global_atts, convert_timeunits, fix_var_metadata, save_variable)

logger = logging.getLogger(__name__)


class OSICmorizer():
    """Cmorizer for OSI-450 datasets."""

    def __init__(self, in_dir, out_dir, cfg, hemisphere):
        self.in_dir = in_dir
        self.out_dir = out_dir
        self.cfg = cfg
        self.hemisphere = hemisphere
        self.min_days = self.cfg['custom'].get('min_days', 50)

    def cmorize(self):
        """Cmorize OSI-450 or OSI-409 dataset."""
        logger.info(
            "Starting cmorization for Tier%s OBS files: %s",
            self.cfg['attributes']['tier'],
            self.cfg['attributes']['dataset_id'])
        logger.info("Input data from: %s", self.in_dir)
        logger.info("Output will be written to: %s", self.out_dir)

        # run the cmorization
        first_run = True
        for var, vals in self.cfg['variables'].items():
            var_info = {}
            for mip in vals['mip']:
                var_info[mip] = self.cfg['cmor_table'].get_variable(mip, var)
            file_pattern = '{0}_{1}_{2}_*.nc'.format(
                vals['raw'], self.hemisphere, vals['grid']
            )
            for year in os.listdir(self.in_dir):
                year = int(year)
                logger.info(
                    "CMORizing var %s for year %s", var, year
                )
                raw_info = {
                    'name': vals['raw'],
                    'file': os.path.join(
                        self.in_dir, str(year), '??', file_pattern)
                }
                self._extract_variable(var_info, raw_info, year, vals['mip'])
                if first_run:
                    sample_file = glob.glob(os.path.join(
                        self.in_dir, str(year), '01', file_pattern))[0]
                    cube = iris.load_cube(
                        sample_file,
                        iris.Constraint(
                            # pylint: disable=cell-var-from-loop
                            cube_func=lambda c: c.var_name == raw_info['name'])
                    )
                    self._create_areacello(cube)
                    first_run = False

    def _extract_variable(self, var_infos, raw_info, year, mips):
        """Extract to all vars."""
        cubes = iris.load(
            raw_info['file'],
            iris.Constraint(cube_func=lambda c: c.var_name == raw_info['name'])
        )
        tracking_ids = self._unify_attributes(cubes)
        cube = cubes.concatenate_cube()
        del cubes
        if tracking_ids:
            cube.attributes['tracking_ids'] = tracking_ids
        cube.coord('projection_x_coordinate').var_name = 'x'
        cube.coord('projection_y_coordinate').var_name = 'y'
        lon_coord = cube.coord('longitude')
        lon_coord.points[lon_coord.points < 0] += 360
        source_cube = cube
        attrs = self.cfg['attributes']
        for mip in mips:
            var_info = var_infos[mip]
            attrs['mip'] = mip
            if var_info.frequency == 'mon':
                cube = monthly_statistics(source_cube)
                cube = self._fill_months(cube)
            elif var_info.frequency == 'day':
                cube = self._fill_days(source_cube, year)
            if not cube:
                continue
            logger.debug(cube)
            fix_var_metadata(cube, var_info)
            convert_timeunits(cube, year)
            set_global_atts(cube, attrs)
            self._try_remove_coord(cube, 'year')
            self._try_remove_coord(cube, 'day_of_year')
            self._try_remove_coord(cube, 'month_number')
            self._try_remove_coord(cube, 'day_of_month')
            save_variable(cube, var_info.short_name, self.out_dir, attrs)
        return cube

    @staticmethod
    def _try_remove_coord(cube, coord):
        try:
            cube.remove_coord(coord)
        except iris.exceptions.CoordinateNotFoundError:
            pass

    @staticmethod
    def _fill_months(cube):
        if cube.coord('time').shape[0] == 12:
            return cube
        cubes = CubeList(cube.slices_over('time'))
        model_cube = cubes[0].copy()
        for month in range(1, 13):
            month_constraint = iris.Constraint(
                # pylint: disable=cell-var-from-loop
                time=lambda cell: cell.point.month == month
            )
            if cubes.extract(month_constraint):
                continue
            cubes.append(
                OSICmorizer._create_nan_cube(model_cube, month, month=True))
        cube = cubes.merge_cube()
        return cube

    def _fill_days(self, cube, year):
        if cube.coord('time').shape[0] < self.min_days:
            logger.warning(
                'Only %s days available. Skip generation of daily files',
                cube.coord('time').shape[0]
            )
            return None
        total_days = 366 if isleap(year) else 365
        if cube.coord('time').shape[0] < total_days:
            cubes = OSICmorizer._add_nan_timesteps(cube, total_days)
            cube = cubes.merge_cube()
            cube.remove_coord('day_of_year')
            del cubes
        return cube

    @staticmethod
    def _add_nan_timesteps(cube, total_days):
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
            nan_cube = OSICmorizer._create_nan_cube(
                model_cube, day_of_year, month=False
            )
            add_day_of_year(nan_cube, 'time')
            cubes.append(nan_cube)
        del model_cube
        return cubes

    @staticmethod
    def _create_nan_cube(model_cube, num, month):
        nan_cube = model_cube.copy(
            np.ma.masked_all(model_cube.shape, dtype=model_cube.dtype)
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

    @staticmethod
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

    def _create_areacello(self, sample_cube):
        if not self.cfg['custom'].get('create_areacello', False):
            return
        var_info = self.cfg['cmor_table'].get_variable('fx', 'areacello')
        lat_coord = sample_cube.coord('latitude')
        self.cfg['attributes']['mip'] = 'fx'
        cube = Cube(
            np.full(
                lat_coord.shape,
                self.cfg['custom']['grid_cell_size'],
                np.float32),
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
        set_global_atts(cube, self.cfg['attributes'])
        save_variable(
            cube, var_info.short_name, self.out_dir,
            self.cfg['attributes'], zlib=True
        )
