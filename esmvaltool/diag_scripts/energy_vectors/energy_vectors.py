"""Energy vectors diagnostic."""
import os
import logging

import iris
from iris.analysis import MEAN
import iris.coord_categorisation as ic

from esmvaltool.diag_scripts.shared import run_diagnostic, group_metadata
from esmvaltool.diag_scripts.shared import names as NAMES
from esmvaltool.diag_scripts.shared.plot import quickplot
from esmvaltool.diag_scripts.energy_vectors.common import (
    low_pass_weights, lanczos_filter
)

logger = logging.getLogger(os.path.basename(__file__))


class EnergyVectors(object):
    """Energy evctor diagnostic class."""

    def __init__(self, config):
        self.cfg = config
        self.window = self.cfg['window']

    def compute(self):
        """Compute energy vectors."""
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        for alias in data:
            logger.info("Processing %s", alias)
            var = group_metadata(data[alias], 'standard_name')
            ua_cube = iris.load_cube(var['eastward_wind'][0]['filename'])
            va_cube = iris.load_cube(var['northward_wind'][0]['filename'])
            self._remove_extra_coords(ua_cube)
            self._remove_extra_coords(va_cube)

            filter_weights = low_pass_weights(
                self.window,
                freq=var['eastward_wind'][0]['frequency']
            )
            assert abs(sum(filter_weights)-1) < 1e-8

            logger.info("Calculating horizontal E-vectors")
            evector_x, evector_y = self._compute_evectors(
                ua_cube, va_cube, filter_weights
            )

            evector_x = self._get_monthly_mean(evector_x)
            evector_y = self._get_monthly_mean(evector_y)

            self._save(alias, evector_x, evector_y)
            self._plot(alias, evector_x, evector_y)

    def _compute_evectors(self, ua_cube, va_cube, filter_weights):
        '''
        Calculate horizontal components of E-vectors

        Calculate from time series cubes of U and V.
        Window and filter_weights required for the Lanczos filter

        Arguments:
        * ua_cube:
            iris cube of U wind data

        * va_cube:
            iris cube of V wind data

        * filter_weights:
            weights for lanczos filtering

        '''
        window_size = len(filter_weights)

        # filter U and V
        ua_cube_filtered = lanczos_filter(ua_cube, filter_weights)
        va_cube_filtered = lanczos_filter(va_cube, filter_weights)

        # get high frequency components by subtracting filtered values
        half_window = window_size//2

        ua_high = ua_cube[half_window:-(half_window)] - ua_cube_filtered
        va_high = va_cube[half_window:-(half_window)] - va_cube_filtered

        # calculate E-vector components, filtering them as we go
        evector_x = lanczos_filter(
            va_high ** 2 - ua_high ** 2, filter_weights
        )
        evector_y = lanczos_filter(-1.0 * ua_high * va_high, filter_weights)

        # Add some appropriate meta data
        units = ua_cube.units ** 2
        filtering = "Lanczos filtering with window width=%s (%s days)" % (
            window_size, self.window)

        evector_x.var_name = 'evector_x'
        evector_x.long_name = 'zonal component of E-Vector'
        evector_x.units = units
        evector_x.attributes['filtering'] = filtering

        evector_y.var_name = 'evector_y'
        evector_y.long_name = 'meridional component of E-Vector'
        evector_y.units = units
        evector_y.attributes['filtering'] = filtering
        return evector_x, evector_y

    @staticmethod
    def _remove_extra_coords(cube):
        cube.remove_coord('year')
        cube.remove_coord('day_of_year')
        cube.remove_coord('month_number')
        cube.remove_coord('day_of_month')

    @staticmethod
    def _get_monthly_mean(vector):
        ic.add_year(vector, 'time')
        ic.add_month_number(vector, 'time')
        return vector.aggregated_by(('month_number', 'year'), MEAN)

    def _save(self, alias, evector_x, evector_y):
        if not self.cfg[NAMES.WRITE_NETCDF]:
            return
        logger.info("Saving results")
        subdir = os.path.join(
            self.cfg[NAMES.WORK_DIR],
            alias,
        )
        os.makedirs(subdir, exist_ok=True)
        iris.save(evector_x, os.path.join(subdir, 'evector_x.nc'))
        iris.save(evector_y, os.path.join(subdir, 'evector_y.nc'))

    def _plot(self, alias, evector_x, evector_y):
        if not self.cfg[NAMES.WRITE_PLOTS]:
            return

        logger.info("Plotting results")
        evector_x = evector_x.collapsed('time', MEAN)
        evector_y = evector_y.collapsed('time', MEAN)
        subdir = os.path.join(
            self.cfg[NAMES.PLOT_DIR],
            alias,
        )
        os.makedirs(subdir, exist_ok=True)
        quickplot(
            evector_x,
            filename=os.path.join(
                subdir,
                'evector_x.{}'.format(self.cfg[NAMES.OUTPUT_FILE_TYPE])
            ),
            **(self.cfg.get(
                'quickplot',
                {'plot_type': 'pcolormesh', 'cmap': 'bwr'}
            ))
        )
        quickplot(
            evector_y,
            filename=os.path.join(
                subdir,
                'evector_y.{}'.format(self.cfg[NAMES.OUTPUT_FILE_TYPE])
            ),
            **(self.cfg.get(
                'quickplot',
                {'plot_type': 'pcolormesh', 'cmap': 'bwr'}
            ))
        )


def main():
    with run_diagnostic() as config:
        EnergyVectors(config).compute()


if __name__ == "__main__":
    main()
