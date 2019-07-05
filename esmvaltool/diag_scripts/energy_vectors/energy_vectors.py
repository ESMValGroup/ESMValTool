import sys
import os
import logging

import numpy as np
import iris
from iris.analysis import MEAN
import iris.coord_categorisation as ic

from esmvaltool.diag_scripts.shared import Variables, Datasets, run_diagnostic
from esmvaltool.diag_scripts.shared import names as NAMES
from .shared import low_pass_weights, lanczos_filter

logger = logging.getLogger(os.path.basename(__file__))


class EnergyVectors(object):

    def __init__(self, config):
        self.cfg = config
        self.datasets = Datasets(self.cfg)
        self.variables = Variables(self.cfg)
        self.window = self.cfg['window']

    def compute(self):

        for ua_path in self.datasets.get_data_list(standard_name='eastward_wind'):
            va_path = self.datasets.get_data(ua_path.replace('_ua', '_va'))

            logger.info("Processing %s", ua_path)

            ua_cube = iris.load_cube(ua_path)
            va_cube = iris.load_cube(va_path)

            filter_weights = low_pass_weights(self.window, freq='6hr')
            assert abs(sum(filter_weights)-1) < 1e-8

            logger.info("Calculating horizontal E-vectors")
            evector_x, evector_y = self._compute_evectors(
                ua_cube, va_cube, filter_weights
            )

            ic.add_month_number(evector_x, 'time', 'month_number')
            ic.add_month_number(evector_y, 'time', 'month_number')

            evector_x = evector_x.aggregated_by(('month_number', 'year'), MEAN)
            evector_y = evector_y.aggregated_by(('month_number', 'year'), MEAN)

            logger.info("Saving results")
            subdir = os.path.join(
                self.cfg[NAMES.WORK_DIR],
                self.datasets.get_info(NAMES.PROJECT, ua_path),
                self.datasets.get_info(NAMES.DATASET, ua_path),
            )
            os.makedirs(subdir, exist_ok=True)
            iris.save(evector_x, os.path.join(subdir, 'evector_x.nc'))
            iris.save(evector_y, os.path.join(subdir, 'evector_x.nc'))

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
        half_window = window_size//2-1
        ua_high = ua_cube[half_window:-(half_window)] - ua_cube_filtered
        va_high = va_cube[half_window:-(half_window)] - va_cube_filtered

        # calculate E-vector components, filtering them as we go
        evector_x = lanczos_filter(
            va_high ** 2 - ua_high ** 2,  filter_weights)
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


def main():
    with run_diagnostic() as config:
        EnergyVectors(config).compute()


if __name__ == "__main__":
    main()
