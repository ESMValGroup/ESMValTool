import os
import sys
import logging

import numpy as np
import iris
import iris.coord_categorisation as ic

from esmvaltool.diag_scripts.shared import run_diagnostic, group_metadata
from esmvaltool.diag_scripts.shared import names as n
from esmvaltool.diag_scripts.energy_vectors.common import low_pass_weights, lanczos_filter

logger = logging.getLogger(os.path.basename(__file__))


class HeatFlux(object):

    def __init__(self, config):
        self.cfg = config
        self.window = self.cfg['window']

    def compute(self):
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        for alias in data:
            logger.info("Processing %s", alias)
            var = group_metadata(data[alias], 'standard_name')
            va_cube = iris.load_cube(var['northward_wind'][0]['filename'])
            ta_cube = iris.load_cube(var['air_temperature'][0]['filename'])
            self._remove_extra_coords(va_cube)
            self._remove_extra_coords(ta_cube)

            filter_weights = low_pass_weights(
                self.window,
                freq=var['northward_wind'][0]['frequency']
            )
            assert abs(sum(filter_weights)-1) < 1e-8

            logger.info("Calculating eddy heat flux")
            heat_flux = self.eddy_heat_flux(ta_cube, va_cube, filter_weights)

            ic.add_month_number(heat_flux, 'time', 'month_number')
            heat_flux = heat_flux.aggregated_by('month_number', iris.analysis.MEAN)
            logger.info("Saving results")
            iris.save(
                heat_flux,
                os.path.join(
                    self.cfg[n.WORK_DIR],
                    'heatflux_{}.nc'.format(alias)
                )
            )

    def eddy_heat_flux(self, va_cube, ta_cube, filter_weights):
        """
        calculate eddy_heat_flux from time series cubes of T and V.
        window and filter_weights required for the Lanczos filter

        Arguments:
        * T_cube:
            iris cube of air temperature data

        * V_cube:
            iris cube of V wind data

        * filter_weights:
            weights for lanczos filtering

        """
        window_size = len(filter_weights)
        va_cube.coord('time').bounds = None
        ta_cube.coord('time').bounds = None

        ta_cube_filtered = lanczos_filter(ta_cube, filter_weights)
        va_cube_filtered = lanczos_filter(va_cube, filter_weights)
        half_window = window_size // 2
        ta_high = ta_cube[half_window:-half_window] - ta_cube_filtered
        va_high = va_cube[half_window:-half_window] - va_cube_filtered

        heat_flux = lanczos_filter(ta_high * va_high.data, filter_weights)

        heat_flux.long_name = "Eddy Heat Flux (V'T')"
        heat_flux.attributes['filtering'] = \
            "Lanczos filtering with window width=%s (%s days)" % (
                len(filter_weights), self.window
            )
        heat_flux._var_name = 'VpTp'
        return heat_flux

    @staticmethod
    def _remove_extra_coords(cube):
        cube.remove_coord('year')
        cube.remove_coord('day_of_year')
        cube.remove_coord('month_number')
        cube.remove_coord('day_of_month')

    def _save(self, alias, heat_flux_cube):
        if not self.cfg[n.WRITE_NETCDF]:
            return
        logger.info("Saving results")
        subdir = os.path.join(
            self.cfg[n.WORK_DIR],
            alias,
        )
        os.makedirs(subdir, exist_ok=True)
        iris.save(heat_flux_cube, os.path.join(subdir, 'heat_flux.nc'))


def main():
    with run_diagnostic() as config:
        HeatFlux(config).compute()


if __name__ == "__main__":
    main()
