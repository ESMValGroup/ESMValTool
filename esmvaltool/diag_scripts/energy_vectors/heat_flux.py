import os
import sys
import logging

import numpy as np
import iris
import iris.coord_categorisation as ic

from esmvaltool.diag_scripts.shared import Variables, Datasets, run_diagnostic
from .common import low_pass_weights, lanczos_filter, IDENTIFY_DATASET

logger = logging.getLogger(os.path.basename(__file__))


class HeatFlux(object):

    def __init__(self, config):
        self.cfg = config
        self.datasets = Datasets(self.cfg)
        self.variables = Variables(self.cfg)
        self.window = self.cfg['window']

    def compute(self):

        va_paths = self.datasets.get_path_list(standard_name='northward_wind')
        for va_path in va_paths:
            va_info = self.datasets.get_dataset_info(va_path)
            logger.info(va_info)
            ta_info = {
                key: va_info[key]
                for key in IDENTIFY_DATASET if key in va_info
            }
            ta_info['standard_name'] =
            ta_path = self.datasets.get_path(**va_info)

            logger.info("Processing %s", ua_path)
            va_cube = iris.load_cube(va_path)
            ta_cube = iris.load_cube(ta_path)

            filter_weights = low_pass_weights(self.window, 1.0 / self.window)
            assert abs(sum(filter_weights)-1) < 1e-8

            ta_cube = ta_cube / H_low
            va_cube = va_cube / H_low

            logger.info("Calculating eddy heat flux")
            heat_flux_cube = self.eddy_heat_flux(ta_cube, va_cube, filter_weights)

            ic.add_month_number(VpTp, 'time', 'month_number')
            heat_flux_cube = heat_flux_cube.aggregated_by('month_number', iris.analysis.MEAN)
            iris.save(heat_flux_cube, outfile)

    def eddy_heat_flux(self, va_cube, ta_cube, filter_weights):
        '''
        calculate eddy_heat_flux from time series cubes of T and V.
        window and filter_weights required for the Lanczos filter

        Arguments:
        * T_cube:
            iris cube of air temperature data

        * V_cube:
            iris cube of V wind data

        * filter_weights:
            weights for lanczos filtering

        '''
        ta_cube_filterd = lanczos_filter(ta_cube, filter_weights)
        va_cube_filtered = lanczos_filter(va_cube, filter_weights)

        ta_high = ta_cube[window//2:-(window//2)] - ta_cube_filterd
        va_high = va_cube[window//2:-(window//2)] - va_cube_filtered

        heat_flux = lanczos_filter(ta_high * va_high, filter_weights)

        heat_flux.long_name = "Eddy Heat Flux (V'T')"
        heat_flux.attributes['filtering'] = \
            "Lanczos filtering with window width=%s (%s days)" % (
                len(filter_weights), self.window
            )
        heat_flux._var_name = 'VpTp'
        return heat_flux

    def _save(self, heat_flux_cube, dataset):
         if not self.cfg[NAMES.WRITE_NETCDF]:
            return
        logger.info("Saving results")
        subdir = os.path.join(
            self.cfg[NAMES.WORK_DIR],
            self.datasets.get_info(NAMES.PROJECT, dataset),
            self.datasets.get_info(NAMES.DATASET, dataset),
        )
        os.makedirs(subdir, exist_ok=True)
        iris.save(heat_flux_cube, os.path.join(subdir, 'heat_flux.nc'))


def main():
    with run_diagnostic() as config:
        HeatFlux(config).compute()


if __name__ == "__main__":
    main()
