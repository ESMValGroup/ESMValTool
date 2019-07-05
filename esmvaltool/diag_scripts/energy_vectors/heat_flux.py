import os
import sys
import logging

import numpy as np
import iris
import iris.coord_categorisation as ic

from esmvaltool.diag_scripts.shared import Variables, Datasets, run_diagnostic
from .shared import low_pass_weights, strip_to_month_plus_window

logger = logging.getLogger(os.path.basename(__file__))


class HeatFlux(object):

    def __init__(self, config):
        self.cfg = config
        self.datasets = Datasets(self.cfg)
        self.variables = Variables(self.cfg)
        self.window = self.cfg['window']

    def compute(self):

        self.variables.get

        filter_weights = low_pass_weights(window, 1.0 / window)
        assert abs(sum(filter_weights)-1) < 1e-8

        T_low = strip_to_month_plus_window(T_low, window, month, year)
        V_low = strip_to_month_plus_window(V_low, window, month, year)
        H_low = strip_to_month_plus_window(H_low, window, month, year)

        T_low = T_low / H_low
        V_low = V_low / H_low

        logger.info("Calculating eddy heat flux")
        VpTp = self.eddy_heat_flux(T_low, V_low, window, filter_weights)

        ic.add_month_number(VpTp, 'time', 'month_number')
        VpTp = VpTp.aggregated_by('month_number', iris.analysis.MEAN)
        iris.save(VpTp, outfile)

    def eddy_heat_flux(self, T_cube, V_cube, window, filter_weights):
        '''
        calculate eddy_heat_flux from time series cubes of T and V.
        window and filter_weights required for the Lanczos filter

        Arguments:
        * T_cube:
            iris cube of air temperature data

        * V_cube:
            iris cube of V wind data

        * window:
            window size for filtering

        * filter_weights:
            weights for lanczos filtering

        '''
        T_cube_filtered = lanczos_filter(T_cube,filter_weights)
        V_cube_filtered = lanczos_filter(V_cube,filter_weights)

        Tp = T_cube[window//2-1:-(window//2-1)] - T_cube_filtered
        Vp = V_cube[window//2-1:-(window//2-1)] - V_cube_filtered

        VpTp_filtered = lanczos_filter(Vp * Tp, filter_weights)

        VpTp_filtered.long_name = "Eddy Heat Flux (V'T')"
        VpTp_filtered.attributes['filtering'] = "Lanczos filtering with window width=%s (%s days)"%(window, window/4.0)
        VpTp_filtered._var_name = 'VpTp'
        return VpTp_filtered


def main():
    with run_diagnostic() as config:
        HeatFlux(config).compute()


if __name__ == "__main__":
    main()