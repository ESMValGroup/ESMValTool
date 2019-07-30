import os
import logging
import string

import iris
import iris.cube
import iris.analysis
import iris.util

import numpy as np

from dask import array as da

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import group_metadata


logger = logging.getLogger(os.path.basename(__file__))

class LandAtmosInteractions(object):
    def __init__(self, config):
        self.cfg = config
        self.filenames = esmvaltool.diag_scripts.shared.Datasets(self.cfg)

    def compute(self):
        data = group_metadata(self.cfg['input_data'].values(), 'dataset')
        for dataset in data:
            # reorder vars in the recipe?
            clt = iris.load_cube(data[dataset][0]['filename'])

            tas = iris.load_cube(data[dataset][1]['filename'])

            hfls = iris.load_cube(data[dataset][2]['filename'])

            hfss = iris.load_cube(data[dataset][3]['filename'])

            rsds = iris.load_cube(data[dataset][4]['filename'])

            rlds = iris.load_cube(data[dataset][5]['filename'])

            hftotal = hfls + hfss

            evaporative_frac = hfls / hftotal

            clim_ef = evaporative_frac.collapsed('time', iris.analysis.MEAN)
            clim_ef.long_name = 'climatological_evaporative_fraction'

            stdv_ef = evaporative_frac.collapsed('time', iris.analysis.STD_DEV)
            stdv_ef.long_name = 'interannual_standard_deviation_of_evaporative_fraction'

            corr_hfls_rad = iris.analysis.stats.pearsonr(hfls, rsds + rlds)
            corr_hfls_rad.long_name = 'interannual_correlation_between_latent_heat_flux_and_total_radiation_at_surface'

            corr_clt_tas = iris.analysis.stats.pearsonr(clt, tas)
            corr_clt_tas.long_name = 'interannual_correlation_between_temperature_and_cloud_cover'

            # add soil moisture to vars




    def main():
        with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
            LandAtmosInteractions(config).compute()