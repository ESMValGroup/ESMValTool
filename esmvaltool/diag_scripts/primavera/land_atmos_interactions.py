import os
import logging

import iris
import iris.cube
import iris.analysis
from iris.analysis import stats
import iris.util

from scipy import signal

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import group_metadata

from esmvalcore.preprocessor._regrid_esmpy import regrid
from esmvalcore.preprocessor._mask import mask_landsea
from esmvalcore.preprocessor._time import extract_season
logger = logging.getLogger(os.path.basename(__file__))


class LandAtmosInteractions(object):
    def __init__(self, config):
        self.cfg = config
        self.filenames = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.output_name = None
        self.target_grid = self.cfg.get('target_grid')
        self.grid_cube = None
        self.sftlf = None

    def compute(self):
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        for alias in data:
            var = group_metadata(data[alias], 'short_name')
            hfls = iris.load_cube(var['hfls'][0]['filename'])
            hfss = iris.load_cube(var['hfss'][0]['filename'])
            rsds = iris.load_cube(var['rsds'][0]['filename'])
            rlds = iris.load_cube(var['rlds'][0]['filename'])
            clt = iris.load_cube(var['clt'][0]['filename'])
            tas = iris.load_cube(var['tas'][0]['filename'])
            mrso = iris.load_cube(var['mrso'][0]['filename'])
            self.sftlf = [var['mrso'][0]['fx_files']['sftlf']]

            self.grid_cube = iris.load_cube(self.target_grid)

            clim_ef, stdv_ef = self.evaporative_fraction_stats(hfls, hfss)
            self.output_name = 'Evaporative_fraction_statistics'
            self.save([clim_ef, stdv_ef], alias, data)

            metrics = self.compute_correlation_metrics(
                hfls, rsds, rlds, clt, tas, mrso
                )
            self.output_name = 'Correlation_metrics'
            self.save(metrics, alias, data)

            # Compute the metrics but with detrended data.
            hfls = hfls.copy(signal.detrend(hfls.lazy_data(), axis=0))
            rsds = rsds.copy(signal.detrend(rsds.lazy_data(), axis=0))
            rlds = rlds.copy(signal.detrend(rlds.lazy_data(), axis=0))
            clt = clt.copy(signal.detrend(clt.lazy_data(), axis=0))
            tas = tas.copy(signal.detrend(tas.lazy_data(), axis=0))
            mrso = mrso.copy(signal.detrend(mrso.lazy_data(), axis=0))
            detrended_metrics = self.compute_correlation_metrics(
                hfls, rsds, rlds, clt, tas, mrso
                )
            self.output_name = 'Detrended_correlation_metrics'
            self.save(detrended_metrics, alias, data)

    def save(self, cubelist, alias, data):
        if self.output_name:
            dataset = data[alias][0]['dataset']
            experiment = data[alias][0]['exp']
            ensemble = data[alias][0]['ensemble']
            start_year = data[alias][0]['start_year']
            end_year = data[alias][0]['end_year']
            filename = '{output}_' \
                       '{dataset}_' \
                       '{experiment}_' \
                       '{ensemble}_' \
                       '{start_year}_' \
                       '{end_year}.nc'.format(output=self.output_name,
                                              dataset=dataset,
                                              experiment=experiment,
                                              ensemble=ensemble,
                                              start_year=start_year,
                                              end_year=end_year)

            iris.save(cubelist, os.path.join(self.cfg[n.WORK_DIR], filename))

    def compute_correlation_metrics(self, hfls, rsds, rlds, clt, tas, mrso):
        corr_hfls_rad = self.heat_radiation_stats(hfls, rsds, rlds)
        corr_clt_tas = self.cloud_temperature_stats(clt, tas)
        corr_hfls_mrso, corr_mrso_mam_jja = self.heat_soil_stats(hfls, mrso)
        return [corr_hfls_rad, corr_clt_tas,
                corr_hfls_mrso, corr_mrso_mam_jja]

    def evaporative_fraction_stats(self, hfls, hfss):
        hftotal = hfls + hfss
        evaporative_frac = hfls / hftotal
        clim = evaporative_frac.collapsed('time', iris.analysis.MEAN)
        clim = self.mask_and_regrid(clim)
        clim.long_name = 'climatological_evaporative_fraction'
        stdv = evaporative_frac.collapsed('time', iris.analysis.STD_DEV)
        stdv = self.mask_and_regrid(stdv)
        stdv.long_name = ('interannual_standard_deviation_'
                          'of_evaporative_fraction')
        return clim, stdv

    def heat_radiation_stats(self, hfls, rsds, rlds):
        rtotal = rsds + rlds
        correlation = stats.pearsonr(hfls, rtotal, 'time')
        correlation = self.mask_and_regrid(correlation)
        correlation.long_name = ('interannual_correlation_between'
                                 '_latent_heat_flux_and_total_radiation'
                                 '_at_surface')
        return correlation

    def cloud_temperature_stats(self, clt, tas):
        correlation = stats.pearsonr(clt, tas, 'time')
        correlation = self.mask_and_regrid(correlation)
        correlation.long_name = ('interannual_correlation_between'
                                 '_temperature_and_cloud_cover')
        return correlation

    def heat_soil_stats(self, hfls, mrso):
        mrso_mam = extract_season(mrso, 'MAM')
        mrso_jja = extract_season(mrso, 'JJA')
        dummy = mrso_jja.copy(mrso_mam.lazy_data())
        correlation = stats.pearsonr(hfls, mrso_jja, 'time')
        correlation = self.mask_and_regrid(correlation)
        correlation.long_name = ('interannual_correlation_between'
                                 '_latent_heat_flux_and_soil_moisture')
        mrso_correlation = stats.pearsonr(dummy, mrso_jja, 'time')
        mrso_correlation = self.mask_and_regrid(mrso_correlation)
        mrso_correlation.long_name = ('interannual_correlation'
                                      '_between_soil_moisture_'
                                      'in_spring_and_soil_moisture'
                                      '_in_summer')
        return correlation, mrso_correlation

    def mask_and_regrid(self, var):
        var = mask_landsea(var, self.sftlf, 'sea')
        var = regrid(var,
                     self.grid_cube,
                     method='area_weighted')
        return var

def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        LandAtmosInteractions(config).compute()


if __name__ == '__main__':
    main()
