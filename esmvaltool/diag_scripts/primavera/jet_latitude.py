import os
import logging

import matplotlib.pyplot as plt

import dask.array as da
import numpy as np

import iris
import iris.cube
import iris.analysis
import iris.util
import iris.coord_categorisation
import iris.quickplot as qp

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import group_metadata

logger = logging.getLogger(os.path.basename(__file__))


class JetLatitude(object):
    def __init__(self, config):
        self.cfg = config
        self.filenames = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.output_name = None
        self.target_grid = self.cfg.get('target_grid')
        self.grid_cube = None
        self.sftlf = None

    def compute(self):
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        lanczos_weight = np.load(
            '/home/users/panos/work_space/jetlat/python/LF_weights.npy')
        for alias in data:
            ua = iris.load_cube(data[alias][0]['filename'])

            ua_filtered = np.apply_along_axis(
                lambda m: np.convolve(m, lanczos_weight, mode='same'),
                axis=ua.coord_dims('time')[0],
                arr=ua.core_data()
            )
            ua_filtered = ua.copy(ua_filtered)

            ua_max = ua_filtered.collapsed('latitude', iris.analysis.MAX)
            ua_max_lat = np.argmax(
                ua_filtered.core_data(),
                axis=ua.coord_dims('latitude')[0]
            )
            ua_max_lat = ua_max.copy(ua_max_lat)

            wind_smooth = self._compute_histogram(ua_max)
            wind_anom = ua_max - wind_smooth

            qp.plot(wind_smooth)
            plt.savefig(os.path.join(
                self.cfg[n.PLOT_DIR], '{}_wind_climatolgy.png'.format(alias)))
            plt.close()
            qp.plot(wind_anom)
            plt.savefig(os.path.join(
                self.cfg[n.PLOT_DIR], '{}_wind_anomaliess.png'.format(alias)))
            plt.close()

            lat_smooth = self._compute_histogram(ua_max_lat)
            qp.plot(lat_smooth)
            plt.savefig(os.path.join(
                self.cfg[n.PLOT_DIR], '{}_lat.png'.format(alias)))
            plt.close()

    def _compute_histogram(self, data):
        iris.coord_categorisation.add_day_of_year(data, 'time')
        iris.coord_categorisation.add_year(data, 'time')
        data.collapsed('year', iris.analysis.MEAN)

        transform = np.fft.rfft(data.data)
        transform[3:np.size(transform)] = 0
        smooth = np.fft.irfft(transform)
        smooth_data = data.copy(smooth)
        logger.debug(smooth_data)
        return smooth_data


def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        JetLatitude(config).compute()


if __name__ == '__main__':
    main()
