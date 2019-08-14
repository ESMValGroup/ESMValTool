import os
import logging

import matplotlib.pyplot as plt

import numpy as np
from scipy import stats

import iris
import iris.cube
import iris.analysis
import iris.util
import iris.coord_categorisation

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
        lanczos_weight = np.load(os.path.expandvars(
            '$HOME/jetlat/python/LF_weights.npy'
        ))
        for alias in data:
            ua = iris.load_cube(data[alias][0]['filename'])
            iris.coord_categorisation.add_season(ua, 'time')

            ua_filtered = np.apply_along_axis(
                lambda m: np.convolve(m, lanczos_weight, mode='same'),
                axis=ua.coord_dims('time')[0],
                arr=ua.core_data()
            )
            ua = ua.copy(ua_filtered)

            wind = ua.collapsed('latitude', iris.analysis.MAX)
            wind.var_name = 'jet'
            wind.standard_name = None
            wind.long_name = 'Jet speed'

            latitude = np.argmax(
                ua.data,
                axis=ua.coord_dims('latitude')[0]
            )
            del ua
            del ua_filtered
            latitude = wind.copy(latitude)
            latitude.var_name = 'lat'
            latitude.standard_name = 'latitude'
            latitude.long_name = 'Jet latitude'
            latitude.units = 'degrees_north'

            lat_bounds = latitude.coord('latitude').bounds
            logger.debug(wind)
            logger.debug(latitude)
            interval = 2.5
            lat_bins = np.arange(
                (lat_bounds.min() // interval) * interval,
                ((lat_bounds.max() // interval) + 1) * interval,
                interval,
            )
            wind_bins = np.arange(0, 41, 1)
            self._compute_var(alias, wind, wind_bins, data[alias][0])
            self._compute_var(alias, latitude, lat_bins + 1.25, data[alias][0])

    def _compute_var(self, alias, data, bins, metadata):
        var_name = data.var_name
        clim = self._smooth_daily_clim(data)
        season_clim = clim.aggregated_by('season', iris.analysis.MEAN)

        anom = data.data
        day_year = data.coord('day_of_year').points
        for day_slice in clim.slices_over('day_of_year'):
            num_day = day_slice.coord('day_of_year').points[0]
            indexes = day_year == num_day
            anom[indexes] = anom[indexes] - day_slice.data
        anom = data.copy(anom)

        for season_slice in season_clim.slices_over('season'):
            season = season_slice.coord('season').points[0].split('|')[0]
            sea_data = anom.extract(iris.Constraint(season=season))
            sea_data = sea_data.copy(sea_data.data + season_slice.data)
            hist = self._compute_histogram(sea_data, bins)
            pdf = self._compute_pdf(sea_data, bins)

            sea_data.var_name = var_name
            if self.cfg[n.WRITE_NETCDF]:
                self._save_cube(
                    hist, '{}_{}hist_{}.nc'.format(alias, var_name, season))
                self._save_cube(
                    pdf, '{}_{}pdf_{}.nc'.format(alias, var_name, season))
            self._plot_histogram(alias, sea_data, hist, pdf, metadata)

    def _smooth_daily_clim(self, data):
        clim = data.aggregated_by('day_of_year', iris.analysis.MEAN)
        clim.remove_coord('time')
        clim.remove_coord('month_number')
        clim.remove_coord('day_of_month')
        clim.remove_coord('year')
        iris.util.promote_aux_coord_to_dim_coord(clim, 'day_of_year')
        clim_fft = np.fft.rfft(clim.data)
        clim_fft[3:np.size(clim_fft)] = 0
        clim_fft = np.fft.irfft(clim_fft)
        return clim.copy(clim_fft)

    def _compute_histogram(self, data, bins):
        hist, _ = np.histogram(data.data, bins=bins)
        cube = iris.cube.Cube(
            hist,
            var_name='n',
            long_name='ocurrences',
            units=1.0
        )
        cube.add_dim_coord(
            iris.coords.DimCoord(
                points=np.array([(bins[x-1] + bins[x]) / 2
                                 for x in range(1, len(bins))]),
                var_name=data.var_name,
                long_name=data.long_name,
                units=data.units,
                bounds=np.array([[bins[x-1], bins[x]]
                                 for x in range(1, len(bins))]),
            ),
            0
        )
        return cube

    def _compute_pdf(self, data, bins):
        kde = stats.gaussian_kde(data.data)
        samples = np.linspace(bins.min(), bins.max(), 100)
        kde.set_bandwidth(bw_method='silverman')
        kde.set_bandwidth(bw_method=kde.factor * 1.06)
        pdf = kde(samples)
        cube = iris.cube.Cube(
            pdf,
            var_name='p',
            long_name='probability',
            units=1.0
        )
        cube.add_dim_coord(
            iris.coords.DimCoord(
                points=samples,
                var_name=data.var_name,
                long_name=data.long_name,
                units=data.units,
            ),
            0
        )
        return cube

    def _plot_histogram(self, alias, anomalies, histogram, pdf, metadata):
        plot_config = self.cfg.get('plot', {})
        season = anomalies.coord('season').points[0]
        plt.figure(figsize=(14, 8), dpi=250)
        x_coord = histogram.coord(anomalies.long_name)
        spacing = x_coord.points[1] - x_coord.points[0]
        plt.bar(
            x_coord.points,
            histogram.data / (spacing * histogram.data.sum()),
            width=spacing, align='center',
            color=[0.4, 0.4, 0.4], alpha=1, edgecolor=[0.5, 0.5, 0.5]
        )
        plt.plot(
            (np.mean(anomalies.data), np.mean(anomalies.data)), (0, 1),
            color=[0.7, 0.7, 0.7], lw=3, ls='--'
        )
        plt.plot(pdf.coord(anomalies.long_name).points, pdf.data,
                 color=[0.8, 0.2, 0.2], lw=3, ls='-')
        plt.xlabel('{0.long_name} ({0.units})'.format(x_coord))
        plt.xlim(x_coord.bounds.min(), x_coord.bounds.max())
        plt.ylabel('Relative Frequency Density')
        y_max = plot_config.get('probability_limit', 0.14)
        plt.ylim(0, y_max)
        plt.yticks(
            np.arange(0, y_max, plot_config.get('probability_tick', 0.1))
        )
        plt.title('{} distribution for {} ({} hPa), {} ({}-{})'.format(
            anomalies.long_name, alias,
            int(anomalies.coord('air_pressure').points[0] / 100), season,
            metadata[n.START_YEAR], metadata[n.END_YEAR],
        ))
        plt.grid()
        plt.savefig(
            os.path.join(
                self.cfg[n.PLOT_DIR],
                '{}_{}_{}.png'.format(alias, anomalies.var_name, season)
            ),
            bbox_inches='tight'
        )

    def _save_cube(self, cube, filename):
        iris.save(
            cube,
            os.path.join(self.cfg[n.WORK_DIR], filename),
            zlib=True
        )


def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        JetLatitude(config).compute()


if __name__ == '__main__':
    main()
