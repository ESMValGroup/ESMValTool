"""Ocean heat content diagnostic"""

import calendar
import logging
import os

import iris
import iris.analysis
import iris.coord_categorisation
import iris.cube
import iris.quickplot as qplt
import iris.util
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from iris.analysis import SUM
from iris.coords import AuxCoord
from iris.cube import CubeList
from matplotlib import colors

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n

matplotlib.use('Agg')  # noqa

logger = logging.getLogger(os.path.basename(__file__))


class OceanHeatContent(object):
    """
    Ocean heat content diagnostic

    Allowed parameters:
    - min_depth: float=0
    - max_depth: float=np.inf

    Parameters
    ----------
    settings_file: str
        Path to the settings file

    """
    HEAT_CAPACITY = 4000
    WATER_DENSITY = 1020

    def __init__(self, config):
        self.cfg = config
        self.datasets = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.variables = esmvaltool.diag_scripts.shared.Variables(self.cfg)

        self.min_depth = self.cfg.get('min_depth', 0.)
        self.max_depth = self.cfg.get('max_depth', np.inf)
        self.compute_2d = self.cfg.get('compute_2D', True)
        self.compute_2d_monthly_clim = self.cfg.get('compute_2D_monthly_clim',
                                                    True)
        self.compute_1d = self.cfg.get('compute_1D', True)
        self.compute_1d_monthly_clim = self.cfg.get('compute_1D_monthly_clim',
                                                    True)

        self.min_color_scale = self.cfg.get('min_color_scale', None)
        self.max_color_scale = self.cfg.get('max_color_scale', None)
        self.log_color_scale = self.cfg.get('log_color_scale', False)

        self.color_map = self.cfg.get('color_map', 'plasma')

        self.color_intervals = self.cfg.get('color_intervals', 100)
        self.color_low = self.cfg.get('color_low', (1, 1, 1))
        self.color_high = self.cfg.get('color_high', (0.5, 0.05, 0.05))
        self.color_bad = self.cfg.get('color_bad', (0.7, 0.7, 0.7))

        self.plot_dpi = self.cfg.get('dpi', 200)

        if not self.max_color_scale:
            if self.max_depth == np.inf:
                depth = 6000 - self.min_depth
            else:
                depth = self.max_depth - self.min_depth
            self.max_color_scale = (
                depth * OceanHeatContent.WATER_DENSITY *
                OceanHeatContent.HEAT_CAPACITY * 300.
                )
        if self.min_color_scale is None:
            if self.max_depth == np.inf:
                depth = 6000 - self.min_depth
            else:
                depth = self.max_depth - self.min_depth
            self.min_color_scale = (
                depth * OceanHeatContent.WATER_DENSITY *
                OceanHeatContent.HEAT_CAPACITY * 260.
            )

        if not self.color_intervals:
            self.color_intervals = 1000
        if self.log_color_scale:
            self.norm = colors.LogNorm(vmin=self.min_color_scale,
                                       vmax=self.max_color_scale)
        else:
            self.norm = None
        if self.color_map == 'custom':
            self.cmap = colors.LinearSegmentedColormap.from_list(
                'ohc', (self.color_low, self.color_high),
                N=self.color_intervals
            )
            self.cmap.set_bad(color=self.color_bad)
        else:
            self.cmap = self.color_map

    def compute(self):
        """Compute diagnostic"""
        logger.info('Computing ocean heat content')
        ohc1d_compare_timeseries = CubeList()
        ohc1d_compare_monthly = CubeList()
        for filename in self.datasets:
            logger.info('Working in %s', filename)
            dataset_info = self.datasets.get_dataset_info(filename)
            thetao = iris.load_cube(filename,
                                    'sea_water_potential_temperature')
            if dataset_info[n.FX_FILES]['areacello']:
                constraint = iris.Constraint(
                    cube_func=lambda c: c.var_name == 'areacello'
                )
                area_cello = iris.load_cube(
                    dataset_info[n.FX_FILES]['areacello'], constraint
                )
            else:
                e1t = iris.load_cube(dataset_info['mesh_file'], 'e1t')
                e2t = iris.load_cube(dataset_info['mesh_file'], 'e2t')
                area_cello = iris.util.squeeze(e1t * e2t)

            self._compute_depth_weights(thetao)

            has_weight = iris.Constraint(depth_weight=lambda x: x > 0)
            thetao = thetao.extract(has_weight)
            depth_weight = thetao.coord('depth_weight').points

            ohc2d = CubeList()
            ohc1d = CubeList()
            final_weight = None
            logger.debug('Starting computation...')
            for time_slice in thetao.slices_over('time'):
                if final_weight is None:
                    index = time_slice.coord_dims('depth')[0]
                    final_weight = \
                        iris.util.broadcast_to_shape(depth_weight,
                                                     time_slice.shape,
                                                     (index,))
                time_slice.convert_units('K')
                ohc = time_slice.collapsed('depth', SUM,
                                           weights=final_weight)

                ohc.remove_coord('year')
                ohc.remove_coord('month_number')
                ohc.units = 'J m-2'
                ohc.var_name = 'ohc'
                ohc.standard_name = None
                ohc.long_name = 'Ocean Heat Content per area unit'
                if self.compute_2d:
                    self._plot(ohc, filename)

                ohc2d.append(ohc)
                if self.compute_1d:
                    ohc = ohc.collapsed(['latitude', 'longitude'], SUM,
                                        weights=area_cello.data)

                    ohc.long_name = 'Ocean Heat Content 1D'
                    ohc1d.append(ohc)

            del thetao

            logger.debug('Merging results...')
            ohc2d = ohc2d.merge_cube()
            if self.compute_2d:
                self._save_netcdf(ohc2d, filename)
            iris.coord_categorisation.add_month_number(ohc2d, 'time')
            iris.coord_categorisation.add_year(ohc2d, 'time')
            self._monthly_2d_clim(filename, ohc2d)

            if self.compute_1d:
                ohc1d = ohc1d.merge_cube()
                iris.coord_categorisation.add_month_number(ohc1d, 'time')
                iris.coord_categorisation.add_year(ohc1d, 'time')
                self._timeseries_1d(filename, ohc1d, ohc1d_compare_timeseries)
                self._monthly_1d_clim(filename, ohc1d, ohc1d_compare_monthly)
        del ohc1d
        del ohc2d
        self._comparison_plot(ohc1d_compare_timeseries, 'timeseries')
        self._comparison_plot(ohc1d_compare_monthly, 'monthly_clim')

    def _monthly_2d_clim(self, filename, ohc2d):
        if not self.compute_2d_monthly_clim:
            return

        ohc_clim = ohc2d.aggregated_by(('month_number'), iris.analysis.MEAN)
        if self.cfg[n.WRITE_NETCDF]:
            new_filename = os.path.basename(filename).replace(
                'thetao', 'ohc2D_monclim'
            )
            netcdf_path = os.path.join(self.cfg[n.WORK_DIR], new_filename)
            iris.save(ohc_clim, netcdf_path, zlib=True)
        if self.cfg[n.WRITE_PLOTS]:
            for month_slice in ohc_clim.slices_over('month_number'):
                month = month_slice.coord('month_number').points[0]
                self._get_clim_plot_filename(filename, month, 2)
                qplt.pcolormesh(
                    month_slice,
                    vmin=self.min_color_scale, vmax=self.max_color_scale,
                    cmap=self.cmap,
                    norm=self.norm,
                )
                plot_path = self._get_clim_plot_filename(
                    filename,
                    month,
                    2
                )
                plt.savefig(plot_path, dpi=self.plot_dpi)
                plt.close()

    def _monthly_1d_clim(self, filename, ohc1d, ohc1d_compare):
        if not self.compute_1d_monthly_clim:
            return

        ohc_clim = ohc1d.aggregated_by(('month_number'), iris.analysis.MEAN)
        if self.cfg[n.WRITE_NETCDF]:
            new_filename = os.path.basename(filename).replace(
                'thetao', 'ohc1D_monclim'
            )
            netcdf_path = os.path.join(self.cfg[n.WORK_DIR], new_filename)
            iris.save(ohc_clim, netcdf_path, zlib=True)
        if self.cfg[n.WRITE_PLOTS]:
            plt.plot(
                ohc_clim.data
            )
            plt.xlabel('Month')
            plt.ylabel('OHC (%s)' % ohc1d.units)
            plt.title(ohc_clim.long_name)
            plt.grid()
            plot_path = self._get_1D_plot_filename(filename, 'monthly_clim')
            plt.savefig(plot_path, dpi=self.plot_dpi)
            plt.close()
        ohc1d_compare.append(ohc_clim)

    def _timeseries_1d(self, filename, ohc1d, ohc1d_compare):
        if self.cfg[n.WRITE_PLOTS]:
            plt.plot(
                ohc1d.data
            )
            plt.xlabel('Month')
            plt.ylabel('OHC (%s)' % ohc1d.units)
            plt.title(ohc1d.long_name)
            plt.grid()
            plot_path = self._get_1D_plot_filename(filename, 'timeseries')
            plt.savefig(plot_path, dpi=self.plot_dpi)
            plt.close()
        ohc1d_compare.append(ohc1d)

    def _comparison_plot(self, ohc_compare, type_of_plot):
        logger.info('Creating comparison plot for %s', type_of_plot)
        for dataset in range(len(ohc_compare)):
            ohc_mean = ohc_compare[dataset].collapsed(
                'time', iris.analysis.MEAN
            )
            ohc_standard = ohc_compare[dataset] - ohc_mean
            plt.figure(1)
            plt.plot(
                ohc_compare[dataset].data,
                label='_'.join([
                    self.datasets.get_info_list(n.DATASET)[dataset],
                    str(self.datasets.get_info_list(n.START_YEAR)[dataset]),
                    str(self.datasets.get_info_list(n.END_YEAR)[dataset])
                ])
            )
            plt.figure(2)
            plt.plot(
                ohc_standard.data,
                label='_'.join([
                    self.datasets.get_info_list(n.DATASET)[dataset],
                    str(self.datasets.get_info_list(n.START_YEAR)[dataset]),
                    str(self.datasets.get_info_list(n.END_YEAR)[dataset])
                ])
            )
        script = self.cfg[n.SCRIPT]
        plt.figure(1)
        plt.title(script.replace('_', ' '))
        plt.xlabel('Month')
        plt.ylabel('OHC (%s)' % ohc_compare[0].units)
        plt.legend()
        plt.grid()
        plt.tight_layout()
        datasets = '_'.join(self.datasets.get_info_list(n.DATASET))
        out_type = self.cfg[n.OUTPUT_FILE_TYPE]
        plt_name = \
            '{script}_{type_of_plot}_comparison_{datasets}.{out_type}'.format(
                script=script,
                type_of_plot=type_of_plot,
                datasets=datasets,
                out_type=out_type
            )
        plt.figure(1).savefig(os.path.join(self.cfg[n.PLOT_DIR], plt_name))
        plt.close(1)
        plt.figure(2)
        plt.title(script.replace('_', ' '))
        plt.xlabel('Month')
        plt.ylabel('OHC (%s)' % ohc_compare[0].units)
        plt.legend()
        plt.grid()
        plt.tight_layout()
        stan_plt_name = \
            '{script}_{type_of_plot}_comparison_standarized_{datasets}' \
            '.{out_type}'.format(
                script=script,
                type_of_plot=type_of_plot,
                datasets=datasets,
                out_type=out_type
            )
        plt.figure(2).savefig(
            os.path.join(self.cfg[n.PLOT_DIR], stan_plt_name)
        )
        plt.close(2)

    def _get_clim_plot_filename(self, filename, month, dimensions):
        dataset = self.datasets.get_info(n.DATASET, filename)
        project = self.datasets.get_info(n.PROJECT, filename)
        ensemble = self.datasets.get_info(n.ENSEMBLE, filename)
        start = self.datasets.get_info(n.START_YEAR, filename)
        end = self.datasets.get_info(n.END_YEAR, filename)
        out_type = self.cfg[n.OUTPUT_FILE_TYPE]
        month = calendar.month_abbr[month]
        plot_filename = 'ohc{dimensions}D_{month}_{project}_{dataset}_' \
                        '{ensemble}_{start}-{end}' \
                        '.{out_type}'.format(dataset=dataset,
                                             project=project,
                                             ensemble=ensemble,
                                             start=start,
                                             end=end,
                                             out_type=out_type,
                                             month=month,
                                             dimensions=dimensions)
        plot_path = os.path.join(
            self._get_plot_folder(filename, 'monthly_clim'),
            plot_filename
        )
        return plot_path

    def _get_1D_plot_filename(self, filename, type_of_plot):
        dataset = self.datasets.get_info(n.DATASET, filename)
        project = self.datasets.get_info(n.PROJECT, filename)
        ensemble = self.datasets.get_info(n.ENSEMBLE, filename)
        start = self.datasets.get_info(n.START_YEAR, filename)
        end = self.datasets.get_info(n.END_YEAR, filename)
        out_type = self.cfg[n.OUTPUT_FILE_TYPE]
        plot_filename = 'ohc1D_{project}_{dataset}_' \
                        '{ensemble}_{start}-{end}' \
                        '.{out_type}'.format(dataset=dataset,
                                             project=project,
                                             ensemble=ensemble,
                                             start=start,
                                             end=end,
                                             out_type=out_type)
        plot_path = os.path.join(
            self._get_plot_folder(
                filename,
                '{type_of_plot}'.format(type_of_plot=type_of_plot)
            ),
            plot_filename
        )
        return plot_path

    def _compute_depth_weights(self, thetao):
        depth = thetao.coord('depth')
        if not depth.has_bounds():
            depth.guess_bounds()
        depth_weight = np.zeros(depth.shape)
        for current_depth in range(depth_weight.size):
            high = depth.bounds[current_depth, 0]
            low = depth.bounds[current_depth, 1]
            if low <= self.min_depth:
                continue
            if high >= self.max_depth:
                continue
            if low > self.max_depth:
                low = self.max_depth
            if high < self.min_depth:
                high = self.min_depth
            size = low - high
            if size < 0:
                size = 0
            depth_weight[current_depth] = (
                size *
                OceanHeatContent.HEAT_CAPACITY *
                OceanHeatContent.WATER_DENSITY
            )
        thetao.add_aux_coord(AuxCoord(var_name='depth_weight',
                                      points=depth_weight),
                             thetao.coord_dims(depth))

    def _save_netcdf(self, ohc2d, filename):
        if self.cfg[n.WRITE_NETCDF]:
            new_filename = os.path.basename(filename).replace('thetao',
                                                              'ohc')
            netcdf_path = os.path.join(self.cfg[n.WORK_DIR],
                                       new_filename)
            iris.save(ohc2d, netcdf_path, zlib=True)

    def _plot(self, ohc2d, filename):
        iris.FUTURE.cell_datetime_objects = True
        if self.cfg[n.WRITE_PLOTS]:
            qplt.pcolormesh(
                ohc2d,
                vmin=self.min_color_scale, vmax=self.max_color_scale,
                cmap=self.cmap,
                norm=self.norm,
            )
            datetime = ohc2d.coord('time').cell(0).point
            plot_path = self._get_plot_name(filename, datetime)
            plt.savefig(plot_path, dpi=self.plot_dpi)
            plt.close()

    def _get_plot_name(self, filename, datetime):
        dataset = self.datasets.get_info(n.DATASET, filename)
        project = self.datasets.get_info(n.PROJECT, filename)
        exp = self.datasets.get_info(n.EXP, filename)
        ensemble = self.datasets.get_info(n.ENSEMBLE, filename)
        time_str = datetime.strftime('%Y%m')
        out_type = self.cfg[n.OUTPUT_FILE_TYPE]
        plot_filename = 'ohc2D_{project}_{exp}_{dataset}_' \
                        '{ensemble}_{time_str}' \
                        '.{out_type}'.format(dataset=dataset,
                                             project=project,
                                             exp=exp,
                                             ensemble=ensemble,
                                             time_str=time_str,
                                             out_type=out_type)

        plot_path = os.path.join(
            self._get_plot_folder(filename, 'timestep'),
            plot_filename)
        return plot_path

    def _get_plot_folder(self, filename, plot_type):
        dataset = self.datasets.get_info(n.DATASET, filename)
        project = self.datasets.get_info(n.PROJECT, filename)
        exp = self.datasets.get_info(n.EXP, filename)
        ensemble = self.datasets.get_info(n.ENSEMBLE, filename)

        plot_path = os.path.join(
            self.cfg[n.PLOT_DIR], plot_type, project, exp, dataset, ensemble
        )
        if not os.path.isdir(plot_path):
            os.makedirs(plot_path)
        return plot_path


if __name__ == '__main__':
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        OceanHeatContent(config).compute()
