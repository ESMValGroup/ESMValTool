"""Blocking diagnostic"""
import os
import logging
import itertools
import calendar

import numpy as np
from numba import vectorize

import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import colors
from matplotlib.lines import Line2D

import iris
import iris.time
import iris.util
import iris.coord_categorisation
import iris.analysis
import iris.analysis.stats
import iris.coords
import iris.quickplot
import cartopy.crs as ccrs

import skill_metrics

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvalcore.preprocessor import regrid

logger = logging.getLogger(os.path.basename(__file__))


class Blocking(object):
    """
    Blocking diagnostic

    Allowed parameters:
    - central_latitude: float=60
    - span: float=20
    - offset: float=5
    - north_threshold: float=-10
    - south_threshold: float=0
    - months: list=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    - max_color_scale: int=15
    - smoothing_window: int=3


    Parameters
    ----------
    settings_file: str
        Path to the settings file

    """

    def __init__(self, conf):
        self.cfg = conf
        self.datasets = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.variables = esmvaltool.diag_scripts.shared.Variables(self.cfg)
        for filename in self.datasets:
            self.reference_dataset = self._get_reference_dataset(
                self.datasets.get_info('reference_dataset', filename)
            )
            break

        self.compute_1d = self.cfg.get('compute_1d', True)
        self.compute_2d = self.cfg.get('compute_2d', True)

        self.span = self.cfg.get('span', 20.0)
        self.north_threshold = self.cfg.get('north_threshold', -10.0)
        self.south_threshold = self.cfg.get('south_threshold', 0.0)
        self.smoothing_window = self.cfg.get('smoothing_window', 0)
        self.persistence = self.cfg.get('persistence', 1)

        # 1D configuration
        self.central_latitude = self.cfg.get('central_latitude', 60.0)
        self.offset = self.cfg.get('offset', 5.0)

        self.max_color_scale = self.cfg.get('max_color_scale', 15)

        self.min_latitude = self.central_latitude - self.span - self.offset
        self.max_latitude = self.central_latitude + self.span + self.offset

        self.blocking_2D = {}
        self._latitude_data = {}

        @vectorize()
        def _get_index(high, central, low,
                       north_distance, south_distance,
                       north_threshold, south_threshold):
            if ((high - central) / north_distance) > north_threshold:
                return 0
            if ((central - low) / south_distance) < south_threshold:
                return 0
            return 1
        self._compute_index = _get_index

    def compute(self):
        """Compute blocking diagnostic"""
        logger.info('Computing blocking')

        logger.info('Reference dataset %s', self.reference_dataset)
        reference_1d, reference_2d = \
            self._get_blocking_indices(self.reference_dataset)

        skills = {'total': {}}
        for month in range(1, 13):
            skills[month] = {}

        datasets = [dataset for dataset in self.datasets
                    if dataset != self.reference_dataset]

        for filename in datasets:
            logger.info('Dataset %s', filename)
            dataset_1d, dataset_2d = self._get_blocking_indices(filename)
            if self.compute_2d:
                cmap = colors.LinearSegmentedColormap.from_list(
                    'mymap',
                    ((0.1, 0.1, 0.7), (1, 1, 1), (0.7, 0.1, 0.09)),
                    N=2 * self.max_color_scale
                )
                projection = ccrs.NorthPolarStereo()
                min_lat = np.min(reference_2d.coord('latitude').bounds)
                max_lat = np.max(reference_2d.coord('latitude').bounds)
                dataset_2d = regrid(dataset_2d, reference_2d, 'linear')
                diff = dataset_2d - reference_2d
                diff.long_name = 'Differences between model and ' \
                                 'reference in blocking index'
                self.plot_differences(
                    filename, diff, cmap, projection, min_lat, max_lat
                )
                taylor_data = skill_metrics.taylor_statistics(
                    np.ravel(dataset_2d.data),
                    np.ravel(reference_2d.data),
                )
                corr = iris.analysis.stats.pearsonr(
                    reference_2d, dataset_2d).data
                taylor_data['ccoef'][1] = corr
                skills['total'][filename] = taylor_data
                logger.debug(dataset_2d)
                for month_slice in dataset_2d.slices_over('month_number'):
                    month = month_slice.coord('month_number').points[0]
                    ref_slice = reference_2d.extract(
                        iris.Constraint(month_number=month)
                    )
                    corr = iris.analysis.stats.pearsonr(
                        ref_slice, month_slice).data
                    taylor_data = skill_metrics.taylor_statistics(
                        np.ravel(month_slice.data),
                        np.ravel(ref_slice.data),
                    )
                    taylor_data['ccoef'][1] = corr
                    skills[month][filename] = taylor_data

        if self.cfg[n.WRITE_PLOTS] and self.compute_2d:
            self.create_comparison_plot(
                datasets, skills['total'], 'blocking2D'
            )
            for month in range(1, 13):
                metrics = skills[month]
                self.create_comparison_plot(
                    datasets, metrics, 'blocking_2D_{}'.format(month)
                )

    def _get_blocking_indices(self, filename):
        result = self._blocking(filename)
        index_1d = self._blocking_1d(filename)
        index_2d = self._blocking_2d(filename, result)
        return (index_1d, index_2d)

    def _get_reference_dataset(self, reference_dataset):
        for filename in self.datasets:
            dataset = self.datasets.get_info(n.DATASET, filename)
            if dataset == reference_dataset:
                return filename
        raise ValueError(
            'Reference dataset "{}" not found'.format(reference_dataset)
        )

    def _blocking(self, filename):
        zg500 = iris.load_cube(filename, 'geopotential_height')
        for coord in zg500.coords():
            coord.points
            if coord.has_bounds():
                coord.bounds
            elif coord.standard_name in ('latitude', 'longitude'):
                coord.guess_bounds()
        iris.coord_categorisation.add_month(zg500, 'time')
        lat = zg500.coord('latitude')
        lat_max = np.max(lat.points)
        lat_min = np.min(lat.points)
        if self.compute_1d:
            lat_values = []
            for displacement in [-self.offset, 0, self.offset]:
                central = self.central_latitude + displacement
                lat_values.append(
                    lat.cell(lat.nearest_neighbour_index(central))
                )

        if self.compute_2d:
            latitudes = lat.points
        else:
            latitudes = (self.central_latitude - self.offset,
                         self.central_latitude,
                         self.central_latitude + self.offset)

        blocking = iris.cube.CubeList()
        self._latitude_data = {}
        total_years = len(set(zg500.coord('year').points))
        block1d_data = None
        for lat_point in latitudes:
            if lat_point + self.span > lat_max:
                continue
            if lat_point - self.span < lat_min:
                continue
            logger.debug('Computing blocking for lat %d', lat_point)
            blocking_index = self._compute_blocking(zg500, lat_point)
            if self.compute_1d and lat_point in lat_values:
                if block1d_data is not None:
                    block1d_data = np.logical_or(
                        block1d_data, blocking_index.data
                    )
                else:
                    block1d_data = blocking_index.data

            blocking_index = blocking_index.aggregated_by(
                'month_number', iris.analysis.SUM) / total_years
            blocking.append(blocking_index)
            blocking_index.long_name = 'Blocking index'
            blocking_index.units = 'Days per month'

        if self.compute_1d:
            blocking_cube = iris.cube.Cube(
                block1d_data.astype(int),
                var_name="blocking",
                attributes=None)

            blocking_cube.add_dim_coord(zg500.coord('time'), (0,))
            blocking_cube.add_dim_coord(zg500.coord('longitude'), (1,))
            blocking_cube.add_aux_coord(iris.coords.AuxCoord.from_coord(
                zg500.coord('latitude').copy([self.central_latitude])))
            iris.coord_categorisation.add_month_number(blocking_cube, 'time')
            result = blocking_cube.aggregated_by(
                'month_number', iris.analysis.SUM
            ) / total_years
            result.remove_coord('time')
            iris.util.promote_aux_coord_to_dim_coord(result, 'month_number')
            self.blocking_cube1d = result

        blocking_cube = blocking.merge_cube()
        return blocking_cube

    def _blocking_1d(self, filename):
        if not self.compute_1d:
            return None
        result = self._smooth_over_longitude(self.blocking_cube1d)
        result.units = 'days per month'
        result.var_name = 'blocking'
        result.long_name = 'Blocking 1D index'

        if self.cfg[n.WRITE_NETCDF]:
            new_filename = os.path.basename(filename).replace('zg',
                                                              'blocking1D')
            netcdf_path = os.path.join(self.cfg[n.WORK_DIR],
                                       new_filename)
            iris.save(result, target=netcdf_path, zlib=True)

        if self.cfg[n.WRITE_PLOTS]:
            iris.coord_categorisation.add_categorised_coord(
                result, 'month', result.coord('month_number'),
                lambda coord, x: calendar.month_abbr[x],
                units='no_unit')
            cmap = colors.LinearSegmentedColormap.from_list('mymap', (
                (1, 1, 1), (0.7, 0.1, 0.09)), N=self.max_color_scale)

            iris.quickplot.pcolormesh(result, coords=('longitude', 'month'),
                                      cmap=cmap, vmin=0,
                                      vmax=self.max_color_scale)
            plt.axis('tight')
            plt.yticks(range(result.coord('month').shape[0]),
                       result.coord('month').points)
            axes = plt.gca()
            axes.set_ylim((result.coord('month').shape[0] - 0.5, -0.5))

            plot_path = self._get_plot_name('blocking1D', filename)
            plt.savefig(plot_path)
            plt.close()
        # result.remove_coord('time')
        # iris.util.promote_aux_coord_to_dim_coord(result, 'month_number')
        return result

    def _get_plot_name(self, name, filename, month=None):
        alias = self.datasets.get_info(n.ALIAS, filename)
        start = self.datasets.get_info(n.START_YEAR, filename)
        end = self.datasets.get_info(n.END_YEAR, filename)

        plot_path = os.path.join(self.cfg[n.PLOT_DIR], alias)
        if not os.path.isdir(plot_path):
            os.makedirs(plot_path)

        if month is not None:
            name = '{}_{:02}'.format(name, month)
        out_type = self.cfg[n.OUTPUT_FILE_TYPE]

        plot_filename = f'{name}_{alias}_{start}-{end}.{out_type}'
        return os.path.join(plot_path, plot_filename)

    def _smooth_over_longitude(self, cube):
        if self.smoothing_window == 0:
            return cube
        logger.debug('Smoothing...')
        smoothed = iris.cube.CubeList()
        for lon_slice in cube.slices_over('longitude'):
            longitude = lon_slice.coord('longitude').points[0]
            lon_window = cube.intersection(
                longitude=(longitude - self.smoothing_window / 2.,
                           longitude + self.smoothing_window / 2.)
            )
            lon_mean = lon_window.collapsed('longitude', iris.analysis.MEAN)
            lon_slice.data[...] = lon_mean.data
            smoothed.append(lon_slice)
        cube = smoothed.merge_cube()
        logger.debug('Smoothing finished!')
        return cube

    def _compute_blocking(self, zg500, central_latitude):
        latitude = zg500.coord('latitude')

        def _get_lat_cell(coord, latitude):
            return coord.cell(coord.nearest_neighbour_index(latitude))

        central_lat = _get_lat_cell(latitude, central_latitude)
        low_lat = _get_lat_cell(latitude, central_latitude - self.span)
        high_lat = _get_lat_cell(latitude, central_latitude + self.span)

        zg_low = self._extract_lat(zg500, low_lat)
        zg_central = self._extract_lat(zg500, central_lat)
        zg_high = self._extract_lat(zg500, high_lat)

        north_distance = high_lat.point - central_lat.point
        south_distance = central_lat.point - low_lat.point

        blocking_index = self._compute_index(
            zg_high, zg_central, zg_low,
            north_distance, south_distance,
            self.north_threshold, self.south_threshold)

        del self._latitude_data[low_lat]

        blocking_cube = self._create_blocking_cube(
            blocking_index, zg500, central_latitude, central_lat.bound)

        if self.persistence > 1:
            self._apply_persistence(blocking_cube)
        return blocking_cube

    def _create_blocking_cube(self, blocking_index, zg500, central_latitude,
                              bounds):
        blocking_cube = iris.cube.Cube(
            blocking_index,
            var_name="blocking",
            units="Days per month",
            long_name="Blocking pattern",
            attributes=None, )
        blocking_cube.add_aux_coord(zg500.coord('month_number'), (0,))
        blocking_cube.add_aux_coord(iris.coords.AuxCoord.from_coord(
            zg500.coord('latitude').copy([central_latitude], bounds=bounds)))
        blocking_cube.add_dim_coord(zg500.coord('time'), (0,))
        blocking_cube.add_dim_coord(zg500.coord('longitude'), (1,))
        return blocking_cube

    def _extract_lat(self, zg500, latitude):
        if latitude not in self._latitude_data:
            lat_data = zg500.extract(iris.Constraint(latitude=latitude)).data
            self._latitude_data[latitude] = lat_data
        else:
            lat_data = self._latitude_data[latitude]
        return lat_data

    def _apply_persistence(self, blocking_cube):
        for point_slice in blocking_cube.slices('time'):
            grouped = ((k, sum(1 for _ in g))
                       for k, g in itertools.groupby(point_slice.data))
            index = 0
            for value, length in grouped:
                if value and length < self.persistence:
                    point_slice.data[index: index + length] = False
                index += length
        return blocking_cube

    def _blocking_2d(self, filename, blocking_index):
        if not self.compute_2d:
            return None
        if self.cfg[n.WRITE_NETCDF]:
            new_filename = os.path.basename(filename).replace('zg',
                                                              'blocking')
            netcdf_path = os.path.join(self.cfg[n.WORK_DIR],
                                       new_filename)
            iris.save(blocking_index, netcdf_path, zlib=True)

        if self.cfg[n.WRITE_PLOTS]:
            projection = ccrs.NorthPolarStereo()
            min_lat = np.min(blocking_index.coord('latitude').bounds)
            max_lat = np.max(blocking_index.coord('latitude').bounds)
            cmap = colors.LinearSegmentedColormap.from_list(
                'mymap',
                ((0.92, 0.92, 0.92), (0.7, 0.1, 0.09)),
                N=self.max_color_scale
            )
            for month_slice in blocking_index.slices_over('month_number'):
                month_number = month_slice.coord('month_number').points[0]
                month_name = calendar.month_name[month_number]
                logger.info('Plotting 2D blocking for ' + month_name)
                month_slice.long_name += ' (' + month_name.title() + ')'

                plt.figure()
                axes = plt.axes(projection=projection)
                axes.set_extent(
                    (-180, 180, min_lat, max_lat),
                    crs=ccrs.PlateCarree()
                )
                iris.quickplot.pcolormesh(
                    month_slice,
                    coords=('longitude', 'latitude'),
                    cmap=cmap, vmin=0, vmax=self.max_color_scale,
                )
                axes.coastlines()
                axes.gridlines(alpha=0.5, linestyle='--')
                theta = np.linspace(0, 2*np.pi, 100)
                center, radius = [0.5, 0.5], 0.5
                verts = np.vstack([np.sin(theta), np.cos(theta)]).T
                circle = mpath.Path(verts * radius + center)
                axes.set_boundary(circle, transform=axes.transAxes)

                plot_path = self._get_plot_name(
                    'blocking2D',
                    filename,
                    month_number
                )
                plt.savefig(plot_path, bbox_inches='tight', pad_inches=0.2,
                            dpi=500)
                plt.close()
        blocking_index.remove_coord('time')
        iris.util.promote_aux_coord_to_dim_coord(
            blocking_index, 'month_number'
        )
        blocking_index.coord('month_number').attributes.clear()
        return blocking_index

    def plot_differences(self, dataset, diff_cube, cmap, projection,
                         min_lat, max_lat):
        plt.figure()
        axes = plt.axes(projection=projection)
        axes.set_extent(
            (-180, 180, min_lat, max_lat),
            crs=ccrs.PlateCarree()
        )
        for diff_month in diff_cube.slices_over('month_number'):
            month_number = diff_month.coord('month_number').points[0]
            month_name = calendar.month_name[month_number]
            logger.info('Plotting 2D blocking for ' + month_name)
            diff_month.long_name += ' (' + month_name.title() + ')'
            iris.quickplot.pcolormesh(
                diff_month,
                coords=('longitude', 'latitude'),
                cmap=cmap,
                vmin=-self.max_color_scale,
                vmax=self.max_color_scale
            )
            axes.coastlines()
            axes.gridlines(alpha=0.5, linestyle='--')
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            axes.set_boundary(circle, transform=axes.transAxes)
            plt.savefig(self._get_plot_name(
                'blocking2Ddiff',
                dataset,
                month_number
            ))
            plt.close()

    def create_comparison_plot(self, datasets, metrics, name):
        logger.debug(metrics)
        sdev = [metrics[datasets[0]]['sdev'][0]] + \
            [m['sdev'][1] for m in metrics.values()]
        crmsd = [metrics[datasets[0]]['crmsd'][0]] + \
            [m['crmsd'][1] for m in metrics.values()]
        ccoef = [metrics[datasets[0]]['ccoef'][0]] + \
            [m['ccoef'][1] for m in metrics.values()]
        skill_metrics.taylor_diagram(
            np.array(sdev),
            np.array(crmsd),
            np.array(ccoef),
            styleOBS='-',
            colOBS='r',
            markerobs='o',
            titleOBS='reference',
        )
        out_type = self.cfg[n.OUTPUT_FILE_TYPE]
        name = '{}.{}'.format(name, out_type)
        plt.savefig(os.path.join(self.cfg[n.PLOT_DIR], name))
        plt.close()

    def _create_dataset_legend(self, datasets):
        handles = []
        for num, filename in enumerate(datasets):
            handles.append(
                Line2D(
                    [0], [0],
                    marker=self._get_marker(num),
                    color='white',
                    markerfacecolor='#000000', markeredgecolor='#000000',
                    label=self.datasets.get_info(n.DATASET, filename)
                )
            )

        plt.legend(
            handles=handles,
            loc='upper left',
            bbox_to_anchor=(1, 1),
            ncol=1,
            frameon=False,
        )

    def _get_marker(self, num):
        return '$\\mathrm{{{0}}}$'.format(chr(num + ord('A')))


def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        Blocking(config).compute()


if __name__ == '__main__':
    main()
