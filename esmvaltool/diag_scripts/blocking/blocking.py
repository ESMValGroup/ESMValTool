"""Blocking diagnostic"""
import os
import logging
import itertools
import calendar
import math

import numpy as np

import matplotlib
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import colors
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter

import iris
import iris.time
import iris.util
import iris.coord_categorisation
import iris.analysis
import iris.analysis.stats
import iris.coords
import iris.quickplot
import cartopy.crs as ccrs

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.preprocessor import regrid

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

        def _get_index(self, high, central, low,
                       north_distance, south_distance):
            if ((high - central) / north_distance) > self.north_threshold:
                return 0
            if ((central - low) / south_distance) < self.south_threshold:
                return 0
            return 1
        self._compute_index = np.vectorize(_get_index, [np.int8])

    def compute(self):
        """Compute blocking diagnostic"""
        logger.info('Computing blocking')

        logger.info('Reference dataset %s', self.reference_dataset)
        reference_1d, reference_2d = \
            self._get_blocking_indices(self.reference_dataset)

        error = {}
        correlation = {}

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

                for diff_slice in diff.slices_over('month_number'):
                    self.plot_differences(
                        filename, diff_slice, cmap, projection,
                        min_lat, max_lat
                    )

                rms = diff.collapsed(
                    ('longitude', 'latitude'),
                    iris.analysis.RMS
                )

                if self.cfg[n.WRITE_NETCDF]:
                    new_filename = os.path.basename(filename).replace(
                        'zg', 'blocking2Drms')
                    netcdf_path = os.path.join(
                        self.cfg[n.WORK_DIR], new_filename)
                    iris.save(rms, target=netcdf_path, zlib=True)

                corr = iris.analysis.stats.pearsonr(
                    reference_2d, dataset_2d,
                    corr_coords=('longitude', 'latitude'),
                    weights=iris.analysis.cartography.area_weights(
                        reference_2d),
                )
                if self.cfg[n.WRITE_NETCDF]:
                    new_filename = os.path.basename(filename).replace(
                        'zg', 'blocking2Dcorr')
                    netcdf_path = os.path.join(
                        self.cfg[n.WORK_DIR], new_filename)
                    iris.save(corr, target=netcdf_path, zlib=True)

                error[filename] = rms
                correlation[filename] = corr
                logger.info('Correlation: %f', corr.data)
                logger.info('RMSE: %f', rms.data)

        if self.cfg[n.WRITE_PLOTS]:
            self.create_comparison_plot(datasets, correlation, error)
            for month in range(1, 13):
                self.create_comparison_plot(
                    datasets, correlation, error, month
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
        self.latitude_data = {}
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
        result.remove_coord('time')
        iris.util.promote_aux_coord_to_dim_coord(result, 'month_number')
        return result

    def _get_plot_name(self, name, filename, month=None):
        dataset = self.datasets.get_info(n.DATASET, filename)
        project = self.datasets.get_info(n.PROJECT, filename)
        ensemble = self.datasets.get_info(n.ENSEMBLE, filename)
        start = self.datasets.get_info(n.START_YEAR, filename)
        end = self.datasets.get_info(n.END_YEAR, filename)

        plot_path = os.path.join(
            self.cfg[n.PLOT_DIR],
            project, dataset)
        if ensemble is not None:
            plot_path = os.path.join(plot_path, ensemble)
        if not os.path.isdir(plot_path):
            os.makedirs(plot_path)

        if ensemble is None:
            ensemble = ''
        else:
            ensemble += '_'
        if month is not None:
            name = '{}_{:02}'.format(name, month)
        out_type = self.cfg[n.OUTPUT_FILE_TYPE]

        plot_filename = '{name}_{project}_{dataset}_' \
                        '{ensemble}{start}-{end}' \
                        '.{out_type}'.format(
                            name=name,
                            dataset=dataset,
                            project=project,
                            ensemble=ensemble,
                            start=start,
                            end=end,
                            out_type=out_type)

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
            self, zg_high, zg_central, zg_low,
            north_distance, south_distance)

        del self.latitude_data[low_lat]

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
        if latitude not in self.latitude_data:
            lat_data = zg500.extract(iris.Constraint(latitude=latitude)).data
            self.latitude_data[latitude] = lat_data
        else:
            lat_data = self.latitude_data[latitude]
        return lat_data

    def _apply_persistence(self, blocking_cube):
        for lon_slice in blocking_cube.slices('longitude'):
            grouped = ((k, sum(1 for _ in g))
                       for k, g in itertools.groupby(lon_slice.data))
            index = 0
            for value, length in grouped:
                if value and length < self.persistence:
                    lon_slice.data[index: index + length] = False
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

    def plot_differences(self, filename, diff_cube, cmap, projection,
                         min_lat, max_lat):
        plt.figure()
        axes = plt.axes(projection=projection)
        axes.set_extent(
            (-180, 180, min_lat, max_lat),
            crs=ccrs.PlateCarree()
        )
        iris.quickplot.pcolormesh(
            diff_cube,
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
            filename,
            diff_cube.coord('month_number').points[0]))
        plt.close()

    def create_comparison_plot(self, datasets, correlation, error,
                               month=None):
        plt.figure()
        ax = plt.gca()
        if month:
            color = 'black'
        else:
            color = [
                '#0000ff', '#7080ff',
                '#00FF00', '#00CC00', '#00AA00',
                '#FFD700', '#CCCC00', '#AAAA00',
                '#A0522D', '#8B0000', '#8B4513',
                '#000090',
                ]

        for num, filename in enumerate(datasets):
            corr = correlation[filename]
            err = error[filename]

            if month is not None:
                corr = corr.extract(iris.Constraint(month_number=month))
                err = err.extract(iris.Constraint(month_number=month))

            plt.scatter(
                corr.data,
                err.data,
                c=color,
                marker=self._get_marker(num),
                zorder=2,
            )

        box = ax.get_position()
        if month is None:
            ax.set_position(
                [box.x0, box.y0 + box.height * 0.20,
                 box.width * 0.80, box.height * 0.80]
            )
        else:
            ax.set_position(
                [box.x0, box.y0,
                 box.width * 0.80, box.height]
            )
        if month is None:
            ax.set_title('Blocking 2D')
        else:
            ax.set_title('Blocking 2D ({})'.format(calendar.month_name[month]))
        ax.set_xlabel('Pearson correlation')
        ax.set_ylabel('Root Mean Square Error (days per month)')
        ax.set_xticks([0], minor=False)
        ax.set_xticks(np.arange(-1, 1.1, 0.25), minor=True)
        ax.tick_params(axis='x', which='major', labelsize=0)
        ax.xaxis.set_minor_formatter(FormatStrFormatter("%.2f"))
        top = math.ceil(ax.get_ylim()[1]) + 1
        ax.set_yticks(np.arange(0, top, 1), minor=False)
        ax.set_yticks(np.arange(0, top - 0.5, 0.5), minor=True)
        ax.grid(True, 'major', linestyle='-', color='black', zorder=0)
        ax.grid(True, 'minor', linestyle=':', color='black', zorder=0)
        plt.ylim(ymin=0)
        plt.xlim(-1, 1)
        if month is None:
            legend = plt.legend(
                handles=[
                    Patch(facecolor=col, label=calendar.month_name[num + 1])
                    for num, col in enumerate(color)
                ],
                loc='upper center',
                bbox_to_anchor=(0.5, -0.15),
                ncol=4,
                frameon=False,
            )
        self._create_dataset_legend(datasets)
        if month is None:
            ax.add_artist(legend)
        out_type = self.cfg[n.OUTPUT_FILE_TYPE]
        if month is None:
            name = 'blocking2D.{}'.format(out_type)
        else:
            name = 'blocking2D_{:02}.{}'.format(month, out_type)
        plt.savefig(os.path.join(
            self.cfg[n.PLOT_DIR],
            name,
        ))
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
