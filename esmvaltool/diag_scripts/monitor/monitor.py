import calendar
import logging
import os

import iris
import iris.coord_categorisation
import matplotlib.pyplot as plt
import numpy as np
from esmvalcore.preprocessor import climate_statistics
from iris.coords import AuxCoord
from mapgenerator.plotting.plotmap import PlotMap
from mapgenerator.plotting.timeseries import PlotSeries
from esmvaltool.diag_scripts.monitor.monitor_base import MonitorBase

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import group_metadata

logger = logging.getLogger(__name__)


class Monitor(MonitorBase):
    def __init__(self, config):
        super().__init__(config)
        self.plots = config.get('plots', {})

    def compute(self):
        for module in ['matplotlib', 'fiona']:
            module_logger = logging.getLogger(module)
            module_logger.setLevel(logging.WARNING)
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        for alias in data:
            variables = group_metadata(data[alias], 'variable_group')

            for var_name, var_info in variables.items():
                logger.info('Plotting variable %s', var_name)
                var_info = var_info[0]
                cube = iris.load_cube(var_info['filename'])
                cube.var_name = self._real_name(var_name)
                cube.attributes['plot_name'] = var_info.get('plot_name', '')

                self.timeseries(cube, var_info)
                self.plot_annual_cycle(cube, var_info)
                self.plot_monthly_climatology(cube, var_info)
                self.plot_seasonal_climatology(cube, var_info)
                self.plot_climatology(cube, var_info)
                self.plot_vertical_profile(cube, var_info)
                self.plot_zonal(cube, var_info)

    @staticmethod
    def _add_month_name(cube):
        if cube.coords('month_number'):
            month_number = cube.coord('month_number')
            points = np.empty(month_number.shape, dtype='|S12')
            for m in range(1, 13):
                points[month_number.points == m] = calendar.month_name[m]
            cube.add_aux_coord(
                AuxCoord(points=points,
                         var_name='month_name',
                         long_name='month_name'),
                cube.coord_dims(month_number))
            points = np.empty(month_number.shape, dtype='|S3')
            for m in range(1, 13):
                points[month_number.points ==
                       m] = calendar.month_name[m].upper()
            cube.add_aux_coord(
                AuxCoord(points=points, var_name='month', long_name='month'),
                cube.coord_dims(month_number))
            return

    def timeseries(self, cube, var_info):
        if 'timeseries' not in self.plots:
            return
        if not cube.coords('year'):
            iris.coord_categorisation.add_year(cube, 'time')
        self.plot_timeseries(cube, var_info, suptitle='Full period')
        if var_info[n.END_YEAR] - var_info[n.START_YEAR] > 75:
            self.plot_timeseries(cube.extract(
                iris.Constraint(
                    year=lambda cell: cell <= (var_info[n.START_YEAR] + 50))),
                                 var_info,
                                 type='start',
                                 suptitle='First 50 years')
            self.plot_timeseries(cube.extract(
                iris.Constraint(
                    year=lambda cell: cell >= (var_info[n.END_YEAR] - 50))),
                                 var_info,
                                 type='end',
                                 suptitle='Last 50 years')

    def plot_annual_cycle(self, cube, var_info):
        if 'annual_cycle' not in self.plots:
            return
        cube = climate_statistics(cube, period='month')
        self._add_month_name(cube)

        plotter = PlotSeries()
        plotter.outdir = self.get_plot_folder(var_info)
        plotter.img_template = self.get_plot_name('annualcycle', var_info,
                                                  None)
        plotter.filefmt = 'svg'
        region_coords = ('shape_id', 'region')
        options = {
            'xlabel': '',
            'xlimits': None,
            'suptitle': 'Annual cycle',
        }
        for region_coord in region_coords:
            if cube.coords(region_coord):
                plotter.multiplot_cube(cube, 'month', region_coord, **options)
                return
        plotter.plot_cube(cube, 'month', **options)

    def plot_monthly_climatology(self, cube, var_info):
        if 'monclim' not in self.plots:
            return

        plot_map = PlotMap()
        maps = self.plots['monclim'].get('maps', ['global'])
        months = self.plots['monclim'].get('months', None)
        plot_size = self.plots['monclim'].get('plot_size', (5, 4))
        columns = self.plots['monclim'].get('columns', 3)
        rows = self.plots['monclim'].get('rows', 4)
        if months:
            cube = cube.extract(
                iris.Constraint(month_number=lambda cell: cell in months))
        self._add_month_name(cube)
        for map_name in maps:
            map_options = self._get_proj_options(map_name)
            variable_options = self._get_variable_options(
                var_info['variable_group'], map_name)
            plt.figure(figsize=(plot_size[0] * columns, plot_size[1] * rows),
                       dpi=120)

            for cube_slice in cube.slices_over('month_number'):
                month = cube_slice.coord('month_number').points[0]
                month_name = cube_slice.coord('month_name').points[0]
                if months:
                    index = months.index(month) + 1
                else:
                    if month == 12:
                        index = 1
                    else:
                        index = month + 1
                plot_map.plot_cube(
                    cube_slice,
                    save=False,
                    subplot=(rows, columns, index),
                    keep_aspect=True,
                    title=month_name.decode(),
                    **{
                        **map_options,
                        **variable_options
                    },
                )
            plt.suptitle('{0.long_name} ({0.units})'.format(cube),
                         fontsize=plot_map.fontsize + 6.,
                         y=1.02)
            plt.savefig(
                self.get_plot_path(f'monclim{map_name}',
                                   var_info,
                                   file_type='png'),
                bbox_inches='tight',
                pad_inches=.2,
            )

    def plot_seasonal_climatology(self, cube, var_info):
        if 'seasonclim' not in self.plots:
            return

        season = {
            12: 'DJF',
            1: 'DJF',
            2: 'DJF',
            3: 'MAM',
            4: 'MAM',
            5: 'MAM',
            6: 'JJA',
            7: 'JJA',
            8: 'JJA',
            9: 'SON',
            10: 'SON',
            11: 'SON'
        }
        points = [season[point] for point in cube.coord('month_number').points]
        cube.add_aux_coord(iris.coords.AuxCoord(points, var_name='season'),
                           cube.coord_dims('month_number'))
        cube = cube.aggregated_by('season', iris.analysis.MEAN)

        plot_map = PlotMap()
        maps = self.plots['seasonclim'].get('maps', ['global'])
        for map_name in maps:
            map_options = self._get_proj_options(map_name)
            variable_options = self._get_variable_options(
                var_info['variable_group'], map_name)
            index = 0
            for cube_slice in cube.slices_over('season'):
                index += 1
                season = cube_slice.coord('season').points[0]
                plot_map.plot_cube(
                    cube_slice,
                    save=False,
                    subplot=(2, 2, index),
                    keep_aspect=True,
                    title=f"Season {season}",
                    **{
                        **map_options,
                        **variable_options
                    },
                )
            plt.tight_layout()
            plt.suptitle(
                'Seasonal climatology\n{0.long_name} ({0.units})'.format(cube),
                fontsize=plot_map.fontsize * 2)
            plt.subplots_adjust(bottom=.05, hspace=.3, left=.1, right=1)
            plt.savefig(
                self.get_plot_path(f'seasonclim{map_name}', var_info, 'png'),
                bbox_inches='tight',
                pad_inches=.2,
            )

    def plot_climatology(self, cube, var_info):
        if 'clim' not in self.plots:
            return

        cube = cube.collapsed('month_number', iris.analysis.MEAN)
        maps = self.plots['clim'].get('maps', ['global'])
        plot_map = PlotMap(loglevel='INFO')
        plot_map.outdir = self.get_plot_folder(var_info)
        for map_name in maps:
            map_options = self._get_proj_options(map_name)

            variable_options = self._get_variable_options(
                var_info['variable_group'], map_name)

            plot_map.plot_cube(cube,
                               save=False,
                               **{
                                   **map_options,
                                   **variable_options
                               })
            plt.suptitle('Climatology',
                         y=map_options.get('suptitle_pos', 0.9),
                         fontsize=plot_map.fontsize + 4)
            plt.savefig(self.get_plot_path(f'clim{map_name}', var_info, 'png'),
                        bbox_inches='tight',
                        pad_inches=.2,
                        dpi=plot_map.dpi)
            plt.close(plt.gcf())

    def plot_vertical_profile(self, cube, var_info):
        if 'vertical_profile' not in self.plots:
            return
        self._add_month_name(cube)
        plotter = PlotSeries(scaley="log")
        plotter.filefmt = 'svg'
        filename = self.get_plot_path('vertprofile', var_info, None)
        plotter.img_template = filename
        plotter.plot_size = (4., 8.)
        plotter.multiplot_cube(
            cube,
            'depth',
            'month_name',
            ncols=3,
            invert=True,
            inverty=True,
            sharex=True,
            sharey=True,
            suptitle='Monthly climatology vertical profiles')

    def plot_zonal(self, cube, var_info):
        if 'zonal' not in self.plots:
            return
        filename = self.get_plot_path('zonal', var_info, None)
        for basin_slice in cube.slices_over('region'):
            logger.info(basin_slice)
            plt.pcolormesh(
                basin_slice.coord('grid_latitude').points,
                basin_slice.coord('depth').points,
                basin_slice.data,
            )
            plt.title('{0.long_name} ({0.units})'.format(cube), fontsize=14)
            plt.suptitle('Zonal profile climatology', fontsize=14)
        plt.subplots_adjust(bottom=.05, hspace=.3, left=.1, right=1)

        plt.savefig(filename, bbox_inches='tight', pad_inches=.2)


def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        Monitor(config).compute()


if __name__ == "__main__":
    main()
