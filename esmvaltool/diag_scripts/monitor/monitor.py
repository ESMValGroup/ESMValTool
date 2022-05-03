"""Diagnostic to plot preprocessor output."""

import calendar
import logging

import iris
import iris.coord_categorisation
import matplotlib.pyplot as plt
import numpy as np
from esmvalcore.preprocessor import climate_statistics
from iris.coords import AuxCoord
from mapgenerator.plotting.plotmap import PlotMap
from mapgenerator.plotting.timeseries import PlotSeries

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.monitor.monitor_base import MonitorBase
from esmvaltool.diag_scripts.shared import group_metadata

logger = logging.getLogger(__name__)


class Monitor(MonitorBase):
    """Diagnostic to plot preprocessor output."""

    def __init__(self, config):
        super().__init__(config)
        self.plots = config.get('plots', {})
        self.has_errors = False

    def compute(self):
        """Plot preprocessed data."""
        for module in ['matplotlib', 'fiona']:
            module_logger = logging.getLogger(module)
            module_logger.setLevel(logging.WARNING)
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        for alias in data:
            variables = group_metadata(data[alias], 'variable_group')

            for var_name, var_info in variables.items():
                logger.info('Plotting variable %s', var_name)
                var_info = var_info[0]
                cubes = iris.load(var_info['filename'])
                if len(cubes) == 1:
                    cube = cubes[0]
                else:
                    for cube in cubes:
                        if cube.var_name == var_name:
                            break
                    else:
                        raise ValueError(
                            f'Can not find cube {var_name} in {cubes}')
                cube.var_name = self._real_name(var_name)
                cube.attributes['plot_name'] = var_info.get('plot_name', '')

                self.timeseries(cube, var_info)
                self.plot_annual_cycle(cube, var_info)
                self.plot_monthly_climatology(cube, var_info)
                self.plot_seasonal_climatology(cube, var_info)
                self.plot_climatology(cube, var_info)
        if self.has_errors:
            raise Exception(
                'Errors detected. Please check log for more details')

    @staticmethod
    def _add_month_name(cube):
        if cube.coords('month_number'):
            month_number = cube.coord('month_number')
            points = np.empty(month_number.shape, dtype='|S12')
            for i in range(1, 13):
                points[month_number.points == i] = calendar.month_name[i]
            cube.add_aux_coord(
                AuxCoord(points=points,
                         var_name='month_name',
                         long_name='month_name'),
                cube.coord_dims(month_number))
            points = np.empty(month_number.shape, dtype='|S3')
            for i in range(1, 13):
                points[month_number.points == i] = str(
                    calendar.month_name[i].upper())
            cube.add_aux_coord(
                AuxCoord(points=points, var_name='month', long_name='month'),
                cube.coord_dims(month_number))
            return

    def timeseries(self, cube, var_info):
        """Plot timeseries according to configuration.

        The key 'timeseries' must be passed to the 'plots' option in the
        configuration.

        Parameters
        ----------
        cube: iris.cube.Cube
            Data to plot. Must be 1D with time or 2D with an extra 'shape_id'
            or 'region' coordinate. In that case, the plot will be a multiple
            one with one figure for each region
        var_info: dict
            Variable's metadata from ESMValTool
        """
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
                                period='start',
                                suptitle='First 50 years')
            self.plot_timeseries(cube.extract(
                iris.Constraint(
                    year=lambda cell: cell >= (var_info[n.END_YEAR] - 50))),
                                 var_info,
                                 period='end',
                                 suptitle='Last 50 years')

    def plot_annual_cycle(self, cube, var_info):
        """Plot the annual cycle according to configuration.

        The key 'annual_cycle' must be passed to the 'plots' option in the
        configuration.

        Parameters
        ----------
        cube: iris.cube.Cube
            Data to plot. Must be 1D with time or 2D with an extra 'shape_id'
            or 'region' coordinate. In that case, the plot will be a multiple
            one with one figure for each region
        var_info: dict
            Variable's metadata from ESMValTool

        Warning
        -------
        The monthly climatology is done inside the function so the users can
        plot both the timeseries and the annual cycle in one go
        """
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
        self.record_plot_provenance(
            self.get_plot_path('annualcycle', var_info, 'svg'),
            var_info,
            'Annual cycle',
        )

    def plot_monthly_climatology(self, cube, var_info):
        """Plot the monthly climatology as a multipanel plot.

        The key 'monclim' must be passed to the 'plots' option in the
        configuration.

        Parameters
        ----------
        cube: iris.cube.Cube
            Data to plot. Must be 3D with latitude, longitude and month_number
        var_info: dict
            Variable's metadata from ESMValTool
        """
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
                if cube_slice.ndim != 2:
                    logger.error(
                        'Climatologies can only be plotted for 2D vars. '
                        'Skipping...')
                    self.has_errors = True
                    return
                self._plot_monthly_cube(plot_map, months, columns, rows,
                                        map_options, variable_options,
                                        cube_slice)
            plt.suptitle(
                'Monthly climatology '
                f'({var_info[n.START_YEAR]}-{var_info[n.END_YEAR]})'
                f'\n{cube.long_name} ({cube.units})',
                fontsize=plot_map.fontsize + 4.,
                y=1.025 - rows * 0.025,
            )
            plt.subplots_adjust(
                top=0.85,
                bottom=.05,
                left=0,
                right=1,
                hspace=.20,
                wspace=.15,
            )
            filename = self.get_plot_path(f'monclim{map_name}',
                                          var_info,
                                          file_type='png')
            plt.savefig(
                filename,
                bbox_inches='tight',
                pad_inches=.2,
            )
            plt.close(plt.gcf())
            self.record_plot_provenance(
                filename,
                var_info,
                'Monthly climatology',
                region=map_name,
            )
        cube.remove_coord('month')
        cube.remove_coord('month_name')

    @staticmethod
    def _plot_monthly_cube(plot_map, months, columns, rows, map_options,
                           variable_options, cube_slice):
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

    def plot_seasonal_climatology(self, cube, var_info):
        """Plot the seasonal climatology as a multipanel plot.

        The key 'seasonclim' must be passed to the 'plots' option in the
        configuration.

        Parameters
        ----------
        cube: iris.cube.Cube
            Data to plot. Must be 3D with latitude, longitude and month_number
            or season
        var_info: dict
            Variable's metadata from ESMValTool

        Warning
        -------
        The seasonal climatology can be done inside the function so the users
        can plot monthly, seasonal and yearly climatologies in one go
        """
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
        if cube.coords('month_number'):
            points = [
                season[point] for point in cube.coord('month_number').points
            ]
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
                if cube_slice.ndim != 2:
                    logger.error(
                        'Climatologies can only be plotted for 2D vars. '
                        'Skipping...')
                    self.has_errors = True
                    return
                plot_map.plot_cube(
                    cube_slice,
                    save=False,
                    subplot=(2, 2, index),
                    keep_aspect=True,
                    title=season,
                    **{
                        **map_options,
                        **variable_options,
                    },
                )
            plt.tight_layout()
            plt.suptitle(
                'Seasonal climatology  '
                f'({var_info[n.START_YEAR]}-{var_info[n.END_YEAR]})\n'
                f'{cube.long_name} ({cube.units})',
                fontsize=plot_map.fontsize + 4,
            )
            plt.subplots_adjust(
                top=0.85,
                bottom=.05,
                left=0,
                right=1,
                hspace=.20,
                wspace=.15,
            )
            filename = self.get_plot_path(f'seasonclim{map_name}', var_info,
                                          'png')
            plt.savefig(
                filename,
                bbox_inches='tight',
                pad_inches=.2,
            )
            plt.close(plt.gcf())
            self.record_plot_provenance(
                filename,
                var_info,
                'Seasonal climatology',
                region=map_name,
            )
        cube.remove_coord('season')

    def plot_climatology(self, cube, var_info):
        """Plot the climatology as a multipanel plot.

        The key 'clim' must be passed to the 'plots' option in the
        configuration.

        Parameters
        ----------
        cube: iris.cube.Cube
            Data to plot. Must be 3D with latitude, longitude and month_number
            or season or 2D with latitude and longitude
        var_info: dict
            Variable's metadata from ESMValTool

        Warning
        -------
        The climatology can be done inside the function from the monthly and
        seasonal climatologies so the users can plot several of them in one go
        """
        if 'clim' not in self.plots:
            return

        if cube.coords('month_number'):
            cube = cube.collapsed('month_number', iris.analysis.MEAN)
        elif cube.coords('season'):
            cube = cube.collapsed('season', iris.analysis.MEAN)
        maps = self.plots['clim'].get('maps', ['global'])
        plot_map = PlotMap(loglevel='INFO')
        plot_map.outdir = self.get_plot_folder(var_info)
        for map_name in maps:
            map_options = self._get_proj_options(map_name)

            variable_options = self._get_variable_options(
                var_info['variable_group'], map_name)
            if cube.ndim != 2:
                logger.error('Climatologies can only be plotted for 2D vars. '
                             'Skipping...')
                self.has_errors = True
                return

            plot_map.plot_cube(cube,
                               save=False,
                               **{
                                   **map_options,
                                   **variable_options
                               })
            plt.suptitle(
                f'Climatology ({var_info[n.START_YEAR]}'
                f'-{var_info[n.END_YEAR]})',
                y=map_options.get('suptitle_pos', 0.95),
                fontsize=plot_map.fontsize + 4)
            filename = self.get_plot_path(f'clim{map_name}', var_info, 'png')
            plt.savefig(filename,
                        bbox_inches='tight',
                        pad_inches=.2,
                        dpi=plot_map.dpi)
            plt.close(plt.gcf())
            self.record_plot_provenance(
                filename,
                var_info,
                'Climatology',
                region=map_name,
            )


def main():
    """Execute diagnostic."""
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        Monitor(config).compute()


if __name__ == "__main__":
    main()
