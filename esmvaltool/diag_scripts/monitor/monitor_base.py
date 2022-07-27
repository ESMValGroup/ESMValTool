"""Base class for monitoring diagnostics."""

import logging
import os

import cartopy
import matplotlib.pyplot as plt
import yaml
from esmvalcore._data_finder import _replace_tags
from iris.analysis import MEAN
from mapgenerator.plotting.timeseries import PlotSeries

from esmvaltool.diag_scripts.shared import ProvenanceLogger, names

logger = logging.getLogger(__name__)


class MonitorBase():
    """Base class for monitoring diagnostic.

    It contains the common methods for path creation, provenance
    recording, option parsing and to create some common plots.

    """

    def __init__(self, config):
        self.cfg = config
        plot_folder = config.get(
            'plot_folder',
            '{plot_dir}/../../{dataset}/{exp}/{modeling_realm}/{real_name}',
        )
        plot_folder = plot_folder.replace('{plot_dir}',
                                          self.cfg[names.PLOT_DIR])
        self.plot_folder = os.path.abspath(plot_folder)
        self.plot_filename = config.get(
            'plot_filename',
            '{plot_type}_{real_name}_{dataset}_{mip}_{exp}_{ensemble}')
        self.plots = config.get('plots', {})
        default_config = os.path.join(os.path.dirname(__file__),
                                      "monitor_config.yml")
        cartopy_data_dir = config.get('cartopy_data_dir', None)
        if cartopy_data_dir:
            cartopy.config['data_dir'] = cartopy_data_dir
        with open(config.get('config_file', default_config)) as config_file:
            self.config = yaml.safe_load(config_file)

    def _add_file_extension(self, filename):
        """Add extension to plot filename."""
        return f"{filename}.{self.cfg['output_file_type']}"

    def _get_proj_options(self, map_name):
        return self.config['maps'][map_name]

    def _get_variable_options(self, variable_group, map_name):
        options = self.config['variables'].get(
            variable_group, self.config['variables']['default'])
        if 'default' not in options:
            variable_options = options
        else:
            variable_options = options['default']
            if map_name in options:
                variable_options = {**variable_options, **options[map_name]}

        if 'bounds' in variable_options:
            if not isinstance(variable_options['bounds'], str):
                variable_options['bounds'] = [
                    float(n) for n in variable_options['bounds']
                ]
        logger.debug(variable_options)
        return variable_options

    def plot_timeseries(self, cube, var_info, period='', **kwargs):
        """Plot timeseries from a cube.

        It also automatically smoothes it for long timeseries of monthly data:
            - Between 10 and 70 years long, it also plots the 12-month rolling
              average along the raw series
            - For more than ten years, it plots the 12-month and 10-years
              rolling averages and not the raw series

        """
        if 'xlimits' not in kwargs:
            kwargs['xlimits'] = 'auto'
        length = cube.coord("year").points.max() - cube.coord(
            "year").points.min()
        filename = self.get_plot_path(f'timeseries{period}', var_info,
                                      add_ext=False)
        caption = ("{} of "
                   f"{var_info[names.LONG_NAME]} of dataset "
                   f"{var_info[names.DATASET]} (project "
                   f"{var_info[names.PROJECT]}) from "
                   f"{var_info[names.START_YEAR]} to "
                   f"{var_info[names.END_YEAR]}.")
        if length < 10 or length * 11 > cube.coord("year").shape[0]:
            self.plot_cube(cube, filename, **kwargs)
            self.record_plot_provenance(
                self._add_file_extension(filename),
                var_info,
                'timeseries',
                period=period,
                caption=caption.format("Time series"),
            )
        elif length < 70:
            self.plot_cube(cube, filename, **kwargs)
            self.record_plot_provenance(
                self._add_file_extension(filename),
                var_info,
                'timeseries',
                period=period,
                caption=caption.format("Time series"),
            )

            # Smoothed time series (12-month running mean)
            plt.gca().set_prop_cycle(None)
            self.plot_cube(cube.rolling_window('time', MEAN, 12),
                           f"{filename}_smoothed_12_months",
                           **kwargs)
            self.record_plot_provenance(
                self._add_file_extension(f"{filename}_smoothed_12_months"),
                var_info,
                'timeseries',
                period=period,
                caption=caption.format(
                    "Smoothed (12-months running mean) time series"),
            )
        else:
            # Smoothed time series (12-month running mean)
            self.plot_cube(cube.rolling_window('time', MEAN, 12),
                           f"{filename}_smoothed_12_months",
                           **kwargs)
            self.record_plot_provenance(
                self._add_file_extension(f"{filename}_smoothed_12_months"),
                var_info,
                'timeseries',
                period=period,
                caption=caption.format(
                    "Smoothed (12-months running mean) time series"),
            )

            # Smoothed time series (10-year running mean)
            self.plot_cube(cube.rolling_window('time', MEAN, 120),
                           f"{filename}_smoothed_10_years",
                           **kwargs)
            self.record_plot_provenance(
                self._add_file_extension(f"{filename}_smoothed_10_years"),
                var_info,
                'timeseries',
                period=period,
                caption=caption.format(
                    "Smoothed (10-years running mean) time series"),
            )

    def record_plot_provenance(self, filename, var_info, plot_type, **kwargs):
        """Write provenance info for a given file."""
        with ProvenanceLogger(self.cfg) as provenance_logger:
            prov = self.get_provenance_record(
                ancestor_files=[var_info['filename']],
                plot_type=plot_type,
                **kwargs,
            )
            provenance_logger.log(filename, prov)

    def plot_cube(self, cube, filename, linestyle='-', **kwargs):
        """Plot a timeseries from a cube.

        Supports multiplot layouts for cubes with extra dimensions
        `shape_id` or `region`.

        """
        plotter = PlotSeries()
        plotter.filefmt = self.cfg['output_file_type']
        plotter.img_template = filename
        region_coords = ('shape_id', 'region')

        for region_coord in region_coords:
            if cube.coords(region_coord):
                if cube.coord(region_coord).shape[0] > 1:
                    plotter.multiplot_cube(cube, 'time', region_coord,
                                           **kwargs)
                    return
        plotter.plot_cube(cube, 'time', linestyle=linestyle, **kwargs)

    @staticmethod
    def get_provenance_record(ancestor_files, **kwargs):
        """Create provenance record for the diagnostic data and plots."""
        record = {
            'authors': [
                'vegas-regidor_javier',
            ],
            'references': [
                'acknow_project',
            ],
            'ancestors': ancestor_files,
            **kwargs
        }
        return record

    def get_plot_path(self, plot_type, var_info, add_ext=True):
        """Get plot full path from variable info.

        Parameters
        ----------
        plot_type: str
            Name of the plot
        var_info: dict
            Variable information from ESMValTool
        add_ext: bool, optional (default: True)
            Add filename extension from configuration file.

        """
        return os.path.join(
            self.get_plot_folder(var_info),
            self.get_plot_name(plot_type, var_info, add_ext=add_ext),
        )

    def get_plot_folder(self, var_info):
        """Get plot storage folder from variable info.

        Parameters
        ----------
        var_info: dict
            Variable information from ESMValTool

        """
        info = {
            'real_name': self._real_name(var_info['variable_group']),
            **var_info
        }
        folder = os.path.expandvars(
            os.path.expanduser(
                list(_replace_tags(self.plot_folder, info))[0]
            )
        )
        if self.plot_folder.startswith('/'):
            folder = '/' + folder
        if not os.path.isdir(folder):
            os.makedirs(folder, exist_ok=True)
        return folder

    def get_plot_name(self, plot_type, var_info, add_ext=True):
        """Get plot filename from variable info.

        Parameters
        ----------
        plot_type: str
            Name of the plot
        var_info: dict
            Variable information from ESMValTool
        add_ext: bool, optional (default: True)
            Add filename extension from configuration file.

        """
        info = {
            "plot_type": plot_type,
            'real_name': self._real_name(var_info['variable_group']),
            **var_info
        }
        file_name = list(_replace_tags(self.plot_filename, info))[0]
        if add_ext:
            file_name = self._add_file_extension(file_name)
        return file_name

    @staticmethod
    def _set_rasterized(axes=None):
        """Rasterize all artists and collection of axes if desired."""
        if axes is None:
            axes = plt.gca()
        if not isinstance(axes, list):
            axes = [axes]
        for single_axes in axes:
            for artist in single_axes.artists:
                artist.set_rasterized(True)
            for collection in single_axes.collections:
                collection.set_rasterized(True)

    @staticmethod
    def _real_name(variable_group):
        for subfix in ('Ymean', 'Ysum', 'mean', 'sum'):
            if variable_group.endswith(subfix):
                variable_group = variable_group.replace(subfix, '')
        return variable_group
