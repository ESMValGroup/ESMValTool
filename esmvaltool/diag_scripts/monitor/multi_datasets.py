#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Monitoring diagnostic to show multiple datasets in one plot (incl. biases).

Description
-----------
This diagnostic can be used to visualize multiple datasets in one plot.

For some plot types, a reference dataset can be defined. For this, use the
facet ``reference_for_monitor_diags: true`` in the definition of the dataset in
the recipe. Note that at most one reference dataset per variable is supported.

Currently supported plot types (use the option ``plots`` to specify them):
    - Time series (data needs to be 1D with single dimension 'time'): for each
      variable separately, all datasets are plotted in one single figure.
    - Maps (data needs to be 2D with dimensions 'latitude', 'longitude'): for
      each variable and dataset, an individual map is plotted. If a reference
      dataset is defined, also include this dataset and a bias plot into the
      figure. Note that if a reference dataset is defined, all input datasets
      need to be given on the same horizontal grid (you can use the
      preprocessor :func:`esmvalcore.preprocessor.regrid` for this).
    - Profiles (data needs to be 2D with dimensions 'latitude',
      'height'/'air_pressure'): for each variable and dataset, an individual
      profile is plotted. If a reference dataset is defined, also include this
      dataset and a bias plot into the figure. Note that if a reference dataset
      is defined, all input datasets need to be given on the same horizontal
      and vertical grid (you can use the preprocessors
      :func:`esmvalcore.preprocessor.regrid` and
      :func:`esmvalcore.preprocessor.extract_levels` for this).

Author
------
Manuel Schlund (DLR, Germany)

Configuration options in recipe
-------------------------------
facet_used_for_labels: str, optional (default: 'dataset')
    Facet used to label different datasets in plot titles and legends. For
    example, ``facet_used_for_labels: dataset`` will use dataset names in plot
    titles and legends; ``facet_used_for_labels: exp`` will use experiments in
    plot titles and legends. In addition, ``facet_used_for_labels`` is used to
    select the correct ``plot_kwargs`` for the different datasets (see
    configuration options for the different plot types below).
figure_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.figure`.
plots: dict, optional
    Plot types plotted by this diagnostic. Dictionary keys must be
    ``timeseries``, ``map``, or ``profile``. Dictionary values are dictionaries
    used as options for the corresponding plot. The allowed options for the
    different plot types are given below.
savefig_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.savefig`. By
    default, uses ``bbox_inches: tight, dpi: 300, orientation: landscape``.

seaborn_settings: dict, optional
    Options for :func:`seaborn.set` (affects all plots).

Configuration options for plot type ``timeseries``
--------------------------------------------------
annual_mean_kwargs: dict, optional
    Optional keyword arguments for :func:`iris.plot.plot` for plotting annual
    means. Use ``annual_mean_kwargs`` to not show annual means.
legend_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.legend`. Use
    ``legend_kwargs: false`` to not show legends.
plot_kwargs: dict, optional
    Optional keyword arguments for :func:`iris.plot.plot`. Dictionary keys are
    elements identified by ``facet_used_for_labels`` or ``default``, e.g.,
    ``CMIP6`` if ``facet_used_for_labels: project`` or ``historical`` if
    ``facet_used_for_labels: exp``. Dictionary values are dictionaries used as
    keyword arguments for :func:`iris.plot.plot`. String arguments can include
    facets in curly brackets which will be derived from the corresponding
    dataset, e.g., ``{project}``, ``{short_name}``, ``{exp}``. Examples:
    ``default: {linestyle: '-', label: '{project}'}, CMIP6: {color: red,
    linestyle: '--'}, OBS: {color: black}``.
pyplot_kwargs: dict, optional
    Optional calls to functions of :mod:`matplotlib.pyplot`. Dictionary keys
    are functions of :mod:`matplotlib.pyplot`. Dictionary values are used as
    single argument for these functions. String arguments can include facets in
    curly brackets which will be derived from the datasets plotted in the
    corresponding plot, e.g., ``{short_name}``, ``{exp}``. Facets like
    ``{project}`` that vary between the different datasets will be transformed
    to something like  ``ambiguous_project``. Examples: ``title: 'Awesome Plot
    of {long_name}', ``xlabel: '{short_name}'``, ``xlim: [0, 5]``.

Configuration options for plot type ``map``
-------------------------------------------
cbar_fontsize: int, optional (default: 10)
    Fontsize used for the colorbar label and ticks.
cbar_label: str, optional (default: '{short_name} [{units}]')
    Colorbar label. Can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``.
cbar_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.colorbar`. By
    default, uses ``orientation: horizontal, aspect: 40``.
gridline_kwargs: dict, optional
    Optional keyword arguments for grid lines. By default, ``color: lightgrey,
    alpha: 0.5`` are used. Use ``gridline_kwargs: false`` to not show grid
    lines.
plot_func: str, optional (default: 'contourf')
    Plot function used to plot the maps. Must be a function of :mod:`iris.plot`
    that supports plotting of 2D cubes with coordinates latitude and longitude.
plot_kwargs: dict, optional
    Optional keyword arguments for the plot function defined by ``plot_func``.
    Dictionary keys are elements identified by ``facet_used_for_labels`` or
    ``default``, e.g., ``CMIP6`` if ``facet_used_for_labels: project`` or
    ``historical`` if ``facet_used_for_labels: exp``. Dictionary values are
    dictionaries used as keyword arguments for :func:`iris.plot.plot`. String
    arguments can include facets in curly brackets which will be derived from
    the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``. Examples: ``default: {linestyle: '-', label: '{project}'},
    CMIP6: {color: red, linestyle: '--'}, OBS: {color: black}``.
projection: str, optional (default: 'Robinson')
    Projection used for the map plot. Needs to be a valid projection class of
    :mod:`cartopy.crs`. Keyword arguments can be specified using the option
    ``projection_kwargs``.
projection_kwargs: dict, optional
    Optional keyword arguments for the projection given by ``projection``. For
    the default projection ``Robinson``, the default keyword arguments
    ``central_longitude: 10`` are used.
pyplot_kwargs: dict, optional
    Optional calls to functions of :mod:`matplotlib.pyplot`. Dictionary keys
    are functions of :mod:`matplotlib.pyplot`. Dictionary values are used as
    single argument for these functions. String arguments can include facets in
    curly brackets which will be derived from the corresponding dataset, e.g.,
    ``{project}``, ``{short_name}``, ``{exp}``.  Examples: ``title: 'Awesome
    Plot of {long_name}', ``xlabel: '{short_name}'``, ``xlim: [0, 5]``.

"""
import logging
from copy import deepcopy
from pathlib import Path
from pprint import pformat

import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
import seaborn as sns
from iris.coord_categorisation import add_year
from matplotlib.gridspec import GridSpec

import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts.monitor.monitor_base import MonitorBase
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(Path(__file__).stem)


class MultiDatasets(MonitorBase):
    """Diagnostic to plot multi-dataset plots."""

    def __init__(self, config):
        """Initialize class member."""
        super().__init__(config)

        # Get default settings
        self.cfg = deepcopy(self.cfg)
        self.cfg.setdefault('facet_used_for_labels', 'dataset')
        self.cfg.setdefault('figure_kwargs', {})
        self.cfg.setdefault('savefig_kwargs', {
            'bbox_inches': 'tight',
            'dpi': 300,
            'orientation': 'landscape',
        })
        self.cfg.setdefault('seaborn_settings', {})
        logger.info("Using facet '%s' to create labels",
                    self.cfg['facet_used_for_labels'])

        # Load input data
        self.input_data = self._load_and_preprocess_data()
        self.grouped_input_data = group_metadata(self.input_data, 'short_name')

        # Check given plot types
        self.supported_plot_types = ['timeseries', 'map', 'profile']
        for (plot_type, plot_options) in self.plots.items():
            if plot_type not in self.supported_plot_types:
                raise ValueError(
                    f"Got unexpected plot type '{plot_type}' for option "
                    f"'plots', expected one of {self.supported_plot_types}")
            if plot_options is None:
                self.plots[plot_type] = {}

        # Check that facet_used_for_labels is present for every dataset
        for dataset in self.input_data:
            if self.cfg['facet_used_for_labels'] not in dataset:
                raise ValueError(
                    f"facet_used_for_labels "
                    f"'{self.cfg['facet_used_for_labels']}' not present for "
                    f"the following dataset:\n{pformat(dataset)}")

        # Load seaborn settings
        sns.set(**self.cfg['seaborn_settings'])

    def _get_label(self, dataset):
        """Get label of dataset."""
        return dataset[self.cfg['facet_used_for_labels']]

    def _get_cbar_kwargs(self, plot_type):
        """Get colorbar kwargs."""
        cbar_kwargs = {'orientation': 'horizontal', 'aspect': 40}
        cbar_kwargs.update(
            self.plots[plot_type].get('cbar_kwargs', {}))
        return cbar_kwargs

    def _get_cbar_label(self, plot_type, dataset):
        """Get colorbar label."""
        cbar_label = self.plots[plot_type].get(
            'cbar_label',
            '{short_name} [{units}]',
        )
        cbar_label = self._fill_facet_placeholders(
            cbar_label,
            dataset,
            f"cbar_label of {plot_type} '{cbar_label}'",
        )
        return cbar_label

    def _get_gridline_kwargs(self):
        """Get gridline kwargs."""
        plot_type = 'map'
        gridline_kwargs = self.plots[plot_type].get('gridline_kwargs', {})
        return gridline_kwargs

    def _get_map_projection(self):
        """Get projection used for map plots."""
        plot_type = 'map'

        # If no projection is specified, use Robinson with a set of default
        # kwargs
        if 'projection' not in self.plots[plot_type]:
            projection = 'Robinson'
            projection_kwargs = {'central_longitude': 10}
        else:
            projection = self.plots[plot_type]['projection']
            projection_kwargs = {}
        projection_kwargs.update(
            self.plots[plot_type].get('projection_kwargs', {})
        )

        # Check if desired projection is valid
        if not hasattr(ccrs, projection):
            raise AttributeError(
                f"Got invalid projection '{projection}' for plotting "
                f"{plot_type}, expected class of cartopy.crs")

        return getattr(ccrs, projection)(**projection_kwargs)

    def _get_plot_kwargs(self, plot_type, dataset):
        """Get keyword arguments for plot functions."""
        all_plot_kwargs = self.plots[plot_type].get('plot_kwargs', {})

        # First get default kwargs, then overwrite them with dataset-specific
        # ones
        plot_kwargs = all_plot_kwargs.get('default', {})
        label = self._get_label(dataset)
        plot_kwargs.update(all_plot_kwargs.get(label, {}))

        # Replace facets with dataset entries for string arguments
        for (key, val) in plot_kwargs.items():
            if isinstance(val, str):
                val = self._fill_facet_placeholders(
                    val,
                    dataset,
                    f"plot_kwargs of {plot_type} '{key}: {val}'",
                )
                plot_kwargs[key] = val

        # Set default label
        if plot_type == 'timeseries':
            plot_kwargs.setdefault('label', label)

        return plot_kwargs

    def _load_and_preprocess_data(self):
        """Load and preprocess data."""
        input_data = list(self.cfg['input_data'].values())
        for dataset in input_data:
            filename = dataset['filename']
            logger.info("Loading %s", filename)
            cube = iris.load_cube(filename)
            if cube.coords('time', dim_coords=True):
                ih.unify_time_coord(cube)
            dataset['cube'] = cube
        return input_data

    def _plot_map_with_ref(self, plot_func, dataset, ref_dataset):
        """Plot map plot for single dataset with a reference dataset."""
        plot_type = 'map'
        logger.debug("Plotting map with reference dataset '%s' for '%s'",
                     self._get_label(ref_dataset), self._get_label(dataset))

        # Make sure that the data has the correct dimensions
        cube = dataset['cube']
        ref_cube = ref_dataset['cube']
        self._check_cube_dimensions(cube, ['latitude', 'longitude'], plot_type)
        self._check_cube_dimensions(ref_cube, ['latitude', 'longitude'],
                                    plot_type)

        # Create single figure with multiple axes
        fig = plt.figure(**self.cfg['figure_kwargs'])
        gridspec = GridSpec(4, 4, figure=fig)

        # Options used for all subplots
        projection = self._get_map_projection()
        plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
        gridline_kwargs = self._get_gridline_kwargs()
        cbar_fontsize = self.plots[plot_type].get('cbar_fontsize', 10)
        cbar_kwargs = self._get_cbar_kwargs(plot_type)
        cbar_label = self._get_cbar_label(plot_type, dataset)

        # Plot dataset (top left)
        axes_data = fig.add_subplot(gridspec[0:2, 0:2], projection=projection)
        plot_kwargs['axes'] = axes_data
        plot_data = plot_func(cube, **plot_kwargs)
        axes_data.coastlines()
        if gridline_kwargs is not False:
            axes_data.gridlines(**gridline_kwargs)
        axes_data.set_title(
            self._get_label(dataset), fontsize=cbar_fontsize + 2)
        cbar_data = plt.colorbar(plot_data, ax=axes_data, **cbar_kwargs)
        cbar_data.set_label(cbar_label, labelpad=0.0, fontsize=cbar_fontsize)
        cbar_data.ax.tick_params(labelsize=cbar_fontsize)

        # Plot reference dataset (top right)
        axes_ref = fig.add_subplot(gridspec[2:4, 0:2], projection=projection)
        plot_kwargs['axes'] = axes_ref
        plot_ref = plot_func(ref_cube, **plot_kwargs)
        axes_ref.coastlines()
        if gridline_kwargs is not False:
            axes_ref.gridlines(**gridline_kwargs)
        axes_ref.set_title(
            self._get_label(ref_dataset), fontsize=cbar_fontsize + 2)
        cbar_ref = plt.colorbar(plot_ref, ax=axes_ref, **cbar_kwargs)
        cbar_ref.set_label(cbar_label, labelpad=0.0, fontsize=cbar_fontsize)
        cbar_ref.ax.tick_params(labelsize=cbar_fontsize)

        # Plot bias
        # bias_cube = cube - ref_cube
        # plot_kwargs = dict(norm=colors.CenteredNorm(), cmap='bwr')
        # plot_kwargs.update(cfg.get('plot_kwargs_bias', {}))
        # map_plot = iris.plot.contourf(bias_cube, axes=axes_list[1][0],
        #                               **plot_kwargs)
        # colorbar = plt.colorbar(map_plot, ax=axes_list[1][0],
        #                         orientation='horizontal', aspect=30)
        # colorbar.set_label(rf"$\Delta${cube.var_name} [{cube.units}]",
        #                    labelpad=0.0, fontsize=fontsize)
        # colorbar.ax.tick_params(labelsize=fontsize)

        # Plot appearance
        # title_kwargs = dict(pad=3.0, fontdict=dict(fontsize=10))
        # title_key = cfg['title_key']
        # axes_list[0][0].set_title(dataset[title_key], **title_kwargs)
        # axes_list[0][1].set_title(ref_dataset[title_key], **title_kwargs)
        # axes_list[1][0].set_title(
        #     f"{dataset[cfg['title_key']]} - {ref_dataset[cfg['title_key']]}",
        #     **title_kwargs)
        # for axes in list(axes_list.ravel()):
        #     axes.gridlines(color='lightgrey', alpha=0.5)
        #     axes.coastlines()
        #     axes.set_global()
        # axes_list[1][1].set_visible(False)

        # # Return figure and basename for plot
        # tag = dataset[cfg['title_key']].replace(' ', '_')
        # tag_ref = ref_dataset[cfg['title_key']].replace(' ', '_')
        # return (fig, f"{tag}_vs_{tag_ref}")

        # Customize plot
        plt.suptitle(dataset['long_name'])
        self._process_pyplot_kwargs(plot_type, dataset)

        return self.get_plot_path(plot_type, dataset)

    def _plot_map_without_ref(self, plot_func, dataset):
        """Plot map plot for single dataset without a reference dataset."""
        plot_type = 'map'
        logger.debug("Plotting map without reference dataset for '%s'",
                     self._get_label(dataset))

        # Make sure that the data has the correct dimensions
        cube = dataset['cube']
        self._check_cube_dimensions(cube, ['latitude', 'longitude'], plot_type)

        # Create plot with desired settings
        fig = plt.figure(**self.cfg['figure_kwargs'])
        axes = fig.add_subplot(projection=self._get_map_projection())
        plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
        plot_kwargs['axes'] = axes
        plot_func(cube, **plot_kwargs)
        axes.coastlines()
        gridline_kwargs = self._get_gridline_kwargs()
        if gridline_kwargs is not False:
            axes.gridlines(**gridline_kwargs)

        # Setup colorbar
        cbar_fontsize = self.plots[plot_type].get('cbar_fontsize', 10)
        colorbar = plt.colorbar(**self._get_cbar_kwargs(plot_type))
        colorbar.set_label(self._get_cbar_label(plot_type, dataset),
                           fontsize=cbar_fontsize)
        colorbar.ax.tick_params(labelsize=cbar_fontsize)

        # Customize plot
        plt.title(self._get_label(dataset))
        self._process_pyplot_kwargs(plot_type, dataset)

        return self.get_plot_path(plot_type, dataset)

    def _process_pyplot_kwargs(self, plot_type, dataset):
        """Process functions for :mod:`matplotlib.pyplot`."""
        pyplot_kwargs = self.plots[plot_type].get('pyplot_kwargs', {})
        for (func, arg) in pyplot_kwargs.items():
            if isinstance(arg, str):
                arg = self._fill_facet_placeholders(
                    arg,
                    dataset,
                    f"pyplot_kwargs of {plot_type} '{func}: {arg}'",
                )
            getattr(plt, func)(arg)

    @staticmethod
    def _check_cube_dimensions(cube, expected_dimensions, plot_type):
        """Check that cube has correct dimensional variables."""
        ndim = len(expected_dimensions)
        if cube.ndim != ndim:
            raise ValueError(
                f"Expected {ndim:d}D cube with dimensions "
                f"{expected_dimensions} for plot {plot_type}, got "
                f"{cube.ndim:d}D cube: {cube.summary(shorten=True)}")

    @staticmethod
    def _fill_facet_placeholders(string, dataset, description):
        """Fill facet placeholders."""
        try:
            string = string.format(**dataset)
        except KeyError as exc:
            raise ValueError(
                f"Not all necessary facets in {description} available for "
                f"dataset\n{pformat(dataset)}") from exc
        return string

    @staticmethod
    def _get_multi_dataset_facets(datasets):
        """Derive common facets for multiple datasets."""
        all_keys = {key for dataset in datasets for key in dataset}
        multi_dataset_facets = {}
        for key in all_keys:
            if all(d.get(key) == datasets[0].get(key) for d in datasets):
                multi_dataset_facets[key] = datasets[0][key]
            else:
                multi_dataset_facets[key] = f'ambiguous_{key}'
        return multi_dataset_facets

    @staticmethod
    def _get_reference_dataset(datasets, short_name):
        """Extract reference dataset."""
        ref_datasets = [d for d in datasets if
                        d.get('reference_for_monitor_diags', False)]
        if len(ref_datasets) > 1:
            raise ValueError(
                f"Expected at most 1 reference dataset (with "
                f"'reference_for_monitor_diags: true' for variable "
                f"'{short_name}', got {len(ref_datasets):d}")
        if ref_datasets:
            return ref_datasets[0]
        return None

    def create_timeseries_plot(self, datasets, short_name):
        """Create time series plot."""
        plot_type = 'timeseries'
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        logger.info("Plotting %s", plot_type)
        fig = plt.figure(**self.cfg['figure_kwargs'])
        axes = fig.add_subplot()

        # Plot all datasets in one single figure
        ancestors = []
        for dataset in datasets:
            ancestors.append(dataset['filename'])
            cube = dataset['cube']
            self._check_cube_dimensions(cube, ['time'], plot_type)

            # Plot original time series
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs['axes'] = axes
            iris.plot.plot(cube, **plot_kwargs)

            # Plot annual means if desired
            annual_mean_kwargs = self.plots[
                plot_type].get('annual_mean_kwargs', {})
            if annual_mean_kwargs is not False:
                logger.debug("Plotting annual means")
                if not cube.coords('year'):
                    add_year(cube, 'time')
                annual_mean_cube = cube.aggregated_by('year',
                                                      iris.analysis.MEAN)
                plot_kwargs.pop('label', None)
                plot_kwargs.update(annual_mean_kwargs)
                iris.plot.plot(annual_mean_cube, **plot_kwargs)

        # Default plot appearance
        multi_dataset_facets = self._get_multi_dataset_facets(datasets)
        plt.title(multi_dataset_facets['long_name'])
        plt.xlabel('Time')
        plt.ylabel(f"{short_name} [{multi_dataset_facets['units']}]")

        # Legend
        legend_kwargs = self.plots[plot_type].get('legend_kwargs', {})
        if legend_kwargs is not False:
            plt.legend(**legend_kwargs)

        # Customize plot appearance
        self._process_pyplot_kwargs(plot_type, multi_dataset_facets)

        # Save plot
        plot_path = self.get_plot_path(plot_type, multi_dataset_facets)
        plt.savefig(plot_path, **self.cfg['savefig_kwargs'])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Provenance tracking
        caption = (f"Time series of {multi_dataset_facets['long_name']} for "
                   f"various datasets.")
        provenance_record = {
            'ancestors': ancestors,
            'authors': ['schlund_manuel'],
            'caption': caption,
            'plot_types': ['line'],
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)

    def create_map_plot(self, datasets, short_name):
        """Create map plot."""
        plot_type = 'map'
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        # Get reference dataset if possible
        ref_dataset = self._get_reference_dataset(datasets, short_name)
        if ref_dataset is None:
            logger.info("Plotting %s without reference dataset", plot_type)
        else:
            logger.info("Plotting %s with reference dataset '%s'", plot_type,
                        self._get_label(ref_dataset))

        # Get plot function
        plot_func = self.plots.get('plot_func', 'contourf')
        if not hasattr(iris.plot, plot_func):
            raise AttributeError(
                f"Got invalid plot function '{plot_func}' for plotting "
                f"{plot_type}, expected function of iris.plot")
        logger.info("Creating %s plots using function '%s'", plot_type,
                    plot_func)
        plot_func = getattr(iris.plot, plot_func)

        # Create a single plot for each dataset (incl. reference dataset if
        # given)
        for dataset in datasets:
            if ref_dataset == dataset:
                continue
            ancestors = [dataset['filename']]
            if ref_dataset is None:
                plot_path = self._plot_map_without_ref(plot_func, dataset)
                caption = (
                    f"Map plot of {dataset['long_name']} of dataset "
                    f"{dataset['dataset']} (project {dataset['project']}) "
                    f"from {dataset['start_year']} to {dataset['end_year']}."
                )
            else:
                plot_path = self._plot_map_with_ref(plot_func, dataset,
                                                    ref_dataset)
                ancestors.append(ref_dataset['filename'])
                caption = (
                    f"Map plot of {dataset['long_name']} of dataset "
                    f"{dataset['dataset']} (project {dataset['project']}) "
                    f"including bias relative to {ref_dataset['dataset']} "
                    f"(project {ref_dataset['project']}) from "
                    f"{dataset['start_year']} to {dataset['end_year']}."
                )

        # Save plot
        plt.savefig(plot_path, **self.cfg['savefig_kwargs'])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Provenance tracking
        provenance_record = {
            'ancestors': ancestors,
            'authors': ['schlund_manuel'],
            'caption': caption,
            'plot_types': ['map'],
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)

    def compute(self):
        """Plot preprocessed data."""
        for (short_name, datasets) in self.grouped_input_data.items():
            logger.info("Processing variable %s", short_name)
            self.create_timeseries_plot(datasets, short_name)
            self.create_map_plot(datasets, short_name)


def main():
    """Run diagnostic."""
    with run_diagnostic() as config:
        MultiDatasets(config).compute()


if __name__ == '__main__':
    main()
