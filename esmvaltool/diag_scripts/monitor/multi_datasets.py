#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Monitoring diagnostic to show multiple datasets in one plot (incl. biases).

Description
-----------
This diagnostic can be used to visualize multiple datasets in one plot.

... reference dataset ...


Currently supported plot types:
    - Time series (data needs to be 1D with single dimension 'time'): for each
      variable separately, all datasets are plotted in one single figure.
    - Maps (data needs to be 2D with dimensions 'latitude', 'longitude'): for
      each variable and dataset, an individual map is plotted. If a reference
      dataset is defined, also include this and a bias plot into the figure.
    - Profiles (data needs to be 2D with dimensions 'latitude',
      'height'/'air_pressure'): for each variable and dataset, an individual
      profile is plotted. If a reference dataset is defined, also include this
      and a bias plot into the figure.

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
    'timeseries', 'map', or 'profile'. Dictionary values are dictionaries used
    as options for the corresponding plot. The allowed options for the
    different plot types are given below.
savefig_kwargs: dict, optional
    Keyword arguments for :func:`matplotlib.pyplot.savefig`.
seaborn_settings: dict, optional
    Options for :func:`seaborn.set` (affects all plots).

Configuration options for timeseries
------------------------------------
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
    to something like  ``ambiguous_project``.  Examples: ``title: 'Awesome Plot
    of {long_name}', ``xlabel: '{short_name}'``, ``xlim: [0, 5]``.

"""
import logging
from copy import deepcopy
from pathlib import Path
from pprint import pformat

import iris
import matplotlib.pyplot as plt
import seaborn as sns
from iris.coord_categorisation import add_year

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
                try:
                    val = val.format(**dataset)
                except KeyError as exc:
                    raise ValueError(
                        f"Not all necessary facets in {plot_type} plot_kwargs "
                        f"'{key}: {val}' available for dataset"
                        f"\n{pformat(dataset)}") from exc
                plot_kwargs[key] = val

        # Set default label
        plot_kwargs.setdefault('label', label)

        return plot_kwargs

    def _load_and_preprocess_data(self):
        """Load and preprocess data."""
        input_data = list(self.cfg['input_data'].values())
        for dataset in input_data:
            filename = dataset['filename']
            logger.info("Loading %s", filename)
            cube = iris.load_cube(filename)
            ih.unify_time_coord(cube)
            dataset['cube'] = cube
        return input_data

    def _process_pyplot_kwargs(self, plot_type, dataset):
        """Process functions for :mod:`matplotlib.pyplot`."""
        pyplot_kwargs = self.plots[plot_type].get('pyplot_kwargs', {})
        for (func, arg) in pyplot_kwargs.items():
            if isinstance(arg, str):
                try:
                    arg = arg.format(**dataset)
                except KeyError as exc:
                    raise ValueError(
                        f"Not all necessary facets in {plot_type} "
                        f"pyplot_kwargs '{func}: {arg}' available for dataset"
                        f"\n{pformat(dataset)}") from exc
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

    def create_timeseries_plot(self, datasets):
        """Create time series plot."""
        plot_type = 'timeseries'
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError("No input data given")

        logger.info("Plotting %s", plot_type)
        plt.figure(**self.cfg['figure_kwargs'])

        # Plot all datasets in single figure
        ancestors = []
        for dataset in datasets:
            ancestors.append(dataset['filename'])
            cube = dataset['cube']
            self._check_cube_dimensions(cube, ['time'], plot_type)

            # Plot original time series
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
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
        plt.ylabel(
            f"{multi_dataset_facets['short_name']} "
            f"[{multi_dataset_facets['units']}]")

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

    def compute(self):
        """Plot preprocessed data."""
        for (short_name, datasets) in self.grouped_input_data.items():
            logger.info("Processing variable %s", short_name)
            self.create_timeseries_plot(datasets)


def main():
    """Run diagnostic."""
    with run_diagnostic() as config:
        MultiDatasets(config).compute()


if __name__ == '__main__':
    main()
