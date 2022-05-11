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
    Optional keyword arguments of :func:`matplotlib.pyplot.figure`.
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
legend_kwargs: dict, optional
    Optional keyword arguments of :func:`matplotlib.pyplot.legend`. Use
    ``legend_kwargs: false`` to not show legends.
plot_kwargs: dict, optional
    Additional keyword arguments for :func:`iris.plot.plot`. Dictionary keys
    are elements identified by ``facet_used_for_labels`` or ``default``, e.g.,
    ``CMIP6`` if ``facet_used_for_labels: project`` or ``historical`` if
    ``facet_used_for_labels: exp``. Dictionary values are dictionaries used as
    keyword arguments for :func:`iris.plot.plot`. Examples: ``default:
    linestyle: '-', CMIP6: {color: red, linestyle: '--'}, OBS: {color:
    black}``.
pyplot_kwargs: dict, optional
    Optional calls to functions of :mod:`matplotlib.pyplot`. Dictionary keys
    are functions of :mod:`matplotlib.pyplot`. Dictionary values are used as
    single  argument for these functions. Examples: ``title: Awesome Plot,
    ``xlabel: x``, ``xlim: [0, 5]``.

"""
import logging
from copy import deepcopy
from pathlib import Path
from pprint import pformat

import iris
import matplotlib.pyplot as plt
import seaborn as sns

import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts.monitor.monitor_base import MonitorBase
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_plot_filename,
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
        self.cfg.setdefault('legend_kwargs', {})
        self.cfg.setdefault('savefig_kwargs', {
            'bbox_inches': 'tight',
            'dpi': 300,
            'orientation': 'landscape',
        })
        self.cfg.setdefault('seaborn_settings', {})
        logger.info("Using key '%s' to create titles for datasets",
                    self.cfg['facet_used_for_labels'])

        # Load input data
        self.input_data = self._load_and_preprocess_data()
        self.grouped_input_data = group_metadata(self.input_data, 'short_name')

        # Check given plot types
        self.supported_plot_types = ['timeseries', 'map', 'profile']
        for plot_type in self.plots:
            if plot_type not in self.supported_plot_types:
                raise ValueError(
                    f"Got unexpected plot type '{plot_type}' for option "
                    f"'plots', expected one of {self.supported_plot_types}")

        # Check that facet_used_for_labels is present for every dataset
        for dataset in self.input_data:
            if self.cfg['facet_used_for_labels'] not in dataset:
                raise ValueError(
                    f"facet_used_for_labels "
                    f"'{self.cfg['facet_used_for_labels']}' not present for "
                    f"the following dataset:\n{pformat(dataset)}")

        # Load seaborn settings
        sns.set(**self.cfg['seaborn_settings'])

    def _get_plot_kwargs(self, plot_type, dataset):
        """Get keyword arguments for plot functions."""
        all_plot_kwargs = self.plots[plot_type].get('plot_kwargs', {})
        plot_kwargs = all_plot_kwargs.get('default', {})
        label = dataset[self.cfg['facet_used_for_labels']]
        plot_kwargs.update(all_plot_kwargs.get(label, {}))
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

    def _process_pyplot_kwargs(self, plot_type):
        """Process functions for :mod:`matplotlib.pyplot`."""
        pyplot_kwargs = self.plots[plot_type].get('pyplot_kwargs', {})
        for (func, arg) in pyplot_kwargs.items():
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

    def create_timeseries_plot(self, datasets):
        """Create time series plot."""
        plot_type = 'timeseries'
        if plot_type not in self.plots:
            return
        plt.figure(**self.cfg['figure_kwargs'])

        # Plot all datasets in single figure
        ancestors = []
        for dataset in datasets:
            ancestors.append(dataset['filename'])
            cube = dataset['cube']
            self._check_cube_dimensions(cube, ['time'], plot_type)

            # Plot monthly means
        #     plot_kwargs = get_plot_kwargs(cfg, dataset)
        #     iris.plot.plot(cube, label=dataset[cfg['facet_used_for_labels']], **plot_kwargs)

        #     # Plot annual means
        #     iris.coord_categorisation.add_year(cube, 'time')
        #     cube = cube.aggregated_by('year', iris.analysis.MEAN)
        #     plot_kwargs['linestyle'] = ':'
        #     iris.plot.plot(cube, **plot_kwargs)

        # # Plot appearance
        # long_name = input_data[0]['long_name']
        # short_name = input_data[0]['short_name']
        # units = cube.units
        # plt.title(f"Global Mean {long_name}")
        # plt.xlabel("Year")
        # plt.ylabel(f"{short_name} [{units}]")
        # plt.legend()

        # # Save plot
        # plot_path = get_plot_filename(short_name, cfg)
        # plt.savefig(plot_path, bbox_inches='tight', orientation='landscape')
        # logger.info("Wrote %s", plot_path)
        # plt.close()

        # # Provenance tracking
        # caption = (f"Monthly mean (solid lines) and annual mean (dashed lines) "
        #            f"time series of {input_data[0]['long_name']} for various "
        #            f"datasets.")
        # provenance_record = {
        #     'ancestors': ancestors,
        #     'authors': ['schlund_manuel'],
        #     'caption': caption,
        #     'plot_types': ['line'],
        # }
        # with ProvenanceLogger(cfg) as provenance_logger:
        #     provenance_logger.log(plot_path, provenance_record)

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
