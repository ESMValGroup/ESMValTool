# -*- coding: utf-8 -*-
"""Lifetime diagnostic to show multiple datasets in one plot.

Description
-----------
This diagnostic can be used to visualize lifetime of multiple datasets in
one plot.

For some plot types, a reference dataset can be defined. For this, use the
facet ``reference_for_lifetime_diags: true`` in the definition of the dataset
in the recipe. Note that at most one reference dataset is supported.

Currently supported plot types (use the option ``plots`` to specify them):
    - Time series (plot type ``timeseries``): all datasets are plotted in one
      single figure.
    - Annual cycle (plot type ``annual_cycle``): all datasets are plotted in
      one single figure.
    - Zonal mean profiles (plot type ``zonalmean``):
      for each dataset, an individual profile is plotted. If a reference
      dataset    is defined, also include this dataset and a bias plot
      into the figure. Note that if a reference dataset is defined, all input
      datasets need to be given on the same horizontal and vertical grid (you
      can use the preprocessors :func:`esmvalcore.preprocessor.regrid` and
      :func:`esmvalcore.preprocessor.extract_levels` for this). Input data
      needs to be 2D with dimensions `latitude`, `height`/`air_pressure`.
    - 1D profiles (plot type ``1d_profile``): all datasets are plotted in one
      single figure.

Author
------
Franziska Winterstein (DLR, Germany)

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
    Optional keyword arguments for :func:`matplotlib.pyplot.figure`. By
    default, uses ``constrained_layout: true``.
plots: dict, optional
    Plot types plotted by this diagnostic (see list above). Dictionary keys
    must be ``timeseries``, ``annual_cycle``, ``map``, ``zonalmean``
    or ``1d_profile``.
    Dictionary values are dictionaries used as options for the corresponding
    plot. The allowed options for the different plot types are given below.
plot_filename: str, optional
    Filename pattern for the plots.
    Defaults to ``{plot_type}_{real_name}_{dataset}_{mip}_{exp}_{ensemble}``.
    All tags (i.e., the entries in curly brackets, e.g., ``{dataset}``, are
    replaced with the corresponding tags).
plot_folder: str, optional
    Path to the folder to store figures. Defaults to
    ``{plot_dir}/../../{dataset}/{exp}/{modeling_realm}/{real_name}``.  All
    tags (i.e., the entries in curly brackets, e.g., ``{dataset}``, are
    replaced with the corresponding tags).  ``{plot_dir}`` is replaced with the
    default ESMValTool plot directory (i.e.,
    ``output_dir/plots/diagnostic_name/script_name/``, see
    :ref:`esmvalcore:user configuration file`).
savefig_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.savefig`. By
    default, uses ``bbox_inches: tight, dpi: 300, orientation: landscape``.
seaborn_settings: dict, optional
    Options for :func:`seaborn.set_theme` (affects all plots). By default, uses
    ``style: ticks``.

Configuration options for plot type ``timeseries``
--------------------------------------------------
annual_mean: str, optional
    Optional switch to turn on annual means to be displayed  'only' or
    additional 'both' to the original timeseries. If not set or set to 'False'
    only the original timeseries is shown.
annual_mean_kwargs: dict, optional
    Optional keyword arguments for :func:`iris.plot.plot` for plotting annual
    means. These keyword arguments update (and potentially overwrite) the
    ``plot_kwargs`` for the annual mean plots. Use ``annual_mean_kwargs`` to
    not show annual means.
gridline_kwargs: dict, optional
    Optional keyword arguments for grid lines. By default, ``color: lightgrey,
    alpha: 0.5`` are used. Use ``gridline_kwargs: false`` to not show grid
    lines.
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
    of {long_name}'``, ``xlabel: '{short_name}'``, ``xlim: [0, 5]``.

Configuration options for plot type ``annual_cycle``
----------------------------------------------------
gridline_kwargs: dict, optional
    Optional keyword arguments for grid lines. By default, ``color: lightgrey,
    alpha: 0.5`` are used. Use ``gridline_kwargs: false`` to not show grid
    lines.
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
    of {long_name}'``, ``xlabel: '{short_name}'``, ``xlim: [0, 5]``.


Configuration options for plot type ``zonalmean``
----------------------------------------------------------
cbar_label: str, optional (default: '{short_name} [{units}]')
    Colorbar label. Can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``.
cbar_label_bias: str, optional (default: 'Δ{short_name} [{units}]')
    Colorbar label for plotting biases. Can include facets in curly brackets
    which will be derived from the corresponding dataset, e.g., ``{project}``,
    ``{short_name}``, ``{exp}``. This option has no effect if no reference
    dataset is given.
cbar_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.colorbar`. By
    default, uses ``orientation: vertical``.
cbar_kwargs_bias: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.colorbar` for
    plotting biases. These keyword arguments update (and potentially overwrite)
    the ``cbar_kwargs`` for the bias plot. This option has no effect if no
    reference dataset is given.
common_cbar: bool, optional (default: False)
    Use a common colorbar for the top panels (i.e., plots of the dataset and
    the corresponding reference dataset) when using a reference dataset. If
    neither ``vmin`` and ``vmix`` nor ``levels`` is given in ``plot_kwargs``,
    the colorbar bounds are inferred from the dataset in the top left panel,
    which might lead to an inappropriate colorbar for the reference dataset
    (top right panel). Thus, the use of the ``plot_kwargs`` ``vmin`` and
    ``vmax`` or ``levels`` is highly recommend when using this ``common_cbar:
    true``. This option has no effect if no reference dataset is given.
fontsize: int, optional (default: 10)
    Fontsize used for ticks, labels and titles. For the latter, use the given
    fontsize plus 2. Does not affect suptitles.
log_y: bool, optional (default: True)
    Use logarithmic Y-axis.
plot_func: str, optional (default: 'contourf')
    Plot function used to plot the profiles. Must be a function of
    :mod:`iris.plot` that supports plotting of 2D cubes with coordinates
    latitude and height/air_pressure.
plot_kwargs: dict, optional
    Optional keyword arguments for the plot function defined by ``plot_func``.
    Dictionary keys are elements identified by ``facet_used_for_labels`` or
    ``default``, e.g., ``CMIP6`` if ``facet_used_for_labels: project`` or
    ``historical`` if ``facet_used_for_labels: exp``. Dictionary values are
    dictionaries used as keyword arguments for the plot function defined by
    ``plot_func``. String arguments can include facets in curly brackets which
    will be derived from the corresponding dataset, e.g., ``{project}``,
    ``{short_name}``, ``{exp}``. Examples: ``default: {levels: 2}, CMIP6:
    {vmin: 200, vmax: 250}``.
plot_kwargs_bias: dict, optional
    Optional keyword arguments for the plot function defined by ``plot_func``
    for plotting biases. These keyword arguments update (and potentially
    overwrite) the ``plot_kwargs`` for the bias plot. This option has no effect
    if no reference dataset is given. See option ``plot_kwargs`` for more
    details. By default, uses ``cmap: bwr``.
pyplot_kwargs: dict, optional
    Optional calls to functions of :mod:`matplotlib.pyplot`. Dictionary keys
    are functions of :mod:`matplotlib.pyplot`. Dictionary values are used as
    single argument for these functions. String arguments can include facets in
    curly brackets which will be derived from the corresponding dataset, e.g.,
    ``{project}``, ``{short_name}``, ``{exp}``.  Examples: ``title: 'Awesome
    Plot of {long_name}'``, ``xlabel: '{short_name}'``, ``xlim: [0, 5]``.
rasterize: bool, optional (default: True)
    If ``True``, use `rasterization
    <https://matplotlib.org/stable/gallery/misc/rasterization_demo.html>`_ for
    profile plots to produce smaller files. This is only relevant for vector
    graphics (e.g., ``output_file_type=pdf,svg,ps``).
show_stats: bool, optional (default: True)
    Show basic statistics on the plots.
show_y_minor_ticklabels: bool, optional (default: False)
    Show tick labels for the minor ticks on the Y axis.
x_pos_stats_avg: float, optional (default: 0.01)
    Text x-position of average (shown on the left) in Axes coordinates. Can be
    adjusted to avoid overlap with the figure. Only relevant if ``show_stats:
    true``.
x_pos_stats_bias: float, optional (default: 0.7)
    Text x-position of bias statistics (shown on the right) in Axes
    coordinates. Can be adjusted to avoid overlap with the figure. Only
    relevant if ``show_stats: true``.

Configuration options for plot type ``1d_profile``
--------------------------------------------------
aspect_ratio: float, optional (default: 1.5)
    Aspect ratio of the plot. The default value results in a slender upright
    plot.
gridline_kwargs: dict, optional
    Optional keyword arguments for grid lines. By default, ``color: lightgrey,
    alpha: 0.5`` are used. Use ``gridline_kwargs: false`` to not show grid
    lines.
legend_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.legend`. Use
    ``legend_kwargs: false`` to not show legends.
log_x: bool, optional (default: False)
    Use logarithmic X-axis. Note that for the logarithmic x axis tickmarks are
    set so that minor tickmarks show up. Setting of individual tickmarks by
    pyplot_kwargs is not recommended in this case.
log_y: bool, optional (default: True)
    Use logarithmic Y-axis.
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
    of {long_name}'``, ``xlabel: '{short_name}'``, ``xlim: [0, 5]``.
show_y_minor_ticklabels: bool, optional (default: False)
    Show tick labels for the minor ticks on the Y axis.

.. hint::

   Extra arguments given to the recipe are ignored, so it is safe to use yaml
   anchors to share the configuration of common arguments with other monitor
   diagnostic script.

"""
import logging
from copy import deepcopy
from pathlib import Path
from pprint import pformat

import iris
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from iris.coord_categorisation import add_year
import iris.common
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import FormatStrFormatter, LogLocator, NullFormatter

from esmvalcore.cmor.fixes import add_model_level
import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts.lifetime.lifetime_base import (
    LifetimeBase)
from esmvaltool.diag_scripts.lifetime.lifetime_func import (
    create_press,
    calculate_gridmassdry,
    calculate_lifetime,
    calculate_reaction_rate,
    calculate_rho,
    climatological_tropopause
)
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    io,
    run_diagnostic,
)


logger = logging.getLogger(Path(__file__).stem)


class CH4Lifetime(LifetimeBase):
    """Diagnostic to plot ch4 lifetime."""

    def __init__(self, config):
        """Initialize class member."""
        super().__init__(config)

        # Get default stettings
        self.cfg = deepcopy(self.cfg)

        self.cfg.setdefault('facet_used_for_labels', 'dataset')
        self.cfg.setdefault('figure_kwargs', {'constrained_layout': True})
        self.cfg.setdefault('savefig_kwargs', {
            'bbox_inches': 'tight',
            'dpi': 300,
            'orientation': 'landscape',
        })
        self.cfg.setdefault('seaborn_settings', {'style': 'ticks'})
        logger.info("Using facet '%s' to create labels",
                    self.cfg['facet_used_for_labels'])
        self.cfg.setdefault('regions', ['TROP'])

        # set default molarmasses
        self.cfg.setdefault('m_air', 28.970)  # [g_air/mol_air]

        # Load input data
        self.input_data_dataset = self._calculate_coefficients()

        # base info
        self.info = {'short_name':
                     f"tau_{self._get_name('reactant').upper()}"}
        oxidants = [ox.upper() for ox in self._get_name('oxidant')]
        self.info['long_name'] = ("Lifetime of"
                                  f" {self._get_name('reactant').upper()}"
                                  " with respect to"
                                  f" {', '.join(oxidants)}")
        self.units = self.cfg['units']
        self.info['units'] = self.cfg['units']

        # Check given plot types and set default settings for them
        self.supported_plot_types = [
            'timeseries',
            'annual_cycle',
            'zonalmean',
            '1d_profile'
        ]
        for (plot_type, plot_options) in self.plots.items():
            if plot_type not in self.supported_plot_types:
                raise ValueError(
                    f"Got unexpected plot type '{plot_type}' for option "
                    f"'plots', expected one of {self.supported_plot_types}")
            if plot_options is None:
                self.plots[plot_type] = {}

            # Default options for the different plot types
            if plot_type == 'timeseries':
                self.plots[plot_type].setdefault('annual_mean', False)
                self.plots[plot_type].setdefault('annual_mean_kwargs', {})
                self.plots[plot_type].setdefault('gridline_kwargs', {})
                self.plots[plot_type].setdefault('legend_kwargs', {})
                self.plots[plot_type].setdefault('plot_kwargs', {})
                self.plots[plot_type].setdefault('pyplot_kwargs', {})
                self.plots[plot_type].setdefault('by_timestep', False)

            if plot_type == 'annual_cycle':
                self.plots[plot_type].setdefault('gridline_kwargs', {})
                self.plots[plot_type].setdefault('legend_kwargs', {})
                self.plots[plot_type].setdefault('plot_kwargs', {})
                self.plots[plot_type].setdefault('pyplot_kwargs', {})

            if plot_type == 'zonalmean':
                self.plots[plot_type].setdefault(
                    'cbar_label', '{short_name} [{units}]')
                self.plots[plot_type].setdefault(
                    'cbar_label_bias', 'Δ{short_name} [{units}]')
                self.plots[plot_type].setdefault(
                    'cbar_kwargs', {'orientation': 'vertical'}
                )
                self.plots[plot_type].setdefault('cbar_kwargs_bias', {})
                self.plots[plot_type].setdefault('common_cbar', False)
                self.plots[plot_type].setdefault('fontsize', 10)
                self.plots[plot_type].setdefault('log_y', True)
                self.plots[plot_type].setdefault('plot_func', 'contourf')
                self.plots[plot_type].setdefault('plot_kwargs', {})
                self.plots[plot_type].setdefault('plot_kwargs_bias', {})
                self.plots[plot_type]['plot_kwargs_bias'].setdefault(
                    'cmap', 'bwr'
                )
                self.plots[plot_type].setdefault('pyplot_kwargs', {})
                self.plots[plot_type].setdefault('rasterize', True)
                self.plots[plot_type].setdefault('show_stats', True)
                self.plots[plot_type].setdefault(
                    'show_y_minor_ticklabels', False
                )
                self.plots[plot_type].setdefault('x_pos_stats_avg', 0.01)
                self.plots[plot_type].setdefault('x_pos_stats_bias', 0.7)

            if plot_type == '1d_profile':
                self.plots[plot_type].setdefault('aspect_ratio', 1.5)
                self.plots[plot_type].setdefault('gridline_kwargs', {})
                self.plots[plot_type].setdefault('legend_kwargs', {})
                self.plots[plot_type].setdefault('log_x', False)
                self.plots[plot_type].setdefault('log_y', True)
                self.plots[plot_type].setdefault('plot_kwargs', {})
                self.plots[plot_type].setdefault('pyplot_kwargs', {})
                self.plots[plot_type].setdefault(
                    'show_y_minor_ticklabels', False
                )

        # Check that facet_used_for_labels is present for every dataset
        for dataset in self.input_data:
            if self.cfg['facet_used_for_labels'] not in dataset:
                raise ValueError(
                    f"facet_used_for_labels "
                    f"'{self.cfg['facet_used_for_labels']}' not present for "
                    f"the following dataset:\n{pformat(dataset)}")

        # Load seaborn settings
        sns.set_theme(**self.cfg['seaborn_settings'])

    def _add_colorbar(self, plot_type, plot_left, plot_right, axes_left,
                      axes_right, dataset_left, dataset_right):
        """Add colorbar(s) for plots."""
        fontsize = self.plots[plot_type]['fontsize']
        cbar_kwargs = self._get_cbar_kwargs(plot_type)
        cbar_label_left = self._get_cbar_label(plot_type, dataset_left)
        cbar_label_right = self._get_cbar_label(plot_type, dataset_right)

        # Create one common colorbar for the top panels
        # Note: Increase aspect ratio for nicer looks
        if self.plots[plot_type]['common_cbar']:
            if 'aspect' in cbar_kwargs:
                cbar_kwargs['aspect'] += 20.0
            cbar = plt.colorbar(plot_left, ax=[axes_left, axes_right],
                                **cbar_kwargs)
            cbar.set_label(cbar_label_left, fontsize=fontsize)
            cbar.ax.tick_params(labelsize=fontsize)

        # Create two separate colorbars for the top panels
        else:
            cbar_left = plt.colorbar(plot_left, ax=axes_left, **cbar_kwargs)
            cbar_left.set_label(cbar_label_left, fontsize=fontsize)
            cbar_left.ax.tick_params(labelsize=fontsize)
            cbar_right = plt.colorbar(plot_right, ax=axes_right, **cbar_kwargs)
            cbar_right.set_label(cbar_label_right, fontsize=fontsize)
            cbar_right.ax.tick_params(labelsize=fontsize)

    def _get_custom_mpl_rc_params(self, plot_type):
        """Get custom matplotlib rcParams."""
        fontsize = self.plots[plot_type]['fontsize']
        custom_rc_params = {
            'axes.titlesize': fontsize + 2.0,
            'axes.labelsize': fontsize,
            'xtick.labelsize': fontsize,
            'ytick.labelsize': fontsize,
        }
        return custom_rc_params

    def _get_label(self, dataset):
        """Get label of dataset."""
        return dataset[self.cfg['facet_used_for_labels']]

    def _get_cbar_kwargs(self, plot_type, bias=False):
        """Get colorbar kwargs."""
        cbar_kwargs = deepcopy(self.plots[plot_type]['cbar_kwargs'])
        if bias:
            cbar_kwargs.update(self.plots[plot_type]['cbar_kwargs_bias'])
        return deepcopy(cbar_kwargs)

    def _get_cbar_label(self, plot_type, dataset, bias=False):
        """Get colorbar label."""
        if bias:
            cbar_label = self.plots[plot_type]['cbar_label_bias']
            descr = f"cbar_label_bias of {plot_type} '{cbar_label}'"
        else:
            cbar_label = self.plots[plot_type]['cbar_label']
            descr = f"cbar_label of {plot_type} '{cbar_label}'"
        cbar_label = self._fill_facet_placeholders(cbar_label, dataset, descr)
        return cbar_label

    def _get_gridline_kwargs(self, plot_type):
        """Get gridline kwargs."""
        gridline_kwargs = self.plots[plot_type]['gridline_kwargs']
        return deepcopy(gridline_kwargs)

    def _get_plot_func(self, plot_type):
        """Get plot function."""
        plot_func = self.plots[plot_type]['plot_func']
        if not hasattr(iris.plot, plot_func):
            raise AttributeError(
                f"Got invalid plot function '{plot_func}' for plotting "
                f"{plot_type}, expected function of iris.plot")
        logger.info("Creating %s plots using function '%s'", plot_type,
                    plot_func)
        return getattr(iris.plot, plot_func)

    def _get_plot_kwargs(self, plot_type, dataset, bias=False):
        """Get keyword arguments for plot functions."""
        all_plot_kwargs = self.plots[plot_type]['plot_kwargs']
        all_plot_kwargs = deepcopy(all_plot_kwargs)

        # First get default kwargs, then overwrite them with dataset-specific
        # ones
        plot_kwargs = all_plot_kwargs.get('default', {})
        label = self._get_label(dataset)
        plot_kwargs.update(all_plot_kwargs.get(label, {}))

        # For bias plots, overwrite the kwargs with bias-specific option
        if bias:
            bias_kwargs = self.plots[plot_type]['plot_kwargs_bias']
            plot_kwargs.update(bias_kwargs)

        # Replace facets with dataset entries for string arguments
        for (key, val) in plot_kwargs.items():
            if isinstance(val, str):
                val = self._fill_facet_placeholders(
                    val,
                    dataset,
                    f"plot_kwargs of {plot_type} '{key}: {val}'",
                )
                plot_kwargs[key] = val

        # Default settings for different plot types
        if plot_type in ('timeseries', 'annual_cycle', '1d_profile'):
            plot_kwargs.setdefault('label', label)

        return deepcopy(plot_kwargs)

    def _calculate_coefficients(self):
        """Load data and calculate coefficients."""
        self.input_data = list(self.cfg['input_data'].values())

        # loops over different variables ch4, oh, ta etc.
        input_data_dataset = {}
        for dataset in self._get_all_datasets(self.input_data):
            input_data_dataset[dataset] = self._get_dataset_data(dataset)[0]
            input_data_dataset[dataset]['dataset_data'] = (
                self._get_dataset_data(dataset))

            variables = {}
            for variable in input_data_dataset[dataset]['dataset_data']:
                filename = variable['filename']

                logger.info("Loading %s", filename)
                cube = iris.load_cube(filename, variable['short_name'])

                # Fix time coordinate if present
                if cube.coords('time', dim_coords=True):
                    ih.unify_time_coord(cube)

                # cube for each variable
                variables[variable['short_name']] = cube

            rho = calculate_rho(variables)

            oxidant = {ox: variables[ox] for ox in self._get_name('oxidant')}
            self._set_oxidant_defaults(oxidant)

            reaction = (
                self._calculate_reaction(oxidant,
                                         rho,
                                         variables['ta'],
                                         self._get_name('reactant')))

            # Set Z-coordinate
            if reaction.coords('air_pressure', dim_coords=True):
                self.z_coord = 'air_pressure'
                z_coord = reaction.coord('air_pressure', dim_coords=True)
                z_coord.attributes['positive'] = 'down'
                z_coord.convert_units('Pa')
            elif reaction.coords('atmosphere_hybrid_sigma_pressure_coordinate',
                                 dim_coords=True):
                self.z_coord = 'atmosphere_hybrid_sigma_pressure_coordinate'
                z_coord = reaction.coord(
                    'atmosphere_hybrid_sigma_pressure_coordinate',
                    dim_coords=True)
                z_coord.attributes['positive'] = 'down'
            else:
                raise NotImplementedError(
                    "Lifetime calculation is not implemented"
                    " for the present type of vertical coordinate."
                )

            if not set(['TROP', 'STRA']).isdisjoint(self.cfg['regions']):

                # calculate climatological tropopause pressure (tp_clim)
                # but only if no tropopause is given by data
                if (
                        'ptp' not in variables
                        and 'tp_i' not in variables):
                    tropopause = climatological_tropopause(
                        variables['ta'][:, 0, :, :])

                # If z_coord is defined as:
                #     - air_pressure, use:
                #          - ptp and air_pressure
                #          - tp_clim and air_pressure
                #     - atmosphere_hybrid_sigma_pressure_coordinate, use:
                #          - tp_i and model_level_number
                #          - ptp and (derived) air_pressure
                #          - tp_clim and (derived) air_pressure
                use_z_coord = 'air_pressure'
                if z_coord.name() == 'air_pressure':
                    tropopause = variables.get('ptp', tropopause)
                elif (z_coord.name()
                      == 'atmosphere_hybrid_sigma_pressure_coordinate'):
                    if 'tp_i' in variables:
                        tropopause = variables['tp_i']
                        use_z_coord = 'model_level_number'
                    else:
                        tropopause = variables.get('ptp', tropopause)

            weight = self._define_weight(variables)

            if reaction.coords('atmosphere_hybrid_sigma_pressure_coordinate',
                               dim_coords=True):
                add_model_level(weight)
                add_model_level(reaction)

            input_data_dataset[dataset]['z_coord'] = z_coord
            input_data_dataset[dataset]['use_z_coord'] = use_z_coord
            input_data_dataset[dataset]['tropopause'] = tropopause
            input_data_dataset[dataset]['variables'] = variables
            input_data_dataset[dataset]['reaction'] = reaction
            input_data_dataset[dataset]['weight'] = weight

        return input_data_dataset

    def _set_oxidant_defaults(self, oxidant):
        """Set the defaults for the oxidant reaction rates."""
        for name, _ in oxidant.items():
            # set default reaction
            if name == 'oh':
                self.cfg['oxidant'][name].setdefault('A', 1.85e-20)
                self.cfg['oxidant'][name].setdefault('ER', 987.0)
                self.cfg['oxidant'][name].setdefault('b', 2.82)
            else:
                if (
                        'A' not in self.cfg['oxidant'][name]
                        or 'ER' not in self.cfg['oxidant'][name]
                ):
                    raise KeyError("No sufficient reaction coefficients"
                                   " are given")
                self.cfg['oxidant'][name].setdefault('b', None)

    def _get_all_datasets(self, input_data):
        """Return (sorted) list of datasets."""
        datasets = []
        for data in input_data:
            datasets.append(data[self.cfg['facet_used_for_labels']])
        datasets = list(set(datasets))
        datasets.sort()
        return datasets

    def _get_dataset_data(self, dataset):
        """Return input data corresponding to the chosen dataset."""
        input_data = []
        for data in self.input_data:
            if data[self.cfg['facet_used_for_labels']] == dataset:
                input_data.append(data)

        return input_data

    def _get_name(self, case='reactant'):
        """Return variable name.

        Return the name of the reactant or the oxidant
        respecively.
        """
        if isinstance(self.cfg[case], dict):
            name = []
            for key in self.cfg[case].keys():
                name.append(key)
        else:
            name = self.cfg[case]

        return name

    def _calculate_reaction(self, oxidant, rho, temp, name_reactant):
        """Calculate product of reaction rate and oxidant."""
        reaction = 0.
        for name, oxid in oxidant.items():
            oxid_in_molec = oxid * rho
            reaction_rate = calculate_reaction_rate(
                temp,
                f"{name_reactant.upper()}+{name.upper()}",
                self.cfg['oxidant'][name]['A'],
                self.cfg['oxidant'][name]['ER'],
                self.cfg['oxidant'][name]['b'])
            reaction = (reaction
                        + reaction_rate * oxid_in_molec)
        # test for correct units
        if not reaction.units == 's-1':
            raise ValueError("The units of the reaction rate is"
                             " not consistent. Check input variables.")

        return reaction

    def _define_weight(self, variables):
        """Define used weights in the lifetime calculation.

        Currently only one weight type is implemented. Any other
        options result in weights = 1. (no weighting)

        weight type:
         mass CH4: mass of CH4 -> convert to kg per gridbox
        """
        if self.cfg['weight_type'] == 'mass CH4':
            weight = self._convert_to_mass('ch4', variables)
        else:
            weight = 1.

        return weight

    def _convert_to_mass(self, varname, variables):
        """Convert to kg per gridbox.

        Used constants:
        m_var Molarmass of current species
        """
        # molar mass constants
        m_var = self.cfg['molarmass']

        if 'grmassdry' in variables:
            grmassdry = variables['grmassdry']
        else:
            hus = variables['hus']
            if 'press' in variables:
                press = variables['press']
            else:
                # create 4D pressure variable from
                # pressure (auxilliary) coordinate
                press = create_press(hus)

            grmassdry = calculate_gridmassdry(press,
                                              hus,
                                              self.z_coord)

        var = variables[varname] * grmassdry * (m_var / self.cfg['m_air'])
        # var = da.multiply(da.multiply(variables[varname], grmassdry),
        #                   da.divide(m_var, self.cfg['m_air']))

        return var

    def plot_zonalmean_with_ref(self, plot_func, region,
                                dataset, ref_dataset):
        """Plot zonal mean profile for single dataset with reference."""
        plot_type = 'zonalmean'
        logger.info("Plotting zonal mean profile with reference dataset"
                    " '%s' for '%s'",
                    self._get_label(ref_dataset), self._get_label(dataset))

        # Make sure that the data has the correct dimensions
        cube = calculate_lifetime(dataset,
                                  plot_type,
                                  region)
        ref_cube = calculate_lifetime(ref_dataset,
                                      plot_type,
                                      region)

        # convert units
        cube.convert_units(self.info['units'])

        # Create single figure with multiple axes
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg['figure_kwargs'])
            gridspec = GridSpec(5, 4, figure=fig,
                                height_ratios=[1.0, 1.0, 0.4, 1.0, 1.0])

            # Options used for all subplots
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            fontsize = self.plots[plot_type]['fontsize']

            # Plot dataset (top left)
            axes_data = fig.add_subplot(gridspec[0:2, 0:2])
            plot_kwargs['axes'] = axes_data
            plot_data = plot_func(cube, **plot_kwargs)
            axes_data.set_title(self._get_label(dataset), pad=3.0)
            z_coord = cube.coord(axis='Z')
            axes_data.set_ylabel(f'{z_coord.long_name} [{z_coord.units}]')
            if self.plots[plot_type]['log_y']:
                axes_data.set_yscale('log')
                axes_data.get_yaxis().set_major_formatter(
                    FormatStrFormatter('%.1f'))
            if self.plots[plot_type]['show_y_minor_ticklabels']:
                axes_data.get_yaxis().set_minor_formatter(
                    FormatStrFormatter('%.1f'))
            else:
                axes_data.get_yaxis().set_minor_formatter(NullFormatter())

            # Plot reference dataset (top right)
            # Note: make sure to use the same vmin and vmax than the top left
            # plot if a common colorbar is desired
            axes_ref = fig.add_subplot(gridspec[0:2, 2:4], sharex=axes_data,
                                       sharey=axes_data)
            plot_kwargs['axes'] = axes_ref
            if self.plots[plot_type]['common_cbar']:
                plot_kwargs.setdefault('vmin', plot_data.get_clim()[0])
                plot_kwargs.setdefault('vmax', plot_data.get_clim()[1])
            plot_ref = plot_func(ref_cube, **plot_kwargs)
            axes_ref.set_title(self._get_label(ref_dataset), pad=3.0)
            plt.setp(axes_ref.get_yticklabels(), visible=False)

            # Add colorbar(s)
            self._add_colorbar(plot_type, plot_data, plot_ref, axes_data,
                               axes_ref, dataset, ref_dataset)

            # Plot bias (bottom center)
            bias_cube = cube - ref_cube
            axes_bias = fig.add_subplot(gridspec[3:5, 1:3], sharex=axes_data,
                                        sharey=axes_data)
            plot_kwargs_bias = self._get_plot_kwargs(plot_type, dataset,
                                                     bias=True)
            plot_kwargs_bias['axes'] = axes_bias
            plot_bias = plot_func(bias_cube, **plot_kwargs_bias)
            axes_bias.set_title(
                f"{self._get_label(dataset)} - {self._get_label(ref_dataset)}",
                pad=3.0,
            )
            axes_bias.set_xlabel('latitude [°N]')
            axes_bias.set_ylabel(f'{z_coord.long_name} [{z_coord.units}]')
            cbar_kwargs_bias = self._get_cbar_kwargs(plot_type, bias=True)
            cbar_bias = fig.colorbar(plot_bias, ax=axes_bias,
                                     **cbar_kwargs_bias)
            cbar_bias.set_label(
                self._get_cbar_label(plot_type, dataset, bias=True),
                fontsize=fontsize,
            )
            cbar_bias.ax.tick_params(labelsize=fontsize)

            # Customize plot
            fig.suptitle(f"{dataset['long_name']} ({dataset['start_year']}-"
                         f"{dataset['end_year']})")
            self._process_pyplot_kwargs(plot_type, dataset)

            # Rasterization
            if self.plots[plot_type]['rasterize']:
                self._set_rasterized([axes_data, axes_ref, axes_bias])

        # File paths
        plot_path = self.get_plot_path(plot_type, dataset)
        netcdf_path = (
            get_diagnostic_filename(Path(plot_path).stem + "_{pos}", self.cfg)
        )
        netcdf_paths = {
            netcdf_path.format(pos='top_left'): cube,
            netcdf_path.format(pos='top_right'): ref_cube,
            netcdf_path.format(pos='bottom'): bias_cube,
        }

        return (plot_path, netcdf_paths)

    def plot_zonalmean_without_ref(self, plot_func, region, dataset,
                                   base_datasets, label):
        """Plot zonal mean profile for single dataset without reference."""
        plot_type = 'zonalmean'
        logger.info("Plotting zonal mean profile without reference dataset"
                    " for '%s'", label)

        # zonalmean lifetime is calculated for each time step
        # (sum over longitude) and then a mean is calculated
        cube = calculate_lifetime(dataset,
                                  plot_type,
                                  region)
        # lifetime is averaged over time
        cube = cube.collapsed(['time'], iris.analysis.MEAN)

        # convert units
        cube.convert_units(self.info['units'])

        # Create plot with desired settings
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg['figure_kwargs'])
            axes = fig.add_subplot()
            plot_kwargs = self._get_plot_kwargs(plot_type, base_datasets)
            plot_kwargs['axes'] = axes
            plot_zonalmean = plot_func(cube, **plot_kwargs)

            # Setup colorbar
            fontsize = self.plots[plot_type]['fontsize']
            colorbar = fig.colorbar(plot_zonalmean, ax=axes,
                                    **self._get_cbar_kwargs(plot_type))
            colorbar.set_label(self._get_cbar_label(
                plot_type, self.info),
                fontsize=fontsize)
            colorbar.ax.tick_params(labelsize=fontsize)

            # Customize plot
            axes.set_title(label)
            fig.suptitle(f"{self.info['long_name']} ({dataset['start_year']}-"
                         f"{dataset['end_year']})")
            axes.set_xlabel('latitude [°N]')
            z_coord = cube.coord(name_or_coord='air_pressure')
            axes.set_ylabel(f'{z_coord.long_name} [{z_coord.units}]')
            if self.plots[plot_type]['log_y']:
                axes.set_yscale('log')
                axes.get_yaxis().set_major_formatter(
                    FormatStrFormatter('%.1f'))
            if self.plots[plot_type]['show_y_minor_ticklabels']:
                axes.get_yaxis().set_minor_formatter(
                    FormatStrFormatter('%.1f'))
            else:
                axes.get_yaxis().set_minor_formatter(NullFormatter())
            self._process_pyplot_kwargs(plot_type, dataset)

            # Rasterization
            if self.plots[plot_type]['rasterize']:
                self._set_rasterized([axes])

        # File paths
        plot_path = self.get_plot_path(plot_type, dataset)
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)

        return (plot_path, {netcdf_path: cube})

    def _process_pyplot_kwargs(self, plot_type, dataset):
        """Process functions for :mod:`matplotlib.pyplot`."""
        pyplot_kwargs = self.plots[plot_type]['pyplot_kwargs']
        for (func, arg) in pyplot_kwargs.items():
            if isinstance(arg, str):
                arg = self._fill_facet_placeholders(
                    arg,
                    dataset,
                    f"pyplot_kwargs of {plot_type} '{func}: {arg}'",
                )
            if arg is None:
                getattr(plt, func)()
            else:
                getattr(plt, func)(arg)

    def _check_cube_dimensions(self, cube, plot_type):
        """Check that cube has correct dimensional variables."""
        expected_dimensions_dict = {
            'annual_cycle': (['month_number'],),
            'map': (['latitude', 'longitude'],),
            'zonalmean': (['latitude', self.z_coord],),
            'timeseries': (['time'],),
            '1d_profile': ([self.z_coord],),

        }
        if plot_type not in expected_dimensions_dict:
            raise NotImplementedError(f"plot_type '{plot_type}' not supported")
        expected_dimensions = expected_dimensions_dict[plot_type]
        for dims in expected_dimensions:
            cube_dims = [cube.coords(dim, dim_coords=True) for dim in dims]
            if all(cube_dims) and cube.ndim == len(dims):
                return dims
        expected_dims_str = ' or '.join(
            [str(dims) for dims in expected_dimensions]
        )
        raise ValueError(
            f"Expected cube that exactly has the dimensional coordinates "
            f"{expected_dims_str}, got {cube.summary(shorten=True)}")

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
                        datasets[d].get('reference_for_monitor_diags', False)]
        if len(ref_datasets) > 1:
            raise ValueError(
                f"Expected at most 1 reference dataset (with "
                f"'reference_for_monitor_diags: true' for variable "
                f"'{short_name}', got {len(ref_datasets):d}")
        if ref_datasets:
            return ref_datasets[0]
        return None

    def create_timeseries_plot(self, region, input_data, base_datasets):
        """Create time series plot."""
        plot_type = 'timeseries'
        if plot_type not in self.plots:
            return

        if not base_datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        logger.info("Plotting %s", plot_type)
        fig = plt.figure(**self.cfg['figure_kwargs'])
        axes = fig.add_subplot()

        # Plot all datasets in one single figure
        ancestors = []
        cubes = {}

        for label, dataset in input_data.items():
            ancestors.extend(variable['filename']
                             for variable in dataset['dataset_data'])

            # call by timestep will take longer,
            # but it is less memory intensive
            if self.plots[plot_type]['by_timestep']:
                slice_dataset = {}
                slice_dataset['z_coord'] = dataset['z_coord']
                slice_dataset['use_z_coord'] = dataset['use_z_coord']
                cube_slices = iris.cube.CubeList()
                for (reaction_slice,
                     weight_slice,
                     tp_slice) in zip(dataset['reaction'].slices_over('time'),
                                      dataset['weight'].slices_over('time'),
                                      dataset['tropopause'].slices_over(
                                          'time')):

                    slice_dataset['reaction'] = reaction_slice
                    slice_dataset['weight'] = weight_slice
                    slice_dataset['tropopause'] = tp_slice

                    cube_slices.append(calculate_lifetime(slice_dataset,
                                                          plot_type,
                                                          region))
                cube = cube_slices.merge_cube()
            else:
                cube = calculate_lifetime(dataset,
                                          plot_type,
                                          region)

            # convert units
            cube.convert_units(self.info['units'])

            cubes[label] = cube
            self._check_cube_dimensions(cube, plot_type)

            # Plot original time series
            plot_kwargs = self._get_plot_kwargs(plot_type,
                                                base_datasets[label])
            plot_kwargs['axes'] = axes

            if self.plots[plot_type]['display_mean'] is not False:
                mean = cube.collapsed('time', iris.analysis.MEAN).data
                plot_kwargs['label'] = (f"{plot_kwargs['label']}"
                                        f" ({mean:.2f})")

            # Plot annual means if desired
            annual_mean = self.plots[plot_type]['annual_mean']
            if annual_mean in [False, 'both']:
                iris.plot.plot(cube, **plot_kwargs)
                plot_kwargs.pop('label', None)
            elif annual_mean in ['both', 'only']:
                logger.debug("Plotting annual means")
                if not cube.coords('year'):
                    add_year(cube, 'time')
                annual_mean_cube = cube.aggregated_by('year',
                                                      iris.analysis.MEAN)

                plot_kwargs.update(self.plots[plot_type]['annual_mean_kwargs'])
                iris.plot.plot(annual_mean_cube, **plot_kwargs)
            else:
                raise ValueError("Unknown option for annual_mean."
                                 "Choose False, 'both', or 'only'.")

        # Default plot appearance
        multi_dataset_facets = self._get_multi_dataset_facets(
            list(base_datasets.values()))
        axes.set_title(f'{self.info["long_name"]} in region {region}')
        axes.set_xlabel('Time')
        axes.set_ylabel(f"{chr(964)}({self._get_name('reactant').upper()})"
                        f" [{self.info['units']}]")
        gridline_kwargs = self._get_gridline_kwargs(plot_type)
        if gridline_kwargs is not False:
            axes.grid(**gridline_kwargs)

        # Legend
        legend_kwargs = self.plots[plot_type]['legend_kwargs']
        if legend_kwargs is not False:
            axes.legend(**legend_kwargs)

        # Customize plot appearance
        self._process_pyplot_kwargs(plot_type, multi_dataset_facets)

        # Save plot
        plot_path = self.get_plot_path(plot_type, multi_dataset_facets)
        fig.savefig(plot_path, **self.cfg['savefig_kwargs'])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save netCDF file
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)
        var_attrs = {
            'short_name': self.info['short_name'],
            'long_name': self.info['long_name'],
            'units': self.info['units']
        }
        io.save_1d_data(cubes, netcdf_path, 'time', var_attrs)

        # Provenance tracking
        caption = (f"Time series of {multi_dataset_facets['long_name']} for "
                   f"various datasets.")
        provenance_record = {
            'ancestors': ancestors,
            'authors': ['schlund_manuel', 'winterstein_franziska'],
            'caption': caption,
            'plot_types': ['line'],
            'long_names': [var_attrs['long_name']],
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)
            provenance_logger.log(netcdf_path, provenance_record)

    def create_annual_cycle_plot(self, region, input_data, base_datasets):
        """Create annual cycle plot."""
        plot_type = 'annual_cycle'
        if plot_type not in self.plots:
            return

        if not base_datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        logger.info("Plotting %s", plot_type)
        fig = plt.figure(**self.cfg['figure_kwargs'])
        axes = fig.add_subplot()

        # Plot all datasets in one single figure
        ancestors = []
        cubes = {}
        for label, dataset in input_data.items():
            ancestors.extend(variable['filename']
                             for variable in dataset['dataset_data'])

            cube = calculate_lifetime(dataset,
                                      plot_type,
                                      region)
            # convert units
            cube.convert_units(self.info['units'])

            cubes[label] = cube
            self._check_cube_dimensions(cube, plot_type)

            # Plot annual cycle
            plot_kwargs = self._get_plot_kwargs(plot_type,
                                                base_datasets[label])
            plot_kwargs['axes'] = axes
            iris.plot.plot(cube, **plot_kwargs)

        # Default plot appearance
        multi_dataset_facets = self._get_multi_dataset_facets(
            list(base_datasets.values()))
        axes.set_title(f'{self.info["long_name"]} in region {region}')
        axes.set_xlabel('Month')
        axes.set_ylabel(f"$\tau$({self._get_name('reactant').upper()})"
                        " [{self.info['units']}]")
        axes.set_xticks(range(1, 13), [str(m) for m in range(1, 13)])
        gridline_kwargs = self._get_gridline_kwargs(plot_type)
        if gridline_kwargs is not False:
            axes.grid(**gridline_kwargs)

        # Legend
        legend_kwargs = self.plots[plot_type]['legend_kwargs']
        if legend_kwargs is not False:
            axes.legend(**legend_kwargs)

        # Customize plot appearance
        self._process_pyplot_kwargs(plot_type, multi_dataset_facets)

        # Save plot
        plot_path = self.get_plot_path(plot_type, multi_dataset_facets)
        fig.savefig(plot_path, **self.cfg['savefig_kwargs'])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save netCDF file
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)
        var_attrs = {
            'short_name': self.info['short_name'],
            'long_name': self.info['long_name'],
            'units': self.info['units']
        }
        io.save_1d_data(cubes, netcdf_path, 'month_number', var_attrs)

        # Provenance tracking
        caption = (f"Annual cycle of {self.info['long_name']} for "
                   f"various datasets.")
        provenance_record = {
            'ancestors': ancestors,
            'authors': ['schlund_manuel', 'winterstein_franziska'],
            'caption': caption,
            'plot_types': ['seas'],
            'long_names': [var_attrs['long_name']],
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)
            provenance_logger.log(netcdf_path, provenance_record)

    def create_zonalmean_plot(self, region, input_data,
                              base_datasets):
        """Create zonal mean profile plot."""
        plot_type = 'zonalmean'
        if plot_type not in self.plots:
            return

        if not base_datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        # Get reference dataset if possible
        ref_dataset = self._get_reference_dataset(base_datasets, 'lifetime')
        if ref_dataset is None:
            logger.info("Plotting %s without reference dataset", plot_type)
        else:
            logger.info("Plotting %s with reference dataset '%s'", plot_type,
                        self._get_label(ref_dataset))

        # Get plot function
        plot_func = self._get_plot_func(plot_type)

        # Create a single plot for each dataset (incl. reference dataset if
        # given)
        ancestors = []
        for label, dataset in input_data.items():
            if dataset == ref_dataset:
                continue
            ancestors.extend(variable['filename']
                             for variable in dataset['dataset_data'])
            if ref_dataset is None:
                (plot_path, netcdf_paths) = (
                    self.plot_zonalmean_without_ref(
                        plot_func, region,
                        dataset, base_datasets[label], label)
                )
                caption = (
                    f"Zonal mean profile of {dataset['long_name']} of dataset "
                    f"{dataset['dataset']} (project {dataset['project']}) "
                    f"from {dataset['start_year']} to {dataset['end_year']}."
                )
            else:
                (plot_path, netcdf_paths) = (
                    self.plot_zonalmean_with_ref(plot_func, region,
                                                 dataset, ref_dataset)
                )
                caption = (
                    f"Zonal mean profile of {dataset['long_name']} of dataset "
                    f"{dataset['dataset']} (project {dataset['project']}) "
                    f"including bias relative to {ref_dataset['dataset']} "
                    f"(project {ref_dataset['project']}) from "
                    f"{dataset['start_year']} to {dataset['end_year']}."
                )
                ancestors.append(ref_dataset['filename'])

            # If statistics are shown add a brief description to the caption
            if self.plots[plot_type]['show_stats']:
                caption += (
                    " The number in the top left corner corresponds to the "
                    "spatial mean (weighted by grid cell areas).")

            # Save plot
            plt.savefig(plot_path, **self.cfg['savefig_kwargs'])
            logger.info("Wrote %s", plot_path)
            plt.close()

            # Save netCDFs
            for (netcdf_path, cube) in netcdf_paths.items():
                io.iris_save(cube, netcdf_path)

            # Provenance tracking
            provenance_record = {
                'ancestors': ancestors,
                'authors': ['schlund_manuel', 'winterstein_franziska'],
                'caption': caption,
                'plot_types': ['vert'],
                'long_names': [dataset['long_name']],
            }
            with ProvenanceLogger(self.cfg) as provenance_logger:
                provenance_logger.log(plot_path, provenance_record)
                for netcdf_path in netcdf_paths:
                    provenance_logger.log(netcdf_path, provenance_record)

    def create_1d_profile_plot(self, region, input_data, base_datasets):
        """Create 1D profile plot."""
        plot_type = '1d_profile'
        if plot_type not in self.plots:
            return

        if not base_datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        logger.info("Plotting %s", plot_type)
        fig = plt.figure(**self.cfg['figure_kwargs'])
        axes = fig.add_subplot()

        # Plot all datasets in one single figure
        ancestors = []
        cubes = {}
        for label, dataset in input_data.items():
            ancestors.extend(variable['filename']
                             for variable in dataset['dataset_data'])

            cube = calculate_lifetime(dataset,
                                      plot_type,
                                      region)
            # convert units
            cube.convert_units(self.info['units'])

            cubes[label] = cube
            self._check_cube_dimensions(cube, plot_type)

            # Plot 1D profile
            plot_kwargs = self._get_plot_kwargs(plot_type,
                                                base_datasets[label])
            plot_kwargs['axes'] = axes

            iris.plot.plot(cube, **plot_kwargs)

        # Default plot appearance
        multi_dataset_facets = self._get_multi_dataset_facets(
            list(base_datasets.values()))
        axes.set_title(f'{self.info["long_name"]} in region {region}')
        axes.set_xlabel(f"$\tau$({self._get_name('reactant').upper()})"
                        f" [{self.info['units']}]")
        z_coord = cube.coord(axis='Z')
        axes.set_ylabel(f'{z_coord.long_name} [{z_coord.units}]')

        # apply logarithmic axes
        if self.plots[plot_type]['log_y']:
            axes.set_yscale('log')
            axes.get_yaxis().set_major_formatter(
                FormatStrFormatter('%.1f'))
        if self.plots[plot_type]['show_y_minor_ticklabels']:
            axes.get_yaxis().set_minor_formatter(
                FormatStrFormatter('%.1f'))
        else:
            axes.get_yaxis().set_minor_formatter(NullFormatter())
        if self.plots[plot_type]['log_x']:
            axes.set_xscale('log')
            # major and minor ticks
            x_major = LogLocator(base=10.0, numticks=12)
            axes.get_xaxis().set_major_locator(x_major)
            x_minor = LogLocator(base=10.0,
                                 subs=np.arange(1.0, 10.0) * 0.1,
                                 numticks=12)

            axes.get_xaxis().set_minor_locator(x_minor)
            axes.get_xaxis().set_minor_formatter(NullFormatter())

        # gridlines
        gridline_kwargs = self._get_gridline_kwargs(plot_type)
        if gridline_kwargs is not False:
            axes.grid(**gridline_kwargs)
        # nicer aspect ratio
        aspect_ratio = self.plots[plot_type]['aspect_ratio']
        axes.set_box_aspect(aspect_ratio)

        # Legend
        legend_kwargs = self.plots[plot_type]['legend_kwargs']
        if legend_kwargs is not False:
            axes.legend(**legend_kwargs)

        # Customize plot appearance
        self._process_pyplot_kwargs(plot_type, multi_dataset_facets)

        # Save plot
        plot_path = self.get_plot_path(plot_type, multi_dataset_facets)
        fig.savefig(plot_path, **self.cfg['savefig_kwargs'])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save netCDF file
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)
        var_attrs = {
            'short_name': self.info['short_name'],
            'long_name': self.info['long_name'],
            'units': self.info['units']
        }
        io.save_1d_data(cubes, netcdf_path, z_coord.standard_name, var_attrs)

        # Provenance tracking
        caption = ("Vertical one-dimensional profile of "
                   f"{self.info['long_name']}"
                   " for various datasets.")
        provenance_record = {
            'ancestors': ancestors,
            'authors': ['schlund_manuel', 'winterstein_franziska'],
            'caption': caption,
            'plot_types': ['line'],
            'long_names': [var_attrs['long_name']],
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)
            provenance_logger.log(netcdf_path, provenance_record)

    def compute(self):
        """Plot preprocessed data."""
        input_data = self.input_data_dataset
        base_datasets = {label: dataset['dataset_data'][0]
                         for label, dataset in input_data.items()}

        # at the moment regions only apply to TROP and STRAT
        for region in self.cfg['regions']:
            logger.info("Plotting lifetime for region %s", region)
            self.create_timeseries_plot(region, input_data, base_datasets)
            self.create_annual_cycle_plot(region, input_data, base_datasets)

        self.create_zonalmean_plot(region, input_data,
                                   base_datasets)
        self.create_1d_profile_plot(region, input_data, base_datasets)


def main():
    """Run diagnostic."""
    with run_diagnostic() as config:
        CH4Lifetime(config).compute()


if __name__ == '__main__':
    main()
