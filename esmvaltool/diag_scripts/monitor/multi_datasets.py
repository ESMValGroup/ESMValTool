#!/usr/bin/env python
"""Monitoring diagnostic to visualize multiple datasets.

Description
-----------
This diagnostic can be used to visualize multiple datasets in a single plot. A
large variety of different plot types is supported.

Currently supported plot types
------------------------------

Plot types can be specified with the :ref:`recipe option
<diag_monitor_multi_datasets_recipe_options>` ``plots``. All plots can be
customized with additional :ref:`plot-specific options
<diag_monitor_multi_datasets_plot_options>`.

- ``1d_profile`` (1D plot): plot 1D profile. Input data needs single dimension
  `altitude` or `air_pressure`.
- ``annual_cycle`` (1D plot): plot annual cycle. Input data needs single
  dimension `month_number`.
- ``benchmarking_annual_cycle`` (1D plot): benchmarking version of
  ``annual_cycle``. Input data needs single dimension `month_number`.
- ``benchmarking_boxplot`` (boxplot): plot boxplot. Input data needs to be
  scalar.
- ``benchmarking_diurnal_cycle`` (1D plot): benchmarking version of
  ``diurnal_cycle``. Input data needs single dimension `hour`.
- ``benchmarking_map`` (2D plot): benchmarking version of ``map``. Input data
  needs dimensions `(longitude, latitude)`.
- ``benchmarking_timeseries`` (1D plot): benchmarking version of
  ``timeseries``. Input data needs single dimension `time`.
- ``benchmarking_zonal`` (2D plot): benchmarking version of
  ``zonal_mean_profile``. Input data needs dimensions `(latitude, altitude)` or
  `latitude, air_pressure`.
- ``diurnal_cycle`` (1D plot): plot diurnal cycle. Input data needs single
  dimension `hour`.
- ``hovmoeller_anncyc_vs_lat_or_lon`` (2D plot): plot Hovmoeller plot. Input
  data needs dimensions `(month_number, latitude)` or `(month_number,
  longitude)`.
- ``hovmoeller_time_vs_lat_or_lon`` (2D plot): plot Hovmoeller plot. Input data
  needs dimensions `(time, latitude)` or `(time, longitude)`.
- ``hovmoeller_z_vs_time`` (2D plot): plot Hovmoeller plot. Input data needs
  dimensions `(time, altitude)` or `(time, air_pressure)`.
- ``map`` (2D plot): plot map plot. Input data needs dimensions `(longitude,
  latitude)`.
- ``timeseries`` (1D plot): plot time series. Input data needs single dimension
  `time`.
- ``variable_vs_lat`` (1D plot): plot variable vs. latitude. Input data needs
  single dimension `latitude`.
- ``zonal_mean_profile`` (2D plot): plot zonal mean profile. Input data needs
  dimensions `(latitude, altitude)` or `latitude, air_pressure`.

Each plot type belongs to one or more classes which define how the data is
plotted:

- Class *1D plot*: for each variable separately, all datasets are plotted in
  one single figure. Input data needs to be 1D.
- Class *2D plot*: for each variable and dataset, an individual figure is
  plotted. Input data needs to be 2D. A single reference dataset can be defined
  by setting the facet ``reference_for_monitor_diags: True`` in the dataset
  definition in the recipe. In this case, three panels are plotted, incl. a
  bias. Note that if a reference dataset is defined, all input datasets need to
  be given on the same horizontal and vertical grid (you can use the
  preprocessors :func:`esmvalcore.preprocessor.regrid` and
  :func:`esmvalcore.preprocessor.extract_levels` for this).
- Class *boxplot*: one figure for all datasets is plotted. Input data needs to
  be 0D.

In addition, for many plot types, a benchmarking version exists
(``benchmarking_*``). These plots include envelopes (1D plots) or hatching (2D
plots) to visualize data ranges (e.g., multi-model percentiles) in addition to
the data itself. For these plots, at least one benchmarking dataset needs to be
defined by setting the facet ``benchmark_dataset: True`` in the dataset
definition in the recipe. Reference datasets are not supported for 2D
benchmarking diagnostics.

Author
------
Manuel Schlund (DLR, Germany)

.. _diag_monitor_multi_datasets_recipe_options:

Configuration options in recipe
-------------------------------
facet_used_for_labels: str, optional (default: 'dataset')
    Facet used to label different datasets in plot titles and legends. For
    example, ``facet_used_for_labels: 'dataset'`` will use dataset names in
    plot titles and legends; ``facet_used_for_labels: 'exp'`` will use
    experiments in plot titles and legends. In addition,
    ``facet_used_for_labels`` is used to select the correct ``plot_kwargs`` for
    the different datasets (see configuration options for the different plot
    types below).
figure_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.figure`. By
    default, uses ``{constrained_layout: True}``.
group_variables_by: str, optional (default: 'short_name')
    Facet or coordinate which is used to create variable groups. For each
    variable group, an individual plot is created. Specifying a coordinate
    allows to create one plot for each point along a dimension. For example,
    when used in combination with the preprocessor function
    :func:`esmvalcore.preprocessor.extract_shape` the `shape_id` coordinate
    can be used to create one plot for each shape.
matplotlib_rc_params: dict, optional
    Optional :class:`matplotlib.RcParams` used to customize matplotlib plots.
    Options given here will be passed to :func:`matplotlib.rc_context` and used
    for all plots produced with this diagnostic. Note: fontsizes specified here
    might be overwritten by the plot-type-specific option ``fontsize`` (see
    below).
plots: dict
    Plot types plotted by this diagnostic (see list above). Dictionary keys
    must be elements of the list above.  Dictionary values are dictionaries
    used as options for the corresponding plot. The allowed options for the
    different plot types are given below.
plot_filename: str, optional
    Filename pattern for the plots. By default, uses
    ``'{plot_type}_{real_name}_{dataset}_{mip}_{exp}_{ensemble}'``.  All tags
    (i.e., the entries in curly brackets, e.g., ``'{dataset}'``, are replaced
    with the corresponding tags).
plot_folder: str, optional
    Path to the folder to store figures. By default, uses
    ``'{plot_dir}/../../{dataset}/{exp}/{modeling_realm}/{real_name}'``.  All
    tags (i.e., the entries in curly brackets, e.g., ``'{dataset}'``, are
    replaced with the corresponding tags). ``'{plot_dir}'`` is replaced with
    the default ESMValTool plot directory (i.e.,
    ``output_dir/plots/diagnostic_name/script_name/``, see
    :ref:`esmvalcore:outputdata`).
savefig_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.savefig`. By
    default, uses ``{bbox_inches: 'tight', dpi: 300, orientation:
    'landscape'}``.
seaborn_settings: dict, optional
    Options for :func:`seaborn.set_theme` (affects all plots). By default, uses
    ``{style: 'ticks'}``.

.. _diag_monitor_multi_datasets_plot_options:

Configuration options for 1D plots
----------------------------------
aspect_ratio: float, optional (default: None)
    Aspect ratio of the plot.
caption: str, optional
    Figure caption used for provenance tracking. Can include facets in curly
    brackets which will be derived from the corresponding dataset, e.g.,
    ``{project}``, ``{short_name}``, ``{exp}``. By default, uses a very basic
    caption.
envelope_kwargs: dict, optional
    Optional keyword arguments for :func:`iris.plot.fill_between`. By default,
    uses ``{alpha: 0.8, facecolor: 'lightblue', linewidth: 0.0, zorder: 1.0}``.
    Only relevant for benchmarking plots.
gridline_kwargs: dict, optional
    Optional keyword arguments for grid lines. By default, uses ``{color:
    'lightgrey', alpha: 0.5}``. Use ``gridline_kwargs: False`` to not show grid
    lines.
hlines: list of dict, optional
    Horizontal lines to show in plot. Each list element corresponds to one
    line, and each list element should contain a dictionary with keywords
    arguments passed to :meth:`matplotlib.axes.Axes.axhline`. Example: ``[{y:
    0}, {y: 1, color: 'red'}]``.
legend_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.legend`. Use
    ``legend_kwargs: False`` to not show legends.
log_x: bool, optional (default: False)
    Use logarithmic X-axis.
log_y: bool, optional (default: False)
    Use logarithmic Y-axis.
plot_kwargs: dict, optional
    Optional keyword arguments for :func:`iris.plot.plot`. Dictionary keys are
    elements identified by ``facet_used_for_labels`` or ``'default'``, e.g.,
    ``'CMIP6'`` if ``facet_used_for_labels: 'project'`` or ``'historical'`` if
    ``facet_used_for_labels: 'exp'``. Dictionary values are dictionaries used
    as keyword arguments for :func:`iris.plot.plot`. String arguments can
    include facets in curly brackets which will be derived from the
    corresponding dataset, e.g., ``{project}``, ``{short_name}``, ``{exp}``.
    Examples: ``{default: {linestyle: '-', label: '{project}'}, CMIP6: {color:
    'red', linestyle: '--'}, OBS: {color: 'black'}}``.
pyplot_kwargs: dict, optional
    Optional calls to functions of :mod:`matplotlib.pyplot`. Dictionary keys
    are functions of :mod:`matplotlib.pyplot`. Dictionary values are used as
    argument(s) for these functions (if values are dictionaries, these are
    interpreted as keyword arguments; otherwise a single argument is assumed).
    String arguments can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``. Examples: ``{title: 'Awesome Plot of {long_name}'}``, ``{xlabel:
    '{short_name}'}``, ``{xlim: [0, 5]}``.
rasterize: bool, optional (default: False)
    If ``True``, use rasterization_ for plots to produce smaller files.  This
    is only relevant for vector graphics (e.g., ``output_file_type: 'pdf'``).
time_format: str, optional (default: None)
    :func:`~datetime.datetime.strftime` format string that is used to format
    the time axis using :class:`matplotlib.dates.DateFormatter`. If ``None``,
    use the default formatting imposed by the iris plotting function.
transpose_axes: bool, optional (default: False)
    Swap X- and Y-axis.
x_major_formatter: str, optional (default: None)
    Format string for :class:`matplotlib.ticker.FormatStrFormatter` used to
    format major tick labels of X-axis.
x_minor_formatter: str, optional (default: None)
    Format string for :class:`matplotlib.ticker.FormatStrFormatter` used to
    format minor tick labels of X-axis.
y_major_formatter: str, optional (default: None)
    Format string for :class:`matplotlib.ticker.FormatStrFormatter` used to
    format major tick labels of Y-axis.
y_minor_formatter: str, optional (default: None)
    Format string for :class:`matplotlib.ticker.FormatStrFormatter` used to
    format minor tick labels of Y-axis.

Configuration options for 2D plots
----------------------------------
aspect_ratio: float, optional (default: None)
    Aspect ratio of the plot.
caption: str, optional
    Figure caption used for provenance tracking. Can include facets in curly
    brackets which will be derived from the corresponding dataset, e.g.,
    ``{project}``, ``{short_name}``, ``{exp}``. By default, uses a very basic
    caption.
cbar_label: str, optional (default: '{short_name} [{units}]')
    Colorbar label. Can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``.
cbar_label_bias: str, optional (default: 'Δ{short_name} [{units}]')
    Colorbar label for plotting biases. Can include facets in curly brackets
    which will be derived from the corresponding dataset, e.g., ``{project}``,
    ``{short_name}``, ``{exp}``. Only relevant for plots including reference
    datasets.
cbar_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.colorbar`. By
    default, uses ``{orientation: 'vertical'}``.
cbar_kwargs_bias: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.colorbar` for
    plotting biases. These keyword arguments update (and potentially overwrite)
    the ``cbar_kwargs`` for the bias plot. Only relevant for plots including
    reference datasets.
common_cbar: bool, optional (default: False)
    Use a common colorbar for the top panels (i.e., plots of the dataset and
    the corresponding reference dataset) when using a reference dataset. If
    neither ``vmin`` and ``vmax`` nor ``levels`` is given in ``plot_kwargs``,
    the colorbar bounds are inferred from the dataset in the top left panel,
    which might lead to an inappropriate colorbar for the reference dataset
    (top right panel). Thus, specify ``vmin`` and ``vmax`` or ``levels`` when
    using ``common_cbar: True``. Only relevant for plots including reference
    datasets.
fontsize: int, optional (default: None)
    Fontsize used for ticks, labels and titles. For the latter, use the given
    fontsize plus 2. Does not affect suptitles. If not given, use default
    matplotlib values. For a more fine-grained definition of fontsizes, use the
    option ``matplotlib_rc_params`` (see above).
gridline_kwargs: dict, optional (default: False)
    Optional keyword arguments for grid lines. Use ``gridline_kwargs: False``
    to not show grid lines.
log_x: bool, optional (default: False)
    Use logarithmic X-axis.
log_y: bool, optional (default: False)
    Use logarithmic Y-axis.
plot_func: str, optional (default: 'contourf')
    Plot function used to plot the maps. Must be a function of :mod:`iris.plot`
    that supports plotting of 2D data.
plot_kwargs: dict, optional
    Optional keyword arguments for the plot function defined by ``plot_func``.
    Dictionary keys are elements identified by ``facet_used_for_labels`` or
    ``'default'``, e.g., ``'CMIP6'`` if ``facet_used_for_labels: 'project'`` or
    ``'historical'`` if ``facet_used_for_labels: 'exp'``. Dictionary values are
    dictionaries used as keyword arguments for the plot function defined by
    ``plot_func``. String arguments can include facets in curly brackets which
    will be derived from the corresponding dataset, e.g., ``{project}``,
    ``{short_name}``, ``{exp}``. Examples: ``{default: {levels: 2}, CMIP6:
    {vmin: 200, vmax: 250}}``. In addition to the normalization_ options
    supported by the plot function, the option ``{norm: 'centered'}`` can be
    specified. In this case, the keywords ``vcenter`` and ``halfrange`` should
    be used instead of ``vmin`` or ``vmax`` (see
    :class:`~matplotlib.colors.CenteredNorm`).
plot_kwargs_bias: dict, optional
    Optional keyword arguments for the plot function defined by ``plot_func``
    for plotting biases. These keyword arguments update (and potentially
    overwrite) the ``plot_kwargs`` for the bias plot. By default, uses ``{cmap:
    'bwr', norm: 'centered'}``. Only relevant for plots including reference
    datasets.
projection: str, optional (default: None)
    Projection used for the plot. Needs to be a valid projection class of
    :mod:`cartopy.crs`. Keyword arguments can be specified using the option
    ``projection_kwargs``. For map plots, ``'Robinson'`` is used as default.
projection_kwargs: dict, optional
    Optional keyword arguments for the projection given by ``projection``. For
    map plots, the default keyword arguments ``{central_longitude: 10}`` are
    used.
pyplot_kwargs: dict, optional
    Optional calls to functions of :mod:`matplotlib.pyplot`. Dictionary keys
    are functions of :mod:`matplotlib.pyplot`. Dictionary values are used as
    argument(s) for these functions (if values are dictionaries, these are
    interpreted as keyword arguments; otherwise a single argument is assumed).
    String arguments can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``. Examples: ``{title: 'Awesome Plot of {long_name}'}``, ``{xlabel:
    '{short_name}'}``, ``{xlim: [0, 5]}``.
rasterize: bool, optional (default: False)
    If ``True``, use rasterization_ for plots to produce smaller files.  This
    is only relevant for vector graphics (e.g., ``output_file_type: 'pdf'``).
show_stats: bool, optional (default: True)
    Show basic statistics on the plots.
time_format: str, optional (default: None)
    :func:`~datetime.datetime.strftime` format string that is used to format
    the time axis using :class:`matplotlib.dates.DateFormatter`. If ``None``,
    use the default formatting imposed by the iris plotting function.
transpose_axes: bool, optional (default: False)
    Swap X- and Y-axis.
x_major_formatter: str, optional (default: None)
    Format string for :class:`matplotlib.ticker.FormatStrFormatter` used to
    format major tick labels of X-axis.
x_minor_formatter: str, optional (default: None)
    Format string for :class:`matplotlib.ticker.FormatStrFormatter` used to
    format minor tick labels of X-axis.
x_pos_stats_avg: float, optional (default: 0.01)
    Text X-position of average (shown on the left) in Axes coordinates. Can be
    adjusted to avoid overlap with the figure. Only relevant if ``show_stats:
    True``.
x_pos_stats_bias: float, optional (default: 0.7)
    Text X-position of bias statistics (shown on the right) in Axes
    coordinates. Can be adjusted to avoid overlap with the figure. Only
    relevant if ``show_stats: True``.
y_major_formatter: str, optional (default: None)
    Format string for :class:`matplotlib.ticker.FormatStrFormatter` used to
    format major tick labels of Y-axis.
y_minor_formatter: str, optional (default: None)
    Format string for :class:`matplotlib.ticker.FormatStrFormatter` used to
    format minor tick labels of Y-axis.

Configuration options for boxplots
----------------------------------
caption: str, optional
    Figure caption used for provenance tracking. Can include facets in curly
    brackets which will be derived from the corresponding dataset, e.g.,
    ``{project}``, ``{short_name}``, ``{exp}``. By default, uses a very basic
    caption.
fontsize: int, optional (default: None)
    Fontsize used for ticks, labels and titles. For the latter, use the given
    fontsize plus 2. Does not affect suptitles. If not given, use default
    matplotlib values. For a more fine-grained definition of fontsizes, use the
    option ``matplotlib_rc_params`` (see above).
plot_kwargs: dict, optional
    Optional keyword arguments for :func:`seaborn.boxplot`.
pyplot_kwargs: dict, optional
    Optional calls to functions of :mod:`matplotlib.pyplot`. Dictionary keys
    are functions of :mod:`matplotlib.pyplot`. Dictionary values are used as
    argument(s) for these functions (if values are dictionaries, these are
    interpreted as keyword arguments; otherwise a single argument is assumed).
    String arguments can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``. Examples: ``{title: 'Awesome Plot of {long_name}'}``, ``{xlabel:
    '{short_name}'}``, ``{xlim: [0, 5]}``.
scatter_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.scatter` (used to
    plot marks for benchmarking dataset). By default, uses ``{marker: 'x', s:
    200, linewidths: 2, color: 'red', zorder: 3}``.
var_order: list of str, optional (default: None)
    Order in which variables are plotted. Must be elements identified via
    ``group_variables_by``. If ``None``, use order in which data is given.

.. _rasterization: https://matplotlib.org/stable/gallery/misc/
   rasterization_demo.html
.. _normalization: https://matplotlib.org/stable/users/explain/colors/
   colormapnorms.html

"""

from __future__ import annotations

import inspect
import logging
import warnings
from collections.abc import Callable, Iterable
from copy import copy, deepcopy
from functools import partial
from pathlib import Path
from pprint import pformat
from typing import Any

import cartopy.crs as ccrs
import dask.array as da
import iris
import iris.analysis
import iris.pandas
import iris.plot
import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from iris.analysis.cartography import area_weights
from iris.coords import AuxCoord, Coord
from iris.cube import Cube, CubeList
from iris.exceptions import ConstraintMismatchError
from matplotlib.axes import Axes
from matplotlib.colors import CenteredNorm
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import (
    AutoMinorLocator,
    FormatStrFormatter,
    LogLocator,
    NullLocator,
)
from sklearn.metrics import r2_score

import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool import ESMValToolDeprecationWarning
from esmvaltool.diag_scripts.monitor.monitor_base import MonitorBase
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    group_metadata,
    io,
    run_diagnostic,
)
from esmvaltool.diag_scripts.shared._base import sorted_metadata

logger = logging.getLogger(Path(__file__).stem)


class MultiDatasets(MonitorBase):
    """Diagnostic to plot multi-dataset plots."""

    @property
    def plot_settings(self) -> dict[str, dict[str, Any]]:
        """Plot settings."""
        default_settings_1d = {
            "aspect_ratio": None,
            "caption": None,
            "envelope_kwargs": {
                "alpha": 0.8,
                "facecolor": "lightblue",
                "linewidth": 0.0,
                "zorder": 1.0,
            },
            "gridline_kwargs": {},
            "hlines": [],
            "legend_kwargs": {},
            "log_x": False,
            "log_y": False,
            "plot_kwargs": {},
            "pyplot_kwargs": {},
            "rasterize": False,
            "time_format": None,
            "transpose_axes": False,
            "x_major_formatter": None,
            "x_minor_formatter": None,
            "y_major_formatter": None,
            "y_minor_formatter": None,
        }
        default_settings_2d = {
            "aspect_ratio": None,
            "caption": None,
            "cbar_label": "{short_name} [{units}]",
            "cbar_label_bias": "Δ{short_name} [{units}]",
            "cbar_kwargs": {"orientation": "vertical"},
            "cbar_kwargs_bias": {},
            "common_cbar": False,
            "fontsize": None,
            "gridline_kwargs": False,
            "log_x": False,
            "log_y": False,
            "plot_func": "contourf",
            "plot_kwargs": {},
            "plot_kwargs_bias": {},
            "projection": None,
            "projection_kwargs": {},
            "pyplot_kwargs": {},
            "rasterize": True,
            "show_stats": True,
            "time_format": None,
            "transpose_axes": False,
            "x_major_formatter": None,
            "x_minor_formatter": None,
            "x_pos_stats_avg": 0.01,
            "x_pos_stats_bias": 0.7,
            "y_major_formatter": None,
            "y_minor_formatter": None,
        }

        # Keys for each plot type:
        # - function: function used to create plot.
        # - coords: expected dimensional coordinates for input data.
        # - provenance: information for provenance tracking.
        # - pyplot_kwargs: additional calls to pyplot functions. Will always be
        #   run, but may be overwritten by user-defined `pyplot_kwargs` given
        #   under `default_settings`.
        # - default_settings: allowed user-defined settings in recipe and their
        #   default values.
        plot_settings = {
            "1d_profile": {
                "function": partial(self.create_1d_plot, "1d_profile"),
                "coords": (["air_pressure"], ["altitude"]),
                "provenance": {
                    "authors": ["schlund_manuel", "winterstein_franziska"],
                    "caption": (
                        "Vertical one-dimensional profile of {long_name} for "
                        "various datasets."
                    ),
                    "plot_types": ["line"],
                },
                "pyplot_kwargs": {},
                "default_settings": {
                    **default_settings_1d,
                    "aspect_ratio": 1.5,
                    "log_y": True,
                    "transpose_axes": True,
                    "y_minor_formatter": "",
                },
            },
            "annual_cycle": {
                "function": partial(self.create_1d_plot, "annual_cycle"),
                "coords": (["month_number"],),
                "provenance": {
                    "authors": ["schlund_manuel"],
                    "caption": (
                        "Annual cycle of {long_name} for various datasets."
                    ),
                    "plot_types": ["seas"],
                },
                "pyplot_kwargs": {
                    "xlabel": "month",
                    "xticks": {
                        "ticks": range(1, 13),
                        "labels": [str(m) for m in range(1, 13)],
                    },
                },
                "default_settings": {**default_settings_1d},
            },
            "benchmarking_annual_cycle": {
                "function": partial(
                    self.create_1d_benchmarking_plot,
                    "benchmarking_annual_cycle",
                ),
                "coords": (["month_number"],),
                "provenance": {
                    "authors": ["lauer_axel", "schlund_manuel"],
                    "caption": (
                        "Annual cycle of {long_name} for various datasets."
                    ),
                    "plot_types": ["seas"],
                },
                "pyplot_kwargs": {
                    "xlabel": "month",
                    "xticks": {
                        "ticks": range(1, 13),
                        "labels": [str(m) for m in range(1, 13)],
                    },
                },
                "default_settings": {**default_settings_1d},
            },
            "benchmarking_boxplot": {
                "function": self.create_benchmarking_boxplot,
                "coords": ([],),
                "provenance": {
                    "authors": ["bock_lisa", "schlund_manuel"],
                    "caption": "Boxplot.",
                    "plot_types": ["box"],
                },
                "pyplot_kwargs": {},
                "default_settings": {
                    "caption": None,
                    "fontsize": None,
                    "plot_kwargs": {},
                    "pyplot_kwargs": {},
                    "scatter_kwargs": {
                        "marker": "x",
                        "s": 200,
                        "linewidths": 2,
                        "color": "red",
                        "zorder": 3,
                    },
                    "var_order": None,
                },
            },
            "benchmarking_diurnal_cycle": {
                "function": partial(
                    self.create_1d_benchmarking_plot,
                    "benchmarking_diurnal_cycle",
                ),
                "coords": (["hour"],),
                "provenance": {
                    "authors": ["lauer_axel", "schlund_manuel"],
                    "caption": (
                        "Diurnal cycle of {long_name} for various datasets."
                    ),
                    "plot_types": ["seas"],
                },
                "pyplot_kwargs": {
                    "xlabel": "hour",
                    "xticks": {
                        "ticks": range(24),
                        "minor": True,
                    },
                },
                "default_settings": {**default_settings_1d},
            },
            "benchmarking_map": {
                "function": partial(
                    self.create_2d_benchmarking_plot,
                    "benchmarking_map",
                ),
                "coords": (["longitude", "latitude"],),
                "provenance": {
                    "authors": ["lauer_axel", "schlund_manuel"],
                    "caption": "Map plot of {long_name} of dataset {alias}.",
                    "plot_types": ["map"],
                },
                "pyplot_kwargs": {},
                "default_settings": {
                    **default_settings_2d,
                    "cbar_kwargs": {"orientation": "horizontal", "aspect": 30},
                    "gridline_kwargs": {},
                    "plot_kwargs": {"default": {"extend": "both"}},
                    "projection": "Robinson",
                    "projection_kwargs": {"central_longitude": 10},
                },
            },
            "benchmarking_timeseries": {
                "function": partial(
                    self.create_1d_benchmarking_plot,
                    "benchmarking_timeseries",
                ),
                "coords": (["time"],),
                "provenance": {
                    "authors": ["lauer_axel", "schlund_manuel"],
                    "caption": (
                        "Time series of {long_name} for various datasets."
                    ),
                    "plot_types": ["line"],
                },
                "pyplot_kwargs": {
                    "xlabel": "time",
                },
                "default_settings": {**default_settings_1d},
            },
            "benchmarking_zonal": {
                "function": partial(
                    self.create_2d_benchmarking_plot,
                    "benchmarking_zonal",
                ),
                "coords": (
                    ["latitude", "air_pressure"],
                    ["latitude", "altitude"],
                ),
                "provenance": {
                    "authors": ["lauer_axel", "schlund_manuel"],
                    "caption": (
                        "Zonal mean profile of {long_name} of dataset {alias}."
                    ),
                    "plot_types": ["vert"],
                },
                "pyplot_kwargs": {
                    "xlabel": "latitude [°N]",
                },
                "default_settings": {
                    **default_settings_2d,
                    "log_y": True,
                    "plot_kwargs": {"default": {"extend": "both"}},
                    "y_minor_formatter": "",
                },
            },
            "diurnal_cycle": {
                "function": partial(self.create_1d_plot, "diurnal_cycle"),
                "coords": (["hour"],),
                "provenance": {
                    "authors": ["lauer_axel", "schlund_manuel"],
                    "caption": (
                        "Diurnal cycle of {long_name} for various datasets."
                    ),
                    "plot_types": ["seas"],
                },
                "pyplot_kwargs": {
                    "xlabel": "hour",
                    "xticks": {
                        "ticks": range(24),
                        "minor": True,
                    },
                },
                "default_settings": {**default_settings_1d},
            },
            "hovmoeller_anncyc_vs_lat_or_lon": {
                "function": partial(
                    self.create_2d_plot,
                    "hovmoeller_anncyc_vs_lat_or_lon",
                ),
                "coords": (
                    ["month_number", "latitude"],
                    ["month_number", "longitude"],
                ),
                "provenance": {
                    "authors": ["hassler_birgit", "schlund_manuel"],
                    "caption": (
                        "Hovmoeller plot of {long_name} of dataset {alias}."
                    ),
                    "plot_types": ["zonal"],
                },
                "pyplot_kwargs": {
                    "xlabel": "month",
                    "xticks": {
                        "ticks": range(1, 13),
                        "labels": [str(m) for m in range(1, 13)],
                    },
                },
                "default_settings": {**default_settings_2d},
            },
            "hovmoeller_time_vs_lat_or_lon": {
                "function": partial(
                    self.create_2d_plot,
                    "hovmoeller_time_vs_lat_or_lon",
                ),
                "coords": (["time", "latitude"], ["time", "longitude"]),
                "provenance": {
                    "authors": [
                        "schlund_manuel",
                        "kraft_jeremy",
                        "lindenlaub_lukas",
                    ],
                    "caption": (
                        "Hovmoeller plot of {long_name} of dataset {alias}."
                    ),
                    "plot_types": ["zonal"],
                },
                "pyplot_kwargs": {
                    "xlabel": "time",
                },
                "default_settings": {**default_settings_2d},
            },
            "hovmoeller_z_vs_time": {
                "function": partial(
                    self.create_2d_plot,
                    "hovmoeller_z_vs_time",
                ),
                "coords": (["time", "air_pressure"], ["time", "altitude"]),
                "provenance": {
                    "authors": [
                        "kuehbacher_birgit",
                        "heuer_helge",
                        "schlund_manuel",
                    ],
                    "caption": (
                        "Hovmoeller plot of {long_name} of dataset {alias}."
                    ),
                    "plot_types": ["vert"],
                },
                "pyplot_kwargs": {
                    "xlabel": "time",
                },
                "default_settings": {
                    **default_settings_2d,
                    "log_y": True,
                    "y_minor_formatter": "",
                },
            },
            "map": {
                "function": partial(self.create_2d_plot, "map"),
                "coords": (["longitude", "latitude"],),
                "provenance": {
                    "authors": ["schlund_manuel"],
                    "caption": "Map plot of {long_name} of dataset {alias}.",
                    "plot_types": ["map"],
                },
                "pyplot_kwargs": {},
                "default_settings": {
                    **default_settings_2d,
                    "cbar_kwargs": {"orientation": "horizontal", "aspect": 30},
                    "gridline_kwargs": {},
                    "projection": "Robinson",
                    "projection_kwargs": {"central_longitude": 10},
                    "x_pos_stats_avg": 0.0,
                    "x_pos_stats_bias": 0.92,
                },
            },
            "timeseries": {
                "function": partial(self.create_1d_plot, "timeseries"),
                "coords": (["time"],),
                "provenance": {
                    "authors": ["schlund_manuel"],
                    "caption": (
                        "Time series of {long_name} for various datasets."
                    ),
                    "plot_types": ["line"],
                },
                "pyplot_kwargs": {
                    "xlabel": "time",
                },
                "default_settings": {**default_settings_1d},
            },
            "variable_vs_lat": {
                "function": partial(self.create_1d_plot, "variable_vs_lat"),
                "coords": (["latitude"],),
                "provenance": {
                    "authors": ["sarauer_ellen", "schlund_manuel"],
                    "caption": (
                        "{long_name} vs. latitude for various datasets."
                    ),
                    "plot_types": ["line"],
                },
                "pyplot_kwargs": {
                    "xlabel": "latitude [°N]",
                },
                "default_settings": {**default_settings_1d},
            },
            "zonal_mean_profile": {
                "function": partial(self.create_2d_plot, "zonal_mean_profile"),
                "coords": (
                    ["latitude", "air_pressure"],
                    ["latitude", "altitude"],
                ),
                "provenance": {
                    "authors": ["schlund_manuel"],
                    "caption": (
                        "Zonal mean profile of {long_name} of dataset {alias}."
                    ),
                    "plot_types": ["vert"],
                },
                "pyplot_kwargs": {
                    "xlabel": "latitude [°N]",
                },
                "default_settings": {
                    **default_settings_2d,
                    "log_y": True,
                    "y_minor_formatter": "",
                },
            },
        }

        return plot_settings

    def __init__(self, cfg: dict) -> None:
        """Initialize class member."""
        super().__init__(cfg)

        # Get default settings
        self.cfg = deepcopy(self.cfg)
        self.cfg.setdefault("facet_used_for_labels", "dataset")
        self.cfg.setdefault("figure_kwargs", {"constrained_layout": True})
        self.cfg.setdefault("group_variables_by", "short_name")
        self.cfg.setdefault("matplotlib_rc_params", {})
        self.cfg.setdefault(
            "savefig_kwargs",
            {
                "bbox_inches": "tight",
                "dpi": 300,
                "orientation": "landscape",
            },
        )
        self.cfg.setdefault("seaborn_settings", {"style": "ticks"})
        logger.info(
            "Using facet '%s' to group variables",
            self.cfg["group_variables_by"],
        )
        logger.info(
            "Using facet '%s' to create labels",
            self.cfg["facet_used_for_labels"],
        )

        # Load input data
        self.input_data = self._load_and_preprocess_data()
        self.grouped_input_data = group_metadata(
            self.input_data,
            self.cfg["group_variables_by"],
            sort=self.cfg["facet_used_for_labels"],
        )

        # Check given plot types and set default settings for them
        for plot_type, plot_options in self.plots.items():
            if plot_type not in self.plot_settings:
                raise ValueError(
                    f"Got unexpected plot type '{plot_type}' for option "
                    f"'plots', expected one of {list(self.plot_settings)}",
                )
            if plot_options is None:
                plot_options = {}
                self.plots[plot_type] = plot_options

            # Only use default projection options if no projection is specified
            if "projection" in plot_options:
                self.plots[plot_type].setdefault("projection_kwargs", {})

            default_settings = self.plot_settings[plot_type][
                "default_settings"
            ]
            for key, val in default_settings.items():
                self.plots[plot_type].setdefault(key, val)

        # Check that facet_used_for_labels is present for every dataset
        for dataset in self.input_data:
            if self.cfg["facet_used_for_labels"] not in dataset:
                raise ValueError(
                    f"facet_used_for_labels "
                    f"'{self.cfg['facet_used_for_labels']}' not present for "
                    f"the following dataset:\n{pformat(dataset)}",
                )

        # Load seaborn settings
        sns.set_theme(**self.cfg["seaborn_settings"])

    def _add_colorbar(
        self,
        plot_type: str,
        plot_1: Any,
        axes_1: Axes,
        dataset_1: dict,
        plot_2: Any | None = None,
        axes_2: Axes | None = None,
        dataset_2: dict | None = None,
        *,
        bias: bool = False,
    ) -> None:
        """Add colorbar(s) for plots."""
        fontsize = (
            self.plots[plot_type]["fontsize"] or mpl.rcParams["axes.labelsize"]
        )
        cbar_kwargs = self._get_cbar_kwargs(plot_type, bias=bias)

        def add_colorbar(
            plot: Any,
            axes: Axes | Iterable[Axes],
            label: str,
        ) -> None:
            """Add single colorbar."""
            cbar = plt.colorbar(plot, ax=axes, **cbar_kwargs)
            cbar.set_label(label, fontsize=fontsize)
            cbar.ax.tick_params(labelsize=fontsize)

        # Single colorbar
        cbar_label_1 = self._get_cbar_label(plot_type, dataset_1, bias=bias)
        if plot_2 is None or axes_2 is None or dataset_2 is None:
            add_colorbar(plot_1, axes_1, cbar_label_1)
            return

        # One common colorbar
        # Note: Increase aspect ratio for nicer looks
        if self.plots[plot_type]["common_cbar"]:
            if "aspect" in cbar_kwargs:
                cbar_kwargs["aspect"] += 20.0
            add_colorbar(plot_1, [axes_1, axes_2], cbar_label_1)
            return

        # Two separate colorbars
        add_colorbar(plot_1, axes_1, cbar_label_1)
        cbar_label_2 = self._get_cbar_label(plot_type, dataset_2)
        add_colorbar(plot_2, axes_2, cbar_label_2)

    def _add_stats(
        self,
        plot_type: str,
        axes: Axes,
        dataset_1: dict,
        dataset_2: dict | None = None,
    ) -> None:
        """Add text to plot that describes basic statistics."""
        cube_1 = dataset_1["cube"]
        if dataset_2 is None:
            cube_2 = None
            label = self._get_label(dataset_1)
        else:
            cube_2 = dataset_2["cube"]
            label = (
                f"{self._get_label(dataset_1)} vs. "
                f"{self._get_label(dataset_2)}"
            )
        coords = [c.name() for c in cube_1.coords(dim_coords=True)]

        # Different options for the different plots types
        fontsize = 6.0
        y_pos = 0.95
        if all(
            [
                "x_pos_stats_avg" in self.plots[plot_type],
                "x_pos_stats_bias" in self.plots[plot_type],
            ],
        ):
            x_pos_bias = self.plots[plot_type]["x_pos_stats_bias"]
            x_pos = self.plots[plot_type]["x_pos_stats_avg"]
        else:
            raise NotImplementedError(f"plot_type '{plot_type}' not supported")

        # Mean
        weights = area_weights(cube_1)
        if cube_2 is None:
            mean = cube_1.collapsed(
                coords,
                iris.analysis.MEAN,
                weights=weights,
            )
            logger.info(
                "Area-weighted mean of %s for %s = %f%s",
                dataset_1["short_name"],
                label,
                mean.data,
                dataset_1["units"],
            )
        else:
            mean = (cube_1 - cube_2).collapsed(
                coords,
                iris.analysis.MEAN,
                weights=weights,
            )
            logger.info(
                "Area-weighted bias of %s for %s = %f%s",
                dataset_1["short_name"],
                label,
                mean.data,
                dataset_1["units"],
            )
        if np.abs(mean.data) >= 0.1:
            mean_val = f"{mean.data:.2f} {cube_1.units}"
        else:
            mean_val = f"{mean.data:.2e} {cube_1.units}"
        axes.text(
            x_pos,
            y_pos,
            mean_val,
            fontsize=fontsize,
            transform=axes.transAxes,
        )

        if cube_2 is None:
            return

        # Weighted RMSE
        rmse = (cube_1 - cube_2).collapsed(
            coords,
            iris.analysis.RMS,
            weights=weights,
        )
        if np.abs(rmse.data) >= 0.1:
            rmse_val = f"{rmse.data:.2f} {cube_1.units}"
        else:
            rmse_val = f"{rmse.data:.2e} {cube_1.units}"
        axes.text(
            x_pos_bias,
            y_pos,
            f"RMSE={rmse_val}",
            fontsize=fontsize,
            transform=axes.transAxes,
        )
        logger.info(
            "Area-weighted RMSE of %s for %s = %f%s",
            dataset_1["short_name"],
            label,
            rmse.data,
            dataset_1["units"],
        )

        # Weighted R2
        mask = np.ma.getmaskarray(cube_1.data).ravel()
        mask |= np.ma.getmaskarray(cube_2.data).ravel()
        cube_data = cube_1.data.ravel()[~mask]
        ref_cube_data = cube_2.data.ravel()[~mask]
        weights = weights.ravel()[~mask]
        r2_val = r2_score(cube_data, ref_cube_data, sample_weight=weights)
        axes.text(
            x_pos_bias,
            y_pos - 0.1,
            rf"R$^2$={r2_val:.2f}",
            fontsize=fontsize,
            transform=axes.transAxes,
        )
        logger.info(
            "Area-weighted R2 of %s for %s = %f",
            dataset_1["short_name"],
            label,
            r2_val,
        )

    def _check_cube_coords(self, cube: Cube, plot_type: str) -> list[str]:
        """Check that cube has correct dimensional coordinates."""
        expected_dimensions = self.plot_settings[plot_type]["coords"]
        for dims in expected_dimensions:
            cube_dims = [cube.coords(dim, dim_coords=True) for dim in dims]
            if all(cube_dims) and cube.ndim == len(dims):
                return dims
        expected_dims_str = " or ".join(
            [str(dims) for dims in expected_dimensions],
        )
        raise ValueError(
            f"Expected cube with dimensional coordinates "
            f"{expected_dims_str}, got {cube.summary(shorten=True)}",
        )

    def _customize_plot(
        self,
        plot_type: str,
        axes: Axes,
        dataset: dict,
    ) -> Axes:
        """Customize plot with user-defined settings."""
        self._process_pyplot_kwargs(
            self.plot_settings[plot_type]["pyplot_kwargs"],
            dataset,
            transpose_axes=self.plots[plot_type]["transpose_axes"],
        )

        # Aspect ratio
        if self.plots[plot_type]["aspect_ratio"] is not None:
            axes.set_box_aspect(self.plots[plot_type]["aspect_ratio"])

        # Axes styles
        if self.plots[plot_type]["log_x"]:
            axes.set_xscale("log")
            axes.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
            x_minor_locator = LogLocator(
                base=10.0,
                subs=np.arange(1.0, 10.0) * 0.1,  # type: ignore
                numticks=12,
            )
        else:
            x_minor_locator = AutoMinorLocator()

        if self.plots[plot_type]["log_y"]:
            axes.set_yscale("log")
            axes.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
            y_minor_locator = LogLocator(
                base=10.0,
                subs=np.arange(1.0, 10.0) * 0.1,
                numticks=12,
            )
        else:
            y_minor_locator = AutoMinorLocator()

        if self.plots[plot_type]["x_major_formatter"] is not None:
            axes.xaxis.set_major_formatter(
                FormatStrFormatter(self.plots[plot_type]["x_major_formatter"]),
            )
        if self.plots[plot_type]["x_minor_formatter"] is not None:
            axes.xaxis.set_minor_locator(x_minor_locator)
            axes.xaxis.set_minor_formatter(
                FormatStrFormatter(self.plots[plot_type]["x_minor_formatter"]),
            )
        else:
            axes.xaxis.set_minor_locator(NullLocator())

        if self.plots[plot_type]["y_major_formatter"] is not None:
            axes.yaxis.set_major_formatter(
                FormatStrFormatter(self.plots[plot_type]["y_major_formatter"]),
            )
        if self.plots[plot_type]["y_minor_formatter"] is not None:
            axes.yaxis.set_minor_locator(y_minor_locator)
            axes.yaxis.set_minor_formatter(
                FormatStrFormatter(self.plots[plot_type]["y_minor_formatter"]),
            )
        else:
            axes.yaxis.set_minor_locator(NullLocator())

        if self.plots[plot_type]["time_format"] is not None:
            if self.plots[plot_type]["transpose_axes"]:
                time_axis = axes.yaxis
            else:
                time_axis = axes.xaxis
            time_axis.set_major_formatter(
                mdates.DateFormatter(self.plots[plot_type]["time_format"]),
            )

        # Gridlines
        gridline_kwargs = self._get_gridline_kwargs(plot_type)
        if gridline_kwargs is not False:
            if "map" in plot_type:
                axes.gridlines(**gridline_kwargs)
            else:
                axes.grid(**gridline_kwargs)

        # Legend
        legend_kwargs = self.plots[plot_type]["legend_kwargs"]
        if legend_kwargs is not False:
            axes.legend(**legend_kwargs)

        # Rasterization
        if self.plots[plot_type]["rasterize"]:
            self._set_rasterized([axes])

        # Further customize plot appearance
        self._process_pyplot_kwargs(
            self.plots[plot_type]["pyplot_kwargs"],
            dataset,
        )

        return axes

    @staticmethod
    def _fill_facet_placeholders(
        string: str,
        dataset: dict,
        description: str,
    ) -> str:
        """Fill facet placeholders."""
        try:
            string = string.format(**dataset)
        except KeyError as exc:
            raise ValueError(
                f"Not all necessary facets in {description} available for "
                f"dataset\n{pformat(dataset)}",
            ) from exc
        return string

    def _get_benchmark_datasets(self, datasets: list[dict]) -> list[dict]:
        """Get dataset to be benchmarked."""
        variable = datasets[0][self.cfg["group_variables_by"]]
        benchmark_datasets = [
            d for d in datasets if d.get("benchmark_dataset", False)
        ]
        if len(benchmark_datasets) >= 1:
            return benchmark_datasets

        raise ValueError(
            f"Expected at least 1 benchmark dataset (with 'benchmark_dataset: "
            f"True' for variable '{variable}'), got "
            f"{len(benchmark_datasets):d}",
        )

    def _get_benchmark_mask(
        self,
        benchmark_dataset: dict,
        percentile_datasets: list[dict],
    ) -> Cube:
        """Create mask for benchmarking cube depending on metric."""
        metric = self._get_benchmark_metric(benchmark_dataset)
        cube = benchmark_dataset["cube"]
        percentile_cubes = [d["cube"] for d in percentile_datasets]

        if metric == "bias":
            maxabs_perc = np.maximum(
                np.abs(percentile_cubes[0].data),  # largest percentile
                np.abs(percentile_cubes[-1].data),  # smallest percentile
            )
            mask = np.where(np.abs(cube.data) >= maxabs_perc, 0, 1)
        elif metric in ("emd", "rmse"):
            mask = np.where(cube.data >= percentile_cubes[0].data, 0, 1)
        elif metric == "pearsonr":
            mask = np.where(cube.data <= percentile_cubes[0].data, 0, 1)
        else:
            raise ValueError(
                f"Could not create benchmarking mask, unknown benchmarking "
                f"metric: '{metric}'",
            )

        return cube.copy(mask)

    def _get_benchmark_metric(self, dataset: dict) -> str:
        """Get benchmarking metric."""
        for metric in ("emd", "pearsonr", "rmse"):
            if dataset["short_name"].startswith(f"{metric}_"):
                return metric
        metric = "bias"
        logger.info(
            "Could not determine metric from short_name, assuming "
            "benchmarking metric = %s",
            metric,
        )
        return metric

    def _get_benchmark_percentiles(self, datasets: list[dict]) -> list[dict]:
        """Get percentile datasets from multi-model statistics preprocessor."""
        percentile_datasets = []
        for dataset in datasets:
            stat = dataset.get("multi_model_statistics")
            if stat is not None and "MultiModelPercentile" in stat:
                dataset["_percentile_int"] = int(
                    stat.replace("MultiModelPercentile", ""),
                )
                percentile_datasets.append(dataset)

        # Sort percentiles in descending order (highest percentile first)
        percentile_datasets = list(
            reversed(sorted_metadata(percentile_datasets, "_percentile_int")),
        )

        # Get number of percentiles expected depending on benchmarking metric
        metric = self._get_benchmark_metric(datasets[0])
        n_percentiles: dict[str, int] = {
            "bias": 2,
            "rmse": 1,
            "pearsonr": 1,
            "emd": 1,
        }
        if metric not in n_percentiles:
            raise ValueError(f"Unknown benchmarking metric: '{metric}'.")

        if len(percentile_datasets) >= n_percentiles[metric]:
            return percentile_datasets

        variable = datasets[0][self.cfg["group_variables_by"]]
        raise ValueError(
            f"Expected at least {n_percentiles[metric]} percentile datasets "
            f"(created with multi-model statistics preprocessor for variable "
            f"'{variable}'), got {len(percentile_datasets):d}",
        )

    def _get_bias_dataset(self, dataset_1: dict, dataset_2: dict) -> dict:
        """Get bias dataset (dataset_1 - dataset_2)."""
        bias_cube = dataset_1["cube"] - dataset_2["cube"]
        bias_cube.metadata = dataset_1["cube"].metadata
        bias_cube.var_name = (
            None
            if bias_cube.var_name is None
            else f"bias_{bias_cube.var_name}"
        )
        bias_cube.long_name = (
            None
            if bias_cube.long_name is None
            else f"Bias in {bias_cube.long_name}"
        )

        dataset_bias = deepcopy(dataset_1)
        dataset_bias["cube"] = bias_cube
        dataset_bias[self.cfg["facet_used_for_labels"]] = (
            f"{self._get_label(dataset_1)} - {self._get_label(dataset_2)}"
        )
        dataset_bias["ancestors"] = [
            dataset_1["filename"],
            dataset_2["filename"],
        ]

        return dataset_bias

    def _get_cbar_kwargs(self, plot_type: str, *, bias: bool = False) -> dict:
        """Get colorbar kwargs."""
        cbar_kwargs = deepcopy(self.plots[plot_type]["cbar_kwargs"])
        if bias:
            cbar_kwargs.update(self.plots[plot_type]["cbar_kwargs_bias"])
        return deepcopy(cbar_kwargs)

    def _get_cbar_label(
        self,
        plot_type: str,
        dataset: dict,
        *,
        bias: bool = False,
    ) -> str:
        """Get colorbar label."""
        if bias:
            cbar_label = self.plots[plot_type]["cbar_label_bias"]
            descr = f"cbar_label_bias of {plot_type} '{cbar_label}'"
        else:
            cbar_label = self.plots[plot_type]["cbar_label"]
            descr = f"cbar_label of {plot_type} '{cbar_label}'"
        cbar_label = self._fill_facet_placeholders(cbar_label, dataset, descr)
        return cbar_label

    def _get_coords_for_2d_plotting(
        self,
        plot_type: str,
        cube: Cube,
    ) -> tuple[Coord, Coord]:
        """Get coordinates to plot 2D data."""
        coords = self._check_cube_coords(cube, plot_type)
        if self.plots[plot_type]["transpose_axes"]:
            x_coord = cube.coord(coords[1], dim_coords=True)
            y_coord = cube.coord(coords[0], dim_coords=True)
        else:
            x_coord = cube.coord(coords[0], dim_coords=True)
            y_coord = cube.coord(coords[1], dim_coords=True)
        return (x_coord, y_coord)

    def _get_custom_mpl_rc_params(self, plot_type: str) -> mpl.RcParams:
        """Get custom matplotlib rcParams."""
        custom_rc_params = copy(self.cfg["matplotlib_rc_params"])
        fontsize = self.plots[plot_type].get("fontsize")
        if fontsize is not None:
            custom_rc_params.update(
                {
                    "axes.titlesize": fontsize + 2.0,
                    "axes.labelsize": fontsize,
                    "xtick.labelsize": fontsize,
                    "ytick.labelsize": fontsize,
                },
            )
        return custom_rc_params

    def _get_gridline_kwargs(self, plot_type: str) -> dict:
        """Get gridline kwargs."""
        gridline_kwargs = self.plots[plot_type]["gridline_kwargs"]
        return deepcopy(gridline_kwargs)

    def _get_label(self, dataset: dict) -> str:
        """Get label of dataset."""
        return dataset[self.cfg["facet_used_for_labels"]]

    @staticmethod
    def _get_multi_dataset_facets(datasets: list[dict]) -> dict:
        """Derive common facets for multiple datasets."""
        all_keys = {key for dataset in datasets for key in dataset}
        multi_dataset_facets = {}
        for key in all_keys:
            if all(d.get(key) == datasets[0].get(key) for d in datasets):
                multi_dataset_facets[key] = datasets[0].get(key)
            else:
                multi_dataset_facets[key] = f"ambiguous_{key}"
        return multi_dataset_facets

    def _get_netcdf_path(
        self,
        plot_path: Path | str,
        suffix: str | None = None,
    ) -> str:
        """Get netCDF path."""
        basename = Path(plot_path).stem
        if suffix is not None:
            basename += suffix
        return get_diagnostic_filename(basename, self.cfg)

    def _get_plot_func(self, plot_type: str) -> Callable:
        """Get plot function."""
        plot_func = self.plots[plot_type]["plot_func"]
        if not hasattr(iris.plot, plot_func):
            raise AttributeError(
                f"Got invalid plot function '{plot_func}' for plotting "
                f"{plot_type}, expected function of iris.plot",
            )
        logger.info(
            "Creating %s plots using function '%s'",
            plot_type,
            plot_func,
        )
        return getattr(iris.plot, plot_func)

    def _get_plot_kwargs(
        self,
        plot_type: str,
        dataset: dict,
        *,
        bias: bool = False,
    ) -> dict:
        """Get keyword arguments for plot functions."""
        all_plot_kwargs = self.plots[plot_type]["plot_kwargs"]
        all_plot_kwargs = deepcopy(all_plot_kwargs)

        # First get default kwargs, then overwrite them with dataset-specific
        # ones
        plot_kwargs = all_plot_kwargs.get("default", {})
        label = self._get_label(dataset)
        plot_kwargs.update(all_plot_kwargs.get(label, {}))

        # For bias plots, overwrite the kwargs with bias-specific option
        if bias:
            bias_kwargs = self.plots[plot_type]["plot_kwargs_bias"]
            bias_kwargs.setdefault("cmap", "bwr")
            bias_kwargs.setdefault("norm", "centered")
            plot_kwargs.update(bias_kwargs)

        # Replace facets with dataset entries for string arguments
        for key, val in plot_kwargs.items():
            if isinstance(val, str):
                val = self._fill_facet_placeholders(
                    val,
                    dataset,
                    f"plot_kwargs of {plot_type} '{key}: {val}'",
                )
                plot_kwargs[key] = val

        # Handle special plot kwargs
        if plot_kwargs.get("norm") == "centered":
            norm = CenteredNorm(
                vcenter=plot_kwargs.pop("vcenter", 0.0),
                halfrange=plot_kwargs.pop("halfrange", None),
            )
            plot_kwargs["norm"] = norm

        return deepcopy(plot_kwargs)

    def _get_projection(self, plot_type: str) -> Any:
        """Get plot projection."""
        projection = self.plots[plot_type]["projection"]
        projection_kwargs = self.plots[plot_type]["projection_kwargs"]

        # Check if desired projection is valid
        if not hasattr(ccrs, projection):
            raise AttributeError(
                f"Got invalid projection '{projection}' for plotting "
                f"{plot_type}, expected class of cartopy.crs",
            )

        return getattr(ccrs, projection)(**projection_kwargs)

    def _get_provenance_record(
        self,
        plot_type: str,
        representative_dataset: dict,
        datasets: list[dict],
    ) -> dict:
        """Get provenance record."""
        provenance_record: dict[str, str | Iterable[str]] = {
            "ancestors": list({a for d in datasets for a in d["ancestors"]}),
            "long_names": list({d["long_name"] for d in datasets}),
        }
        provenance_record.update(self.plot_settings[plot_type]["provenance"])
        if self.plots[plot_type]["caption"] is not None:
            provenance_record["caption"] = self.plots[plot_type]["caption"]
        for prov_key, prov_val in provenance_record.items():
            if isinstance(prov_val, str):
                provenance_record[prov_key] = self._fill_facet_placeholders(
                    prov_val,
                    representative_dataset,
                    f"provenance {prov_key} of {plot_type} '{prov_val}'",
                )
        return provenance_record

    def _get_reference_dataset(self, datasets: list[dict]) -> dict | None:
        """Extract reference dataset."""
        variable = datasets[0][self.cfg["group_variables_by"]]
        ref_datasets = [
            d for d in datasets if d.get("reference_for_monitor_diags", False)
        ]
        if len(ref_datasets) > 1:
            raise ValueError(
                f"Expected at most 1 reference dataset (with "
                f"'reference_for_monitor_diags: True' for variable "
                f"'{variable}', got {len(ref_datasets):d}",
            )
        if ref_datasets:
            return ref_datasets[0]
        return None

    def _load_and_preprocess_data(self) -> list[dict]:
        """Load and preprocess data."""
        input_data = list(self.cfg["input_data"].values())

        if not input_data:
            raise ValueError("No input data given")

        slices = not any(
            self.cfg["group_variables_by"] in ds for ds in input_data
        )
        datasets = []
        for dataset in input_data:
            filename = dataset["filename"]
            logger.info("Loading %s", filename)
            cubes = iris.load(filename)
            if len(cubes) == 1:
                cube: Cube = cubes[0]  # type: ignore
            else:
                var_name = dataset["short_name"]
                try:
                    cube = cubes.extract_cube(
                        iris.NameConstraint(var_name=var_name),
                    )  # type: ignore
                except ConstraintMismatchError as exc:
                    var_names = [c.var_name for c in cubes]
                    raise ValueError(
                        f"Cannot load data: multiple variables ({var_names}) "
                        f"are available in file {filename}, but not the "
                        f"requested '{var_name}'",
                    ) from exc

            # Fix time coordinate if present
            if cube.coords("time", dim_coords=True):
                ih.unify_time_coord(cube)

            # Add scalar latitude and longitude coordinates if these are not
            # present (necessary for calculation of area weights). The exact
            # values for the points/bounds of these coordinates do not matter
            # since they don't change the weights.
            if not cube.coords("latitude"):
                lon_coord = AuxCoord(
                    0.0,
                    bounds=[-90.0, 90.0],
                    var_name="lat",
                    standard_name="latitude",
                    long_name="latitude",
                    units="degrees_north",
                )
                cube.add_aux_coord(lon_coord, ())
            if not cube.coords("longitude"):
                lon_coord = AuxCoord(
                    180.0,
                    bounds=[0.0, 360.0],
                    var_name="lon",
                    standard_name="longitude",
                    long_name="longitude",
                    units="degrees_east",
                )
                cube.add_aux_coord(lon_coord, ())

            # Fix Z-coordinate if present
            if cube.coords("air_pressure", dim_coords=True):
                z_coord = cube.coord("air_pressure", dim_coords=True)
                z_coord.attributes["positive"] = "down"
                z_coord.convert_units("hPa")
            elif cube.coords("altitude", dim_coords=True):
                z_coord = cube.coord("altitude")
                z_coord.attributes["positive"] = "up"

            # Save ancestors
            dataset["ancestors"] = [filename]

            if slices:
                slice_coord_name = self.cfg["group_variables_by"]
                for subcube in cube.slices_over([slice_coord_name]):
                    dataset_copy = deepcopy(dataset)
                    dataset_copy["cube"] = subcube
                    dataset_copy[slice_coord_name] = subcube.coord(
                        slice_coord_name,
                    ).points[0]
                    datasets.append(dataset_copy)
            else:
                dataset_copy = deepcopy(dataset)
                dataset_copy["cube"] = cube
                datasets.append(dataset_copy)
        return datasets

    def _plot_1d_data(
        self,
        plot_type: str,
        datasets: list[dict],
        axes: Axes,
    ) -> None:
        """Plot 1D data."""
        # Plot all datasets in one single figure
        coord_label = "unkown coordinate"
        for dataset in datasets:
            label = self._get_label(dataset)
            cube = dataset["cube"]
            coords = self._check_cube_coords(cube, plot_type)
            coord = cube.coord(coords[0], dim_coords=True)
            coord_label = f"{coord.name()} [{coord.units}]"

            # Actual plot
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs.setdefault("label", label)
            plot_kwargs["axes"] = axes
            if self.plots[plot_type]["transpose_axes"]:
                iris.plot.plot(cube, coord, **plot_kwargs)
            else:
                iris.plot.plot(coord, cube, **plot_kwargs)

        # Plot horizontal lines
        for hline_kwargs in self.plots[plot_type]["hlines"]:
            axes.axhline(**hline_kwargs)

        # Title and axis labels
        multi_dataset_facets = self._get_multi_dataset_facets(datasets)
        axes.set_title(multi_dataset_facets["long_name"])
        var_label = (
            f"{multi_dataset_facets[self.cfg['group_variables_by']]} "
            f"[{multi_dataset_facets['units']}]"
        )
        if self.plots[plot_type]["transpose_axes"]:
            axes.set_xlabel(var_label)
            axes.set_ylabel(coord_label)
        else:
            axes.set_xlabel(coord_label)
            axes.set_ylabel(var_label)

        # Customize plot with user-defined settings
        self._customize_plot(plot_type, axes, multi_dataset_facets)

    def _plot_2d(self, plot_type: str, cube: Cube, **plot_kwargs: Any) -> Any:
        """Plot 2D data (plain plotting, no changes in plot appearance)."""
        plot_func = self._get_plot_func(plot_type)

        # Setup coordinates
        coords = self._get_coords_for_2d_plotting(plot_type, cube)
        plot_kwargs["coords"] = coords

        # Fix Cartopy bug if necessary (see
        # https://github.com/SciTools/cartopy/issues/2457 and
        # https://github.com/SciTools/cartopy/issues/2468)
        fix_cartopy_bug = all(
            [
                self.plots[plot_type]["projection"] == "Robinson",
                plot_func is iris.plot.contourf,
            ],
        )
        if fix_cartopy_bug:
            plot_kwargs["transform_first"] = True
            npx = da if cube.has_lazy_data() else np
            cube = cube.copy(npx.ma.filled(cube.core_data(), np.nan))  # type: ignore

        return plot_func(cube, **plot_kwargs)

    def _plot_2d_data(
        self,
        plot_type: str,
        dataset: dict,
        axes: Axes,
        *,
        bias: bool = False,
        **additional_plot_kwargs: Any,
    ) -> Any:
        """Plot 2D data."""
        fig: Figure = axes.get_figure()  # type: ignore

        # Some options are not supported for map plots
        if "map" in plot_type:
            self.plots[plot_type]["aspect_ratio"] = None

        # Plot data
        cube = dataset["cube"]
        plot_kwargs = self._get_plot_kwargs(plot_type, dataset, bias=bias)
        plot_kwargs.update(additional_plot_kwargs)
        plot_kwargs["axes"] = axes
        plot_output = self._plot_2d(plot_type, cube, **plot_kwargs)

        # Show coastlines for map plots
        if "map" in plot_type:
            axes.coastlines()  # type: ignore

        # Title and axis labels
        fig.suptitle(dataset["long_name"])
        axes.set_title(self._get_label(dataset))
        (x_coord, y_coord) = self._get_coords_for_2d_plotting(plot_type, cube)
        axes.set_xlabel(f"{x_coord.name()} [{x_coord.units}]")  # type: ignore
        axes.set_ylabel(f"{y_coord.name()} [{y_coord.units}]")  # type: ignore

        # Customize plot with user-defined settings
        self._customize_plot(plot_type, axes, dataset)

        return plot_output

    def _plot_2d_data_1_panel(
        self,
        plot_type: str,
        dataset: dict,
    ) -> tuple[Figure, Axes]:
        """Plot 2D data (single panel)."""
        fig = plt.figure(**self.cfg["figure_kwargs"])
        subplot_kwargs = {}
        if self.plots[plot_type]["projection"] is not None:
            subplot_kwargs["projection"] = self._get_projection(plot_type)
        axes = fig.add_subplot(**subplot_kwargs)

        # Plot data
        plot_output = self._plot_2d_data(plot_type, dataset, axes)
        self._add_colorbar(plot_type, plot_output, axes, dataset)

        # Show statistics if desired
        if self.plots[plot_type]["show_stats"]:
            self._add_stats(plot_type, axes, dataset)

        return fig, axes

    def _plot_2d_data_3_panel(
        self,
        plot_type: str,
        dataset_left: dict,
        dataset_right: dict,
    ) -> Figure:
        """Plot 2D data (three panels)."""
        # Create single figure with multiple axes
        fig = plt.figure(**self.cfg["figure_kwargs"])
        gridspec = GridSpec(
            5,
            4,
            figure=fig,
            height_ratios=[1.0, 1.0, 0.4, 1.0, 1.0],
        )
        subplot_kwargs = {}
        if self.plots[plot_type]["projection"] is not None:
            subplot_kwargs["projection"] = self._get_projection(plot_type)

        # Plot top left panel
        axes_left = fig.add_subplot(gridspec[0:2, 0:2], **subplot_kwargs)
        plot_left = self._plot_2d_data(plot_type, dataset_left, axes_left)

        # Plot top right panel
        # Note: make sure to use the same vmin and vmax than the top left plot
        # if a common colorbar is desired
        axes_right = fig.add_subplot(
            gridspec[0:2, 2:4],
            sharex=axes_left,
            sharey=axes_left,
            **subplot_kwargs,
        )
        additional_plot_kwargs = {}
        if self.plots[plot_type]["common_cbar"]:
            additional_plot_kwargs.setdefault("vmin", plot_left.get_clim()[0])
            additional_plot_kwargs.setdefault("vmax", plot_left.get_clim()[1])
        plot_right = self._plot_2d_data(
            plot_type,
            dataset_right,
            axes_right,
            **additional_plot_kwargs,
        )

        # Colorbar(s) for top plots
        self._add_colorbar(
            plot_type,
            plot_left,
            axes_left,
            dataset_left,
            plot_right,
            axes_right,
            dataset_right,
        )

        # Plot bottom panel (bias) including colorbar
        axes_bottom = fig.add_subplot(
            gridspec[3:5, 1:3],
            sharex=axes_left,
            sharey=axes_left,
            **subplot_kwargs,
        )
        dataset_bottom = self._get_bias_dataset(dataset_left, dataset_right)
        plot_bottom = self._plot_2d_data(
            plot_type,
            dataset_bottom,
            axes_bottom,
            bias=True,
        )
        self._add_colorbar(
            plot_type,
            plot_bottom,
            axes_bottom,
            dataset_bottom,
            bias=True,
        )

        # Show statistics if desired
        if self.plots[plot_type]["show_stats"]:
            self._add_stats(plot_type, axes_left, dataset_left)
            self._add_stats(plot_type, axes_right, dataset_right)
            self._add_stats(
                plot_type,
                axes_bottom,
                dataset_left,
                dataset_right,
            )

        # Hide superfluous axis labels
        plt.setp(axes_right.get_yticklabels(which="both"), visible=False)
        axes_left.set_xlabel("")
        axes_right.set_xlabel("")
        axes_right.set_ylabel("")

        return fig

    def _plot_benchmarking_boxplot(
        self,
        dframe: pd.DataFrame,
        benchmark_datasets: dict[str, dict],
    ) -> None:
        """Plot benchmarking boxplot."""
        plot_type = "benchmarking_boxplot"
        plot_kwargs = dict(self.plots[plot_type]["plot_kwargs"])
        scatter_kwargs = dict(self.plots[plot_type]["scatter_kwargs"])

        # Create plot with desired settings
        fig = plt.figure(**self.cfg["figure_kwargs"])
        fig.suptitle("Benchmarking Boxplot")

        for idx, var_key in enumerate(benchmark_datasets):
            benchmark_dataset = benchmark_datasets[var_key]
            cube_to_benchmark = benchmark_dataset["cube"]

            # Create separate plot for each variable
            axes = plt.subplot(1, len(benchmark_datasets), idx + 1)
            plot_kwargs["ax"] = axes
            df_to_plot = dframe[
                (dframe["Variable"] == var_key)
                & (dframe["Benchmark Dataset"] == False)  # noqa: E712
            ].drop("Benchmark Dataset", axis=1)
            plot_boxplot = sns.boxplot(df_to_plot, **plot_kwargs)
            plot_boxplot.set(xticklabels=[])

            # Plot data to benchmark
            plt.scatter(0, cube_to_benchmark.data, **scatter_kwargs)

            # Plot appearance
            plt.xlabel(var_key)
            if cube_to_benchmark.units != 1:
                plt.ylabel(cube_to_benchmark.units)
            self._process_pyplot_kwargs(
                self.plots[plot_type]["pyplot_kwargs"],
                benchmark_dataset,
            )

    def _process_pyplot_kwargs(
        self,
        pyplot_kwargs: dict[str, Any],
        dataset: dict,
        *,
        transpose_axes: bool = False,
    ) -> None:
        """Process functions for :mod:`matplotlib.pyplot`."""
        for func, arg in pyplot_kwargs.items():
            if transpose_axes and func.startswith("x"):
                func = func.replace("x", "y", 1)
            elif transpose_axes and func.startswith("y"):
                func = func.replace("y", "x", 1)
            if isinstance(arg, str):
                arg = self._fill_facet_placeholders(
                    arg,
                    dataset,
                    f"pyplot_kwargs '{func}: {arg}'",
                )
            if arg is None:
                getattr(plt, func)()
            elif isinstance(arg, dict):
                getattr(plt, func)(**arg)
            else:
                getattr(plt, func)(arg)

    def _save_1d_data(
        self,
        plot_type: str,
        datasets: list[dict],
        fig: Figure,
    ) -> None:
        """Save 1D plot and netCDF files."""
        multi_dataset_facets = self._get_multi_dataset_facets(datasets)

        # Save plot file
        plot_path = self._save_plot(plot_type, multi_dataset_facets, fig)

        # Save netCDF file
        cubes: dict[str, Cube] = {
            self._get_label(d): d["cube"] for d in datasets
        }
        netcdf_path = self._get_netcdf_path(plot_path)
        coord_name = datasets[0]["cube"].coord(dim_coords=True).name()
        var_attrs = {
            n: datasets[0][n] for n in ("short_name", "long_name", "units")
        }
        io.save_1d_data(cubes, netcdf_path, coord_name, var_attrs)

        # Provenance tracking
        provenance_record = self._get_provenance_record(
            plot_type,
            multi_dataset_facets,
            datasets,
        )
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)
            provenance_logger.log(netcdf_path, provenance_record)

    def _save_boxplot_data(
        self,
        plot_type: str,
        dframe: pd.DataFrame,
        benchmark_datasets: dict[str, dict],
    ) -> None:
        """Save plot and netCDF files for boxplot."""
        # Save plot
        all_vars = "_".join(benchmark_datasets)
        all_datasets = [d for g in self.grouped_input_data.values() for d in g]
        multi_dataset_facets = self._get_multi_dataset_facets(all_datasets)
        multi_dataset_facets["variable_group"] = all_vars
        multi_dataset_facets["short_name"] = all_vars
        plot_path = self.get_plot_path(plot_type, multi_dataset_facets)
        plt.savefig(plot_path, **self.cfg["savefig_kwargs"])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save netCDF file
        dframe = dframe.drop("Benchmark Dataset", axis=1)
        cubes_to_save = CubeList()
        for var_key, dataset in benchmark_datasets.items():
            df_single_var = dframe[dframe["Variable"] == var_key].drop(
                "Variable",
                axis=1,
            )
            cube: Cube = iris.pandas.as_cubes(
                df_single_var,
                aux_coord_cols=["Dataset"],
            )[0]  # type: ignore
            cube.var_name = var_key
            cube.long_name = dataset["long_name"]
            if dataset["standard_name"]:
                cube.standard_name = dataset["standard_name"]
                cube.units = dataset["units"]
            cube.remove_coord("unknown")

            # Avoid invalid dtype for dataset coordinate when saving netCDF
            dataset_list = cube.coord("Dataset").points.astype(str)
            cube.remove_coord("Dataset")
            cube.add_aux_coord(AuxCoord(dataset_list, var_name="dataset"), 0)

            cube.attributes["benchmark_dataset"] = self._get_label(dataset)
            cubes_to_save.append(cube)
        netcdf_path = self._get_netcdf_path(plot_path)
        io.iris_save(cubes_to_save, netcdf_path)

        # Provenance tracking
        ancestors = [
            d["filename"] for g in self.grouped_input_data.values() for d in g
        ]
        long_names = sorted(group_metadata(all_datasets, "long_name"))
        provenance_record = {
            "ancestors": ancestors,
            "long_names": long_names,
        }
        provenance_record.update(self.plot_settings[plot_type]["provenance"])
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)
            provenance_logger.log(netcdf_path, provenance_record)

    def _save_data(
        self,
        plot_type: str,
        representative_dataset: dict,
        datasets: dict[str, dict],
        fig: Figure,
    ) -> None:
        """Save plot and netCDF files."""
        # Save single plot file
        plot_path = self._save_plot(plot_type, representative_dataset, fig)
        provenance_record = self._get_provenance_record(
            plot_type,
            representative_dataset,
            list(datasets.values()),
        )
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)

        # Save one netCDF file per dataset
        for label, dataset in datasets.items():
            netcdf_path = self._get_netcdf_path(plot_path, suffix=label)
            io.iris_save(dataset["cube"], netcdf_path)
            provenance_record = self._get_provenance_record(
                plot_type,
                dataset,
                [dataset],
            )
            provenance_record["ancestors"] = dataset["ancestors"]
            with ProvenanceLogger(self.cfg) as provenance_logger:
                provenance_logger.log(netcdf_path, provenance_record)

    def _save_plot(
        self,
        plot_type: str,
        dataset: dict,
        fig: Figure,
    ) -> str:
        """Save plot file."""
        plot_path = self.get_plot_path(plot_type, dataset)
        fig.savefig(plot_path, **self.cfg["savefig_kwargs"])
        logger.info("Wrote %s", plot_path)
        plt.close()
        return plot_path

    def create_1d_benchmarking_plot(
        self,
        plot_type: str,
        datasets: list[dict],
    ) -> None:
        """Create 1D x vs. y benchmarking plot (lines or markers)."""
        fig = plt.figure(**self.cfg["figure_kwargs"])
        axes = fig.add_subplot()

        # Some options are not supported for benchmarking plots
        self.plots[plot_type]["transpose_axes"] = False

        # Plot benchmarking datasets
        benchmark_datasets = self._get_benchmark_datasets(datasets)
        self._plot_1d_data(plot_type, benchmark_datasets, axes)

        # Plot envelope using percentile datasets
        percentile_datasets = self._get_benchmark_percentiles(datasets)
        ylims = axes.get_ylim()
        max_percentile_cube = percentile_datasets[0]["cube"]
        coords = self._check_cube_coords(max_percentile_cube, plot_type)
        coord = max_percentile_cube.coord(coords[0], dim_coords=True)
        if len(percentile_datasets) > 1:
            min_percentile_cube = percentile_datasets[-1]["cube"]
            self._check_cube_coords(min_percentile_cube, plot_type)
        else:
            min_data = np.full(
                max_percentile_cube.shape,
                ylims[0],
                dtype=max_percentile_cube.dtype,
            )
            min_percentile_cube = max_percentile_cube.copy(min_data)
        envelope_kwargs = dict(self.plots[plot_type]["envelope_kwargs"])
        envelope_kwargs["axes"] = axes
        iris.plot.fill_between(
            coord,
            min_percentile_cube,
            max_percentile_cube,
            **envelope_kwargs,
        )

        # Customize plot with user-defined settings
        multi_dataset_facets = self._get_multi_dataset_facets(datasets)
        self._customize_plot(plot_type, axes, multi_dataset_facets)

        # Save data
        self._save_1d_data(plot_type, datasets, fig)

    def create_1d_plot(self, plot_type: str, datasets: list[dict]) -> None:
        """Create 1D x vs. y plot (lines or markers)."""
        fig = plt.figure(**self.cfg["figure_kwargs"])
        axes = fig.add_subplot()
        self._plot_1d_data(plot_type, datasets, axes)
        self._save_1d_data(plot_type, datasets, fig)

    def create_2d_benchmarking_plot(
        self,
        plot_type: str,
        datasets: list[dict],
    ) -> None:
        """Create 2D benchmarking plot."""
        benchmark_datasets = self._get_benchmark_datasets(datasets)
        percentile_datasets = self._get_benchmark_percentiles(datasets)

        # Some options are not supported for benchmarking plots
        self.plots[plot_type]["legend_kwargs"] = False
        self.plots[plot_type]["show_stats"] = False
        self.plots[plot_type]["plot_func"] = "contourf"

        # Create one plot per benchmark dataset
        for dataset in benchmark_datasets:
            fig, axes = self._plot_2d_data_1_panel(plot_type, dataset)

            # Apply hatching (dots) to all points which are not outside the
            # percentile range (the defintion of "outside" depends on the
            # metric)
            hatching_cube = self._get_benchmark_mask(
                dataset,
                percentile_datasets,
            )
            hatching_plot_kwargs = {
                "axes": axes,
                "colors": "none",
                "hatches": ["......"],
                "levels": [0.5, 1.5],
            }
            plot_hatching = self._plot_2d(
                plot_type,
                hatching_cube,
                **hatching_plot_kwargs,
            )
            plot_hatching.set_edgecolor("black")
            plot_hatching.set_linewidth(0.0)

            # Save plot and netCDF files
            save_datasets: dict[str, dict] = {}
            for save_dataset in [dataset, *percentile_datasets]:
                if "_percentile_int" in save_dataset:
                    save_key = f"p{save_dataset['_percentile_int']}"
                else:
                    save_key = ""
                save_datasets[save_key] = save_dataset
            self._save_data(plot_type, dataset, save_datasets, fig)

    def create_2d_plot(self, plot_type: str, datasets: list[dict]) -> None:
        """Create 2D plot."""
        dataset_ref = self._get_reference_dataset(datasets)
        if dataset_ref is not None:
            logger.info(
                "Using reference dataset %s",
                self._get_label(dataset_ref),
            )

        # Some options are not supported for 2D plots
        self.plots[plot_type]["legend_kwargs"] = False

        # Create one plot per (non-reference) dataset
        for dataset in datasets:
            if dataset == dataset_ref:
                continue
            if dataset_ref is None:
                fig, _ = self._plot_2d_data_1_panel(plot_type, dataset)
                save_datasets = {"": dataset}
            else:
                fig = self._plot_2d_data_3_panel(
                    plot_type,
                    dataset,
                    dataset_ref,
                )
                save_datasets = {
                    "_top_left": dataset,
                    "_top_right": dataset_ref,
                    "_bottom": self._get_bias_dataset(dataset, dataset_ref),
                }
            self._save_data(plot_type, dataset, save_datasets, fig)

    def create_benchmarking_boxplot(self) -> None:
        """Create boxplot."""
        plot_type = "benchmarking_boxplot"

        # Respect desired variable order
        given_variables = list(self.grouped_input_data)
        if var_order := self.plots[plot_type]["var_order"]:
            if len(set(var_order)) != len(var_order):
                raise ValueError(
                    f"List of variables given by `var_order` ({var_order}) "
                    f"contains duplicates",
                )
            if set(var_order) != set(given_variables):
                raise ValueError(
                    f"List of variables given by `var_order` ({var_order}) "
                    f"does not agree with given variables ({given_variables})",
                )
        else:
            var_order = given_variables

        # Collect data
        data_idx = 0
        dframe = pd.DataFrame(
            columns=["Variable", "Dataset", "Value", "Benchmark Dataset"],
        )
        all_benchmark_datasets: dict[str, dict] = {}
        for var_key in var_order:
            datasets = self.grouped_input_data[var_key]
            logger.info("Processing variable %s", var_key)

            # Get dataset to be benchmarked
            benchmark_datasets = self._get_benchmark_datasets(datasets)
            if len(benchmark_datasets) > 1:
                raise ValueError(
                    f"Plot {plot_type} only supports a single dataset with "
                    f"'benchmark_dataset: True' for variable '{var_key}', got "
                    f"{len(benchmark_datasets):d}",
                )
            benchmark_dataset = benchmark_datasets[0]
            all_benchmark_datasets[var_key] = benchmark_dataset

            # Collect all datasets in a single data frame
            for dataset in datasets:
                cube = dataset["cube"]
                self._check_cube_coords(cube, plot_type)
                dframe.loc[data_idx] = {
                    "Variable": var_key,
                    "Dataset": self._get_label(dataset),
                    "Value": cube.data.astype(float),
                    "Benchmark Dataset": dataset == benchmark_dataset,
                }
                data_idx = data_idx + 1
        dframe["Value"] = dframe["Value"].astype(float)

        # Create plot
        self._plot_benchmarking_boxplot(dframe, all_benchmark_datasets)

        # Save plot
        self._save_boxplot_data(plot_type, dframe, all_benchmark_datasets)

    def compute(self) -> None:
        """Plot preprocessed data."""
        for plot_type in self.plots:
            plot_settings = self.plot_settings[plot_type]
            plot_function = plot_settings["function"]
            mpl_rc_params = self._get_custom_mpl_rc_params(plot_type)
            logger.info("Plotting %s", plot_type)

            # Handle deprecations
            if "show_x_minor_ticks" in self.plots[plot_type]:
                msg = (
                    f"The option `show_x_minor_ticks` for plot type "
                    f"`{plot_type}` has been deprecated in ESMValTool version "
                    f"2.13.0 and is scheduled for removal in version 2.15.0. "
                    f"Please use the option `x_minor_formatter: null/''` (for "
                    f"`show_x_minor_ticks: false/true`) instead."
                )
                warnings.warn(msg, ESMValToolDeprecationWarning, stacklevel=2)
                show_ticks = self.plots[plot_type].pop("show_x_minor_ticks")
                self.plots[plot_type]["x_minor_formatter"] = (
                    "" if show_ticks else None
                )
            if "show_y_minor_ticks" in self.plots[plot_type]:
                msg = (
                    f"The option `show_y_minor_ticks` for plot type "
                    f"`{plot_type}` has been deprecated in ESMValTool version "
                    f"2.13.0 and is scheduled for removal in version 2.15.0. "
                    f"Please use the option `y_minor_formatter: null/''` (for "
                    f"`show_y_minor_ticks: false/true`) instead."
                )
                warnings.warn(msg, ESMValToolDeprecationWarning, stacklevel=2)
                show_ticks = self.plots[plot_type].pop("show_y_minor_ticks")
                self.plots[plot_type]["y_minor_formatter"] = (
                    "" if show_ticks else None
                )
            if "show_y_minor_ticklabels" in self.plots[plot_type]:
                msg = (
                    f"The option `show_y_minor_ticklabels` for plot type "
                    f"`{plot_type}` has been deprecated in ESMValTool version "
                    f"2.13.0 and is scheduled for removal in version 2.15.0. "
                    f"Please use the option `y_minor_formatter: '%.1f'` "
                    f"instead."
                )
                warnings.warn(msg, ESMValToolDeprecationWarning, stacklevel=2)
                self.plots[plot_type].pop("show_y_minor_ticklabels")
                self.plots[plot_type]["y_minor_formatter"] = "%.1f"
            if "time_on" in self.plots[plot_type]:
                msg = (
                    f"The option `time_on` for plot type `{plot_type}` has "
                    f"been deprecated in ESMValTool version 2.13.0 and is "
                    f"scheduled for removal in version 2.15.0. Please use the "
                    f"option `transpose_axes: true/false` instead."
                )
                warnings.warn(msg, ESMValToolDeprecationWarning, stacklevel=2)
                time_on = self.plots[plot_type].pop("time_on")
                self.plots[plot_type]["transpose_axes"] = time_on == "y-axis"

            # Inspect plot function to determine arguments
            plot_parameters = inspect.signature(plot_function).parameters

            # Plot types where only one plot in total is created
            if not plot_parameters:
                with mpl.rc_context(mpl_rc_params):
                    plot_function()

            # Plot types where multiple plots might be created
            else:
                for var_key, datasets in self.grouped_input_data.items():
                    logger.info("Processing variable %s", var_key)
                    with mpl.rc_context(mpl_rc_params):
                        plot_function(datasets)


def main(cfg: dict) -> None:
    """Run diagnostic."""
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Using DEFAULT_SPHERICAL_EARTH_RADIUS",
            category=UserWarning,
            module="iris",
        )
        MultiDatasets(cfg).compute()


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
