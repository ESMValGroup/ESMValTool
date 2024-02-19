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
    - Time series (plot type ``timeseries``): for each variable separately, all
      datasets are plotted in one single figure. Input data needs to be 1D with
      single dimension `time`.
    - Annual cycle (plot type ``annual_cycle``): for each variable separately,
      all datasets are plotted in one single figure. Input data needs to be 1D
      with single dimension `month_number`.
    - Maps (plot type ``map``): for each variable and dataset, an individual
      map is plotted. If a reference dataset is defined, also include this
      dataset and a bias plot into the figure. Note that if a reference dataset
      is defined, all input datasets need to be given on the same horizontal
      grid (you can use the preprocessor :func:`esmvalcore.preprocessor.regrid`
      for this). Input data needs to be 2D with dimensions `latitude`,
      `longitude`.
    - Zonal mean profiles (plot type ``zonal_mean_profile``):
      for each variable and dataset, an individual profile is plotted. If a
      reference dataset is defined, also include this dataset and a bias plot
      into the figure. Note that if a reference dataset is defined, all input
      datasets need to be given on the same horizontal and vertical grid (you
      can use the preprocessors :func:`esmvalcore.preprocessor.regrid` and
      :func:`esmvalcore.preprocessor.extract_levels` for this). Input data
      needs to be 2D with dimensions `latitude`, `altitude`/`air_pressure`.

      .. warning::

          The plot_type ``profile`` for zonal mean profiles has been deprecated
          in ESMValTool version 2.9.0 and is scheduled for removal in version
          2.11.0. Please use plot type ``zonal_mean_profile`` instead. This is
          an exact replacement.

    - 1D profiles (plot type ``1d_profile``): for each variable separately, all
      datasets are plotted in one single figure. Input data needs to be 1D with
      single dimension `altitude` / `air_pressure`
    - Variable vs. latitude plot (plot type ``variable_vs_lat``):
      for each variable separately, all datasets are plotted in one
      single figure. Input data needs to be 1D with single
      dimension `latitude`.
    - Hovmoeller Z vs. time (plot type ``hovmoeller_z_vs_time``): for each
      variable and dataset, an individual figure is plotted. If a reference
      dataset is defined, also include this dataset and a bias plot into the
      figure. Note that if a reference dataset is defined, all input datasets
      need to be given on the same temporal and vertical grid (you can use
      the preprocessors :func:`esmvalcore.preprocessor.regrid_time` and
      :func:`esmvalcore.preprocessor.extract_levels` for this). Input data
      needs to be 2D with dimensions `time`, `altitude`/`air_pressure`.
    - Hovmoeller time vs. latitude or longitude (plot type
      ``hovmoeller_time_vs_lat_or_lon``): for each variable and dataset, an
      individual figure is plotted. If a reference dataset is defined, also
      include this dataset and a bias plot into the figure. Note that if a
      reference dataset is defined, all input datasets need to be given on the
      same temporal and horizontal grid (you can use the preprocessors
      :func:`esmvalcore.preprocessor.regrid_time` and
      :func:`esmvalcore.preprocessor.regrid` for this). Input data
      needs to be 2D with dimensions `time`, `latitude`/`longitude`.

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
    Optional keyword arguments for :func:`matplotlib.pyplot.figure`. By
    default, uses ``constrained_layout: true``.
group_variables_by: str, optional (default: 'short_name')
    Facet which is used to create variable groups. For each variable group, an
    individual plot is created.
plots: dict, optional
    Plot types plotted by this diagnostic (see list above). Dictionary keys
    must be ``timeseries``, ``annual_cycle``, ``map``, ``zonal_mean_profile``,
    ``1d_profile``, ``variable_vs_lat``, ``hovmoeller_z_vs_time``,
    ``hovmoeller_time_vs_lat_or_lon``. Dictionary values are dictionaries used
    as options for the corresponding plot. The allowed options for the
    different plot types are given below.
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
time_format: str, optional (default: None)
    :func:`~datetime.datetime.strftime` format string that is used to format
    the time axis using :class:`matplotlib.dates.DateFormatter`. If ``None``,
    use the default formatting imposed by the iris plotting function.

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

Configuration options for plot type ``map``
-------------------------------------------
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
    default, uses ``orientation: horizontal, aspect: 30``.
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
    dictionaries used as keyword arguments for the plot function defined by
    ``plot_func``. String arguments can include facets in curly brackets which
    will be derived from the corresponding dataset, e.g., ``{project}``,
    ``{short_name}``, ``{exp}``. Examples: ``default: {levels: 2}, CMIP6:
    {vmin: 200, vmax: 250}``. In addition to the normalization_ options
    supported by the plot function, the option ``norm: centered`` can be
    specified. In this case, the keywords ``vcenter`` and ``halfrange`` should
    be used instead of ``vmin`` or ``vmax`` (see
    :class:`~matplotlib.colors.CenteredNorm`).
plot_kwargs_bias: dict, optional
    Optional keyword arguments for the plot function defined by ``plot_func``
    for plotting biases. These keyword arguments update (and potentially
    overwrite) the ``plot_kwargs`` for the bias plot. This option has no effect
    if no reference dataset is given. See option ``plot_kwargs`` for more
    details. By default, uses ``cmap: bwr`` and ``norm: centered``.
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
    Plot of {long_name}'``, ``xlabel: '{short_name}'``, ``xlim: [0, 5]``.
rasterize: bool, optional (default: True)
    If ``True``, use rasterization_ for map plots to produce smaller files.
    This is only relevant for vector graphics (e.g., ``output_file_type:
    pdf,svg,ps``).
show_stats: bool, optional (default: True)
    Show basic statistics on the plots.
x_pos_stats_avg: float, optional (default: 0.0)
    Text x-position of average (shown on the left) in Axes coordinates. Can be
    adjusted to avoid overlap with the figure. Only relevant if ``show_stats:
    true``.
x_pos_stats_bias: float, optional (default: 0.92)
    Text x-position of bias statistics (shown on the right) in Axes
    coordinates. Can be adjusted to avoid overlap with the figure. Only
    relevant if ``show_stats: true``.

Configuration options for plot type ``zonal_mean_profile``
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
    latitude and altitude/air_pressure.
plot_kwargs: dict, optional
    Optional keyword arguments for the plot function defined by ``plot_func``.
    Dictionary keys are elements identified by ``facet_used_for_labels`` or
    ``default``, e.g., ``CMIP6`` if ``facet_used_for_labels: project`` or
    ``historical`` if ``facet_used_for_labels: exp``. Dictionary values are
    dictionaries used as keyword arguments for the plot function defined by
    ``plot_func``. String arguments can include facets in curly brackets which
    will be derived from the corresponding dataset, e.g., ``{project}``,
    ``{short_name}``, ``{exp}``. Examples: ``default: {levels: 2}, CMIP6:
    {vmin: 200, vmax: 250}``. In addition to the normalization_ options
    supported by the plot function, the option ``norm: centered`` can be
    specified. In this case, the keywords ``vcenter`` and ``halfrange`` should
    be used instead of ``vmin`` or ``vmax`` (see
    :class:`~matplotlib.colors.CenteredNorm`).
plot_kwargs_bias: dict, optional
    Optional keyword arguments for the plot function defined by ``plot_func``
    for plotting biases. These keyword arguments update (and potentially
    overwrite) the ``plot_kwargs`` for the bias plot. This option has no effect
    if no reference dataset is given. See option ``plot_kwargs`` for more
    details. By default, uses ``cmap: bwr`` and ``norm: centered``.
pyplot_kwargs: dict, optional
    Optional calls to functions of :mod:`matplotlib.pyplot`. Dictionary keys
    are functions of :mod:`matplotlib.pyplot`. Dictionary values are used as
    single argument for these functions. String arguments can include facets in
    curly brackets which will be derived from the corresponding dataset, e.g.,
    ``{project}``, ``{short_name}``, ``{exp}``.  Examples: ``title: 'Awesome
    Plot of {long_name}'``, ``xlabel: '{short_name}'``, ``xlim: [0, 5]``.
rasterize: bool, optional (default: True)
    If ``True``, use rasterization_ for profile plots to produce smaller files.
    This is only relevant for vector graphics (e.g., ``output_file_type:
    pdf,svg,ps``).
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

Configuration options for plot type ``variable_vs_lat``
-------------------------------------------------------
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

Configuration options for plot type ``hovmoeller_z_vs_time``
------------------------------------------------------------
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
    latitude and altitude/air_pressure.
plot_kwargs: dict, optional
    Optional keyword arguments for the plot function defined by ``plot_func``.
    Dictionary keys are elements identified by ``facet_used_for_labels`` or
    ``default``, e.g., ``CMIP6`` if ``facet_used_for_labels: project`` or
    ``historical`` if ``facet_used_for_labels: exp``. Dictionary values are
    dictionaries used as keyword arguments for the plot function defined by
    ``plot_func``. String arguments can include facets in curly brackets which
    will be derived from the corresponding dataset, e.g., ``{project}``,
    ``{short_name}``, ``{exp}``. Examples: ``default: {levels: 2}, CMIP6:
    {vmin: 200, vmax: 250}``. In addition to the normalization_ options
    supported by the plot function, the option ``norm: centered`` can be
    specified. In this case, the keywords ``vcenter`` and ``halfrange`` should
    be used instead of ``vmin`` or ``vmax`` (see
    :class:`~matplotlib.colors.CenteredNorm`).
plot_kwargs_bias: dict, optional
    Optional keyword arguments for the plot function defined by ``plot_func``
    for plotting biases. These keyword arguments update (and potentially
    overwrite) the ``plot_kwargs`` for the bias plot. This option has no effect
    if no reference dataset is given. See option ``plot_kwargs`` for more
    details. By default, uses ``cmap: bwr`` and ``norm: centered``.
pyplot_kwargs: dict, optional
    Optional calls to functions of :mod:`matplotlib.pyplot`. Dictionary keys
    are functions of :mod:`matplotlib.pyplot`. Dictionary values are used as
    single argument for these functions. String arguments can include facets in
    curly brackets which will be derived from the corresponding dataset, e.g.,
    ``{project}``, ``{short_name}``, ``{exp}``.  Examples: ``title: 'Awesome
    Plot of {long_name}'``, ``xlabel: '{short_name}'``, ``xlim: [0, 5]``.
rasterize: bool, optional (default: True)
    If ``True``, use rasterization_ for profile plots to produce smaller files.
    This is only relevant for vector graphics (e.g., ``output_file_type:
    pdf,svg,ps``).
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
time_format: str, optional (default: None)
    :func:`~datetime.datetime.strftime` format string that is used to format
    the time axis using :class:`matplotlib.dates.DateFormatter`. If ``None``,
    use the default formatting imposed by the iris plotting function.

Configuration options for plot type ``hovmoeller_time_vs_lat_or_lon``
---------------------------------------------------------------------
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
    {vmin: 200, vmax: 250}``. In addition to the normalization_ options
    supported by the plot function, the option ``norm: centered`` can be
    specified. In this case, the keywords ``vcenter`` and ``halfrange`` should
    be used instead of ``vmin`` or ``vmax`` (see
    :class:`~matplotlib.colors.CenteredNorm`).
plot_kwargs_bias: dict, optional
    Optional keyword arguments for the plot function defined by ``plot_func``
    for plotting biases. These keyword arguments update (and potentially
    overwrite) the ``plot_kwargs`` for the bias plot. This option has no effect
    if no reference dataset is given. See option ``plot_kwargs`` for more
    details. By default, uses ``cmap: bwr`` and ``norm: centered``.
pyplot_kwargs: dict, optional
    Optional calls to functions of :mod:`matplotlib.pyplot`. Dictionary keys
    are functions of :mod:`matplotlib.pyplot`. Dictionary values are used as
    single argument for these functions. String arguments can include facets in
    curly brackets which will be derived from the corresponding dataset, e.g.,
    ``{project}``, ``{short_name}``, ``{exp}``.  Examples: ``title: 'Awesome
    Plot of {long_name}'``, ``xlabel: '{short_name}'``, ``xlim: [0, 5]``.
rasterize: bool, optional (default: True)
    If ``True``, use rasterization_ for profile plots to produce smaller files.
    This is only relevant for vector graphics (e.g., ``output_file_type:
    pdf,svg,ps``).
show_y_minor_ticks: bool, optional (default: True)
    Show minor ticks for time on the Y axis.
show_x_minor_ticks: bool, optional (default: True)
    Show minor ticks for latitude or longitude on the X axis.
time_format: str, optional (default: None)
    :func:`~datetime.datetime.strftime` format string that is used to format
    the time axis using :class:`matplotlib.dates.DateFormatter`. If ``None``,
    use the default formatting imposed by the iris plotting function.

.. hint::

   Extra arguments given to the recipe are ignored, so it is safe to use yaml
   anchors to share the configuration of common arguments with other monitor
   diagnostic script.

.. _rasterization: https://matplotlib.org/stable/gallery/misc/
   rasterization_demo.html
.. _normalization: https://matplotlib.org/stable/users/explain/colors/
   colormapnorms.html

"""
import logging
import warnings
from copy import deepcopy
from pathlib import Path
from pprint import pformat

import cartopy.crs as ccrs
import iris
import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from iris.analysis.cartography import area_weights
from iris.coord_categorisation import add_year
from iris.coords import AuxCoord
from matplotlib.colors import CenteredNorm
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import (
    AutoMinorLocator,
    FormatStrFormatter,
    LogLocator,
    NullFormatter,
)
from sklearn.metrics import r2_score

import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts.monitor.monitor_base import MonitorBase
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    group_metadata,
    io,
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
        self.cfg.setdefault('figure_kwargs', {'constrained_layout': True})
        self.cfg.setdefault('group_variables_by', 'short_name')
        self.cfg.setdefault('savefig_kwargs', {
            'bbox_inches': 'tight',
            'dpi': 300,
            'orientation': 'landscape',
        })
        self.cfg.setdefault('seaborn_settings', {'style': 'ticks'})
        logger.info("Using facet '%s' to group variables",
                    self.cfg['group_variables_by'])
        logger.info("Using facet '%s' to create labels",
                    self.cfg['facet_used_for_labels'])

        # Load input data
        self.input_data = self._load_and_preprocess_data()
        self.grouped_input_data = group_metadata(
            self.input_data,
            self.cfg['group_variables_by'],
            sort=self.cfg['facet_used_for_labels'],
        )

        if 'profile' in self.plots:
            logger.warning("The plot_type ``profile`` for zonal mean profiles"
                           " has been deprecated in ESMValTool version 2.9.0"
                           " and is scheduled for removal in version 2.11.0."
                           " Please use plot type ``zonal_mean_profile``"
                           " instead. This is an exact replacement.")
            if 'zonal_mean_profile' in self.plots:
                raise ValueError(
                    "Both ``profile`` and ``zonal_mean_profile`` is used."
                    " Please use ``zonal_mean_profile`` only.")
            self.plots['zonal_mean_profile'] = self.plots.pop('profile')

        # Check given plot types and set default settings for them
        self.supported_plot_types = [
            'timeseries',
            'annual_cycle',
            'map',
            'zonal_mean_profile',
            '1d_profile',
            'variable_vs_lat',
            'hovmoeller_z_vs_time',
            'hovmoeller_time_vs_lat_or_lon',
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
                self.plots[plot_type].setdefault('annual_mean_kwargs', {})
                self.plots[plot_type].setdefault('gridline_kwargs', {})
                self.plots[plot_type].setdefault('legend_kwargs', {})
                self.plots[plot_type].setdefault('plot_kwargs', {})
                self.plots[plot_type].setdefault('pyplot_kwargs', {})
                self.plots[plot_type].setdefault('time_format', None)

            elif plot_type == 'annual_cycle':
                self.plots[plot_type].setdefault('gridline_kwargs', {})
                self.plots[plot_type].setdefault('legend_kwargs', {})
                self.plots[plot_type].setdefault('plot_kwargs', {})
                self.plots[plot_type].setdefault('pyplot_kwargs', {})

            elif plot_type == 'map':
                self.plots[plot_type].setdefault(
                    'cbar_label', '{short_name} [{units}]')
                self.plots[plot_type].setdefault(
                    'cbar_label_bias', 'Δ{short_name} [{units}]')
                self.plots[plot_type].setdefault(
                    'cbar_kwargs', {'orientation': 'horizontal', 'aspect': 30}
                )
                self.plots[plot_type].setdefault('cbar_kwargs_bias', {})
                self.plots[plot_type].setdefault('common_cbar', False)
                self.plots[plot_type].setdefault('fontsize', 10)
                self.plots[plot_type].setdefault('gridline_kwargs', {})
                self.plots[plot_type].setdefault('plot_func', 'contourf')
                self.plots[plot_type].setdefault('plot_kwargs', {})
                self.plots[plot_type].setdefault('plot_kwargs_bias', {})
                self.plots[plot_type]['plot_kwargs_bias'].setdefault(
                    'cmap', 'bwr'
                )
                self.plots[plot_type]['plot_kwargs_bias'].setdefault(
                    'norm', 'centered'
                )
                if 'projection' not in self.plots[plot_type]:
                    self.plots[plot_type].setdefault('projection', 'Robinson')
                    self.plots[plot_type].setdefault(
                        'projection_kwargs', {'central_longitude': 10}
                    )
                else:
                    self.plots[plot_type].setdefault('projection_kwargs', {})
                self.plots[plot_type].setdefault('pyplot_kwargs', {})
                self.plots[plot_type].setdefault('rasterize', True)
                self.plots[plot_type].setdefault('show_stats', True)
                self.plots[plot_type].setdefault('x_pos_stats_avg', 0.0)
                self.plots[plot_type].setdefault('x_pos_stats_bias', 0.92)

            elif plot_type == 'zonal_mean_profile':
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
                self.plots[plot_type]['plot_kwargs_bias'].setdefault(
                    'norm', 'centered'
                )
                self.plots[plot_type].setdefault('pyplot_kwargs', {})
                self.plots[plot_type].setdefault('rasterize', True)
                self.plots[plot_type].setdefault('show_stats', True)
                self.plots[plot_type].setdefault(
                    'show_y_minor_ticklabels', False
                )
                self.plots[plot_type].setdefault('x_pos_stats_avg', 0.01)
                self.plots[plot_type].setdefault('x_pos_stats_bias', 0.7)

            elif plot_type == '1d_profile':
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
            elif plot_type == 'variable_vs_lat':
                self.plots[plot_type].setdefault('gridline_kwargs', {})
                self.plots[plot_type].setdefault('legend_kwargs', {})
                self.plots[plot_type].setdefault('plot_kwargs', {})
                self.plots[plot_type].setdefault('pyplot_kwargs', {})

            elif plot_type == 'hovmoeller_z_vs_time':
                self.plots[plot_type].setdefault('cbar_label',
                                                 '{short_name} [{units}]')
                self.plots[plot_type].setdefault('cbar_label_bias',
                                                 'Δ{short_name} [{units}]')
                self.plots[plot_type].setdefault('cbar_kwargs',
                                                 {'orientation': 'vertical'})
                self.plots[plot_type].setdefault('cbar_kwargs_bias', {})
                self.plots[plot_type].setdefault('common_cbar', False)
                self.plots[plot_type].setdefault('fontsize', 10)
                self.plots[plot_type].setdefault('log_y', True)
                self.plots[plot_type].setdefault('plot_func', 'contourf')
                self.plots[plot_type].setdefault('plot_kwargs', {})
                self.plots[plot_type].setdefault('plot_kwargs_bias', {})
                self.plots[plot_type]['plot_kwargs_bias'].setdefault(
                    'cmap', 'bwr')
                self.plots[plot_type]['plot_kwargs_bias'].setdefault(
                    'norm', 'centered'
                )
                self.plots[plot_type].setdefault('pyplot_kwargs', {})
                self.plots[plot_type].setdefault('rasterize', True)
                self.plots[plot_type].setdefault('show_stats', True)
                self.plots[plot_type].setdefault('show_y_minor_ticklabels',
                                                 False)
                self.plots[plot_type].setdefault('time_format', None)
                self.plots[plot_type].setdefault('x_pos_stats_avg', 0.01)
                self.plots[plot_type].setdefault('x_pos_stats_bias', 0.7)

            elif plot_type == 'hovmoeller_time_vs_lat_or_lon':
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
                self.plots[plot_type].setdefault('plot_func', 'contourf')
                self.plots[plot_type].setdefault('plot_kwargs', {})
                self.plots[plot_type].setdefault('plot_kwargs_bias', {})
                self.plots[plot_type]['plot_kwargs_bias'].setdefault(
                    'cmap', 'bwr'
                )
                self.plots[plot_type]['plot_kwargs_bias'].setdefault(
                    'norm', 'centered'
                )
                self.plots[plot_type].setdefault('pyplot_kwargs', {})
                self.plots[plot_type].setdefault('rasterize', True)
                self.plots[plot_type].setdefault(
                    'show_y_minor_ticks', True
                )
                self.plots[plot_type].setdefault(
                    'show_x_minor_ticks', True
                )
                self.plots[plot_type].setdefault('time_format', None)

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

    def _add_stats(self, plot_type, axes, dim_coords, dataset,
                   ref_dataset=None):
        """Add text to plot that describes basic statistics."""
        if not self.plots[plot_type]['show_stats']:
            return

        # Extract cube(s)
        cube = dataset['cube']
        if ref_dataset is None:
            ref_cube = None
            label = self._get_label(dataset)
        else:
            ref_cube = ref_dataset['cube']
            label = (f'{self._get_label(dataset)} vs. '
                     f'{self._get_label(ref_dataset)}')

        # Different options for the different plots types
        fontsize = 6.0
        y_pos = 0.95
        if all([
                'x_pos_stats_avg' in self.plots[plot_type],
                'x_pos_stats_bias' in self.plots[plot_type],
        ]):
            x_pos_bias = self.plots[plot_type]['x_pos_stats_bias']
            x_pos = self.plots[plot_type]['x_pos_stats_avg']
        else:
            raise NotImplementedError(f"plot_type '{plot_type}' not supported")

        # For zonal_mean_profile plots add scalar longitude coordinate
        # (necessary for calculation of area weights). The exact values for the
        # points/bounds of this coordinate do not matter since they don't
        # change the weights.
        if not cube.coords('longitude'):
            lon_coord = AuxCoord(
                180.0,
                bounds=[0.0, 360.0],
                var_name='lon',
                standard_name='longitude',
                long_name='longitude',
                units='degrees_east',
            )
            cube.add_aux_coord(lon_coord, ())

        # Mean
        weights = area_weights(cube)
        if ref_cube is None:
            mean = cube.collapsed(dim_coords, iris.analysis.MEAN,
                                  weights=weights)
            logger.info(
                "Area-weighted mean of %s for %s = %f%s",
                dataset['short_name'],
                label,
                mean.data,
                dataset['units'],
            )
        else:
            mean = (cube - ref_cube).collapsed(dim_coords, iris.analysis.MEAN,
                                               weights=weights)
            logger.info(
                "Area-weighted bias of %s for %s = %f%s",
                dataset['short_name'],
                label,
                mean.data,
                dataset['units'],
            )
        if np.abs(mean.data) >= 0.1:
            mean_val = f"{mean.data:.2f} {cube.units}"
        else:
            mean_val = f"{mean.data:.2e} {cube.units}"
        axes.text(
            x_pos, y_pos, mean_val, fontsize=fontsize, transform=axes.transAxes
        )
        if ref_cube is None:
            return

        # Weighted RMSE
        rmse = (cube - ref_cube).collapsed(dim_coords, iris.analysis.RMS,
                                           weights=weights)
        if np.abs(rmse.data) >= 0.1:
            rmse_val = f"{rmse.data:.2f} {cube.units}"
        else:
            rmse_val = f"{rmse.data:.2e} {cube.units}"
        axes.text(x_pos_bias, y_pos, f"RMSE={rmse_val}",
                  fontsize=fontsize, transform=axes.transAxes)
        logger.info(
            "Area-weighted RMSE of %s for %s = %f%s",
            dataset['short_name'],
            label,
            rmse.data,
            dataset['units'],
        )

        # Weighted R2
        mask = np.ma.getmaskarray(cube.data).ravel()
        mask |= np.ma.getmaskarray(ref_cube.data).ravel()
        cube_data = cube.data.ravel()[~mask]
        ref_cube_data = ref_cube.data.ravel()[~mask]
        weights = weights.ravel()[~mask]
        r2_val = r2_score(cube_data, ref_cube_data, sample_weight=weights)
        axes.text(x_pos_bias, y_pos - 0.1, rf"R$^2$={r2_val:.2f}",
                  fontsize=fontsize, transform=axes.transAxes)
        logger.info(
            "Area-weighted R2 of %s for %s = %f",
            dataset['short_name'],
            label,
            r2_val,
        )

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

    def _get_map_projection(self):
        """Get projection used for map plots."""
        plot_type = 'map'
        projection = self.plots[plot_type]['projection']
        projection_kwargs = self.plots[plot_type]['projection_kwargs']

        # Check if desired projection is valid
        if not hasattr(ccrs, projection):
            raise AttributeError(
                f"Got invalid projection '{projection}' for plotting "
                f"{plot_type}, expected class of cartopy.crs")

        return getattr(ccrs, projection)(**projection_kwargs)

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
        if plot_type in ('timeseries', 'annual_cycle', '1d_profile',
                         'variable_vs_lat'):
            plot_kwargs.setdefault('label', label)

        if plot_kwargs.get('norm') == 'centered':
            norm = CenteredNorm(
                vcenter=plot_kwargs.pop('vcenter', 0.0),
                halfrange=plot_kwargs.pop('halfrange', None),
            )
            plot_kwargs['norm'] = norm

        return deepcopy(plot_kwargs)

    def _load_and_preprocess_data(self):
        """Load and preprocess data."""
        input_data = list(self.cfg['input_data'].values())

        for dataset in input_data:
            filename = dataset['filename']
            logger.info("Loading %s", filename)
            cube = iris.load_cube(filename)

            # Fix time coordinate if present
            if cube.coords('time', dim_coords=True):
                ih.unify_time_coord(cube)

            # Fix Z-coordinate if present
            if cube.coords('air_pressure', dim_coords=True):
                z_coord = cube.coord('air_pressure', dim_coords=True)
                z_coord.attributes['positive'] = 'down'
                z_coord.convert_units('hPa')
            elif cube.coords('altitude', dim_coords=True):
                z_coord = cube.coord('altitude')
                z_coord.attributes['positive'] = 'up'

            dataset['cube'] = cube

        return input_data

    def _plot_map_with_ref(self, plot_func, dataset, ref_dataset):
        """Plot map plot for single dataset with a reference dataset."""
        plot_type = 'map'
        logger.info("Plotting map with reference dataset '%s' for '%s'",
                    self._get_label(ref_dataset), self._get_label(dataset))

        # Make sure that the data has the correct dimensions
        cube = dataset['cube']
        ref_cube = ref_dataset['cube']
        dim_coords_dat = self._check_cube_dimensions(cube, plot_type)
        dim_coords_ref = self._check_cube_dimensions(ref_cube, plot_type)

        # Create single figure with multiple axes
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg['figure_kwargs'])
            gridspec = GridSpec(5, 4, figure=fig,
                                height_ratios=[1.0, 1.0, 0.4, 1.0, 1.0])

            # Options used for all subplots
            projection = self._get_map_projection()
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            gridline_kwargs = self._get_gridline_kwargs(plot_type)
            fontsize = self.plots[plot_type]['fontsize']

            # Plot dataset (top left)
            axes_data = fig.add_subplot(gridspec[0:2, 0:2],
                                        projection=projection)
            plot_kwargs['axes'] = axes_data
            plot_data = plot_func(cube, **plot_kwargs)
            axes_data.coastlines()
            if gridline_kwargs is not False:
                axes_data.gridlines(**gridline_kwargs)
            axes_data.set_title(self._get_label(dataset), pad=3.0)
            self._add_stats(plot_type, axes_data, dim_coords_dat, dataset)
            self._process_pyplot_kwargs(plot_type, dataset)

            # Plot reference dataset (top right)
            # Note: make sure to use the same vmin and vmax than the top left
            # plot if a common cpltolorbar is desired
            axes_ref = fig.add_subplot(gridspec[0:2, 2:4],
                                       projection=projection)
            plot_kwargs['axes'] = axes_ref
            if self.plots[plot_type]['common_cbar']:
                plot_kwargs.setdefault('vmin', plot_data.get_clim()[0])
                plot_kwargs.setdefault('vmax', plot_data.get_clim()[1])
            plot_ref = plot_func(ref_cube, **plot_kwargs)
            axes_ref.coastlines()
            if gridline_kwargs is not False:
                axes_ref.gridlines(**gridline_kwargs)
            axes_ref.set_title(self._get_label(ref_dataset), pad=3.0)
            self._add_stats(plot_type, axes_ref, dim_coords_ref, ref_dataset)
            self._process_pyplot_kwargs(plot_type, ref_dataset)

            # Add colorbar(s)
            self._add_colorbar(plot_type, plot_data, plot_ref, axes_data,
                               axes_ref, dataset, ref_dataset)

            # Plot bias (bottom center)
            bias_cube = cube - ref_cube
            axes_bias = fig.add_subplot(gridspec[3:5, 1:3],
                                        projection=projection)
            plot_kwargs_bias = self._get_plot_kwargs(plot_type, dataset,
                                                     bias=True)
            plot_kwargs_bias['axes'] = axes_bias
            plot_bias = plot_func(bias_cube, **plot_kwargs_bias)
            axes_bias.coastlines()
            if gridline_kwargs is not False:
                axes_bias.gridlines(**gridline_kwargs)
            axes_bias.set_title(
                f"{self._get_label(dataset)} - {self._get_label(ref_dataset)}",
                pad=3.0,
            )
            cbar_kwargs_bias = self._get_cbar_kwargs(plot_type, bias=True)
            cbar_bias = fig.colorbar(plot_bias, ax=axes_bias,
                                     **cbar_kwargs_bias)
            cbar_bias.set_label(
                self._get_cbar_label(plot_type, dataset, bias=True),
                fontsize=fontsize,
            )
            cbar_bias.ax.tick_params(labelsize=fontsize)
            self._add_stats(plot_type, axes_bias, dim_coords_dat, dataset,
                            ref_dataset)

            # Customize plot
            fig.suptitle(dataset['long_name'])
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

    def _plot_map_without_ref(self, plot_func, dataset):
        """Plot map plot for single dataset without a reference dataset."""
        plot_type = 'map'
        logger.info("Plotting map without reference dataset for '%s'",
                    self._get_label(dataset))

        # Make sure that the data has the correct dimensions
        cube = dataset['cube']
        dim_coords_dat = self._check_cube_dimensions(cube, plot_type)

        # Create plot with desired settings
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg['figure_kwargs'])
            axes = fig.add_subplot(projection=self._get_map_projection())
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs['axes'] = axes
            plot_map = plot_func(cube, **plot_kwargs)
            axes.coastlines()
            gridline_kwargs = self._get_gridline_kwargs(plot_type)
            if gridline_kwargs is not False:
                axes.gridlines(**gridline_kwargs)

            # Print statistics if desired
            self._add_stats(plot_type, axes, dim_coords_dat, dataset)

            # Setup colorbar
            fontsize = self.plots[plot_type]['fontsize']
            colorbar = fig.colorbar(plot_map, ax=axes,
                                    **self._get_cbar_kwargs(plot_type))
            colorbar.set_label(self._get_cbar_label(plot_type, dataset),
                               fontsize=fontsize)
            colorbar.ax.tick_params(labelsize=fontsize)

            # Customize plot
            axes.set_title(self._get_label(dataset))
            fig.suptitle(dataset['long_name'])
            self._process_pyplot_kwargs(plot_type, dataset)

            # Rasterization
            if self.plots[plot_type]['rasterize']:
                self._set_rasterized([axes])

        # File paths
        plot_path = self.get_plot_path(plot_type, dataset)
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)

        return (plot_path, {netcdf_path: cube})

    def _plot_zonal_mean_profile_with_ref(self, plot_func, dataset,
                                          ref_dataset):
        """Plot zonal mean profile for single dataset with reference."""
        plot_type = 'zonal_mean_profile'
        logger.info("Plotting zonal mean profile with reference dataset"
                    " '%s' for '%s'",
                    self._get_label(ref_dataset), self._get_label(dataset))

        # Make sure that the data has the correct dimensions
        cube = dataset['cube']
        ref_cube = ref_dataset['cube']
        dim_coords_dat = self._check_cube_dimensions(cube, plot_type)
        dim_coords_ref = self._check_cube_dimensions(ref_cube, plot_type)

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
            self._add_stats(plot_type, axes_data, dim_coords_dat, dataset)
            self._process_pyplot_kwargs(plot_type, dataset)

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
            self._add_stats(plot_type, axes_ref, dim_coords_ref, ref_dataset)
            self._process_pyplot_kwargs(plot_type, ref_dataset)

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
            self._add_stats(plot_type, axes_bias, dim_coords_dat, dataset,
                            ref_dataset)

            # Customize plot
            fig.suptitle(dataset['long_name'])
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

    def _plot_zonal_mean_profile_without_ref(self, plot_func, dataset):
        """Plot zonal mean profile for single dataset without reference."""
        plot_type = 'zonal_mean_profile'
        logger.info("Plotting zonal mean profile without reference dataset"
                    " for '%s'",
                    self._get_label(dataset))

        # Make sure that the data has the correct dimensions
        cube = dataset['cube']
        dim_coords_dat = self._check_cube_dimensions(cube, plot_type)

        # Create plot with desired settings
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg['figure_kwargs'])
            axes = fig.add_subplot()
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs['axes'] = axes
            plot_zonal_mean_profile = plot_func(cube, **plot_kwargs)

            # Print statistics if desired
            self._add_stats(plot_type, axes, dim_coords_dat, dataset)

            # Setup colorbar
            fontsize = self.plots[plot_type]['fontsize']
            colorbar = fig.colorbar(plot_zonal_mean_profile, ax=axes,
                                    **self._get_cbar_kwargs(plot_type))
            colorbar.set_label(self._get_cbar_label(plot_type, dataset),
                               fontsize=fontsize)
            colorbar.ax.tick_params(labelsize=fontsize)

            # Customize plot
            axes.set_title(self._get_label(dataset))
            fig.suptitle(dataset['long_name'])
            axes.set_xlabel('latitude [°N]')
            z_coord = cube.coord(axis='Z')
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

    def _plot_hovmoeller_z_vs_time_without_ref(self, plot_func, dataset):
        """Plot Hovmoeller Z vs. time for single dataset without reference."""
        plot_type = 'hovmoeller_z_vs_time'
        logger.info(
            "Plotting Hovmoeller Z vs. time without reference dataset"
            " for '%s'", self._get_label(dataset))

        # Make sure that the data has the correct dimensions
        cube = dataset['cube']
        dim_coords_dat = self._check_cube_dimensions(cube, plot_type)

        # Create plot with desired settings
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg['figure_kwargs'])
            axes = fig.add_subplot()
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs['axes'] = axes
            plot_hovmoeller = plot_func(cube, **plot_kwargs)

            # Print statistics if desired
            self._add_stats(plot_type, axes, dim_coords_dat, dataset)

            # Setup colorbar
            fontsize = self.plots[plot_type]['fontsize']
            colorbar = fig.colorbar(plot_hovmoeller,
                                    ax=axes,
                                    **self._get_cbar_kwargs(plot_type))
            colorbar.set_label(self._get_cbar_label(plot_type, dataset),
                               fontsize=fontsize)
            colorbar.ax.tick_params(labelsize=fontsize)

            # Customize plot
            axes.set_title(self._get_label(dataset))
            fig.suptitle(dataset['long_name'])
            z_coord = cube.coord(axis='Z')
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
            if self.plots[plot_type]['time_format'] is not None:
                axes.get_xaxis().set_major_formatter(
                    mdates.DateFormatter(self.plots[plot_type]['time_format']))
            axes.set_xlabel('time')
            self._process_pyplot_kwargs(plot_type, dataset)

            # Rasterization
            if self.plots[plot_type]['rasterize']:
                self._set_rasterized([axes])

        # File paths
        plot_path = self.get_plot_path(plot_type, dataset)
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)

        return (plot_path, {netcdf_path: cube})

    def _plot_hovmoeller_z_vs_time_with_ref(self, plot_func, dataset,
                                            ref_dataset):
        """Plot Hovmoeller Z vs. time for single dataset with reference."""
        plot_type = 'hovmoeller_z_vs_time'
        logger.info(
            "Plotting Hovmoeller z vs. time with reference dataset"
            " '%s' for '%s'", self._get_label(ref_dataset),
            self._get_label(dataset))

        # Make sure that the data has the correct dimensions
        cube = dataset['cube']
        ref_cube = ref_dataset['cube']
        dim_coords_dat = self._check_cube_dimensions(cube, plot_type)
        dim_coords_ref = self._check_cube_dimensions(ref_cube, plot_type)

        # Create single figure with multiple axes
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg['figure_kwargs'])
            gridspec = GridSpec(5,
                                4,
                                figure=fig,
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
            if self.plots[plot_type]['time_format'] is not None:
                axes_data.get_xaxis().set_major_formatter(
                    mdates.DateFormatter(self.plots[plot_type]['time_format']))
            self._add_stats(plot_type, axes_data, dim_coords_dat, dataset)
            self._process_pyplot_kwargs(plot_type, dataset)

            # Plot reference dataset (top right)
            # Note: make sure to use the same vmin and vmax than the top left
            # plot if a common colorbar is desired
            axes_ref = fig.add_subplot(gridspec[0:2, 2:4],
                                       sharex=axes_data,
                                       sharey=axes_data)
            plot_kwargs['axes'] = axes_ref
            if self.plots[plot_type]['common_cbar']:
                plot_kwargs.setdefault('vmin', plot_data.get_clim()[0])
                plot_kwargs.setdefault('vmax', plot_data.get_clim()[1])
            plot_ref = plot_func(ref_cube, **plot_kwargs)
            axes_ref.set_title(self._get_label(ref_dataset), pad=3.0)
            plt.setp(axes_ref.get_yticklabels(), visible=False)
            self._add_stats(plot_type, axes_ref, dim_coords_ref, ref_dataset)
            self._process_pyplot_kwargs(plot_type, ref_dataset)

            # Add colorbar(s)
            self._add_colorbar(plot_type, plot_data, plot_ref, axes_data,
                               axes_ref, dataset, ref_dataset)

            # Plot bias (bottom center)
            bias_cube = cube - ref_cube
            axes_bias = fig.add_subplot(gridspec[3:5, 1:3],
                                        sharex=axes_data,
                                        sharey=axes_data)
            plot_kwargs_bias = self._get_plot_kwargs(plot_type,
                                                     dataset,
                                                     bias=True)
            plot_kwargs_bias['axes'] = axes_bias
            plot_bias = plot_func(bias_cube, **plot_kwargs_bias)
            axes_bias.set_title(
                f"{self._get_label(dataset)} - {self._get_label(ref_dataset)}",
                pad=3.0,
            )
            axes_bias.set_xlabel('time')
            axes_bias.set_ylabel(f'{z_coord.long_name} [{z_coord.units}]')
            cbar_kwargs_bias = self._get_cbar_kwargs(plot_type, bias=True)
            cbar_bias = fig.colorbar(plot_bias,
                                     ax=axes_bias,
                                     **cbar_kwargs_bias)
            cbar_bias.set_label(
                self._get_cbar_label(plot_type, dataset, bias=True),
                fontsize=fontsize,
            )
            cbar_bias.ax.tick_params(labelsize=fontsize)
            self._add_stats(plot_type, axes_bias, dim_coords_dat, dataset,
                            ref_dataset)

            # Customize plot
            fig.suptitle(dataset['long_name'])
            self._process_pyplot_kwargs(plot_type, dataset)

            # Rasterization
            if self.plots[plot_type]['rasterize']:
                self._set_rasterized([axes_data, axes_ref, axes_bias])

        # File paths
        plot_path = self.get_plot_path(plot_type, dataset)
        netcdf_path = (get_diagnostic_filename(
            Path(plot_path).stem + "_{pos}", self.cfg))
        netcdf_paths = {
            netcdf_path.format(pos='top_left'): cube,
            netcdf_path.format(pos='top_right'): ref_cube,
            netcdf_path.format(pos='bottom'): bias_cube,
        }

        return (plot_path, netcdf_paths)

    def _plot_hovmoeller_time_vs_lat_or_lon_with_ref(self, plot_func, dataset,
                                                     ref_dataset):
        """Plot the hovmoeller profile for single dataset with reference."""
        plot_type = 'hovmoeller_time_vs_lat_or_lon'
        logger.info("Plotting Hovmoeller plots with reference dataset"
                    " '%s' for '%s'",
                    self._get_label(ref_dataset), self._get_label(dataset))

        # Make sure that the data has the correct dimensions
        cube = dataset['cube']
        ref_cube = ref_dataset['cube']
        dim_coords_dat = self._check_cube_dimensions(cube, plot_type)
        self._check_cube_dimensions(ref_cube, plot_type)

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
            coord_names = [coord[0].name() for coord in cube.dim_coords]
            if coord_names[0] == "time":
                coord_names.reverse()
            plot_kwargs['coords'] = coord_names
            plot_data = plot_func(cube, **plot_kwargs)
            axes_data.set_title(self._get_label(dataset), pad=3.0)
            axes_data.set_ylabel('time')
            if self.plots[plot_type]['time_format'] is not None:
                axes_data.get_yaxis().set_major_formatter(mdates.DateFormatter(
                    self.plots[plot_type]['time_format']))
            if self.plots[plot_type]['show_y_minor_ticks']:
                axes_data.get_yaxis().set_minor_locator(AutoMinorLocator())
            if self.plots[plot_type]['show_x_minor_ticks']:
                axes_data.get_xaxis().set_minor_locator(AutoMinorLocator())
            self._process_pyplot_kwargs(plot_type, dataset)

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
            self._process_pyplot_kwargs(plot_type, ref_dataset)

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
            plot_kwargs_bias['coords'] = coord_names
            plot_bias = plot_func(bias_cube, **plot_kwargs_bias)
            axes_bias.set_title(
                f"{self._get_label(dataset)} - {self._get_label(ref_dataset)}",
                pad=3.0,
            )
            axes_bias.set_ylabel('time')
            if 'latitude' in dim_coords_dat:
                axes_bias.set_xlabel('latitude [°N]')
            elif 'longitude' in dim_coords_dat:
                axes_bias.set_xlabel('longitude [°E]')
            cbar_kwargs_bias = self._get_cbar_kwargs(plot_type, bias=True)
            cbar_bias = fig.colorbar(plot_bias, ax=axes_bias,
                                     **cbar_kwargs_bias)
            cbar_bias.set_label(
                self._get_cbar_label(plot_type, dataset, bias=True),
                fontsize=fontsize,
            )
            cbar_bias.ax.tick_params(labelsize=fontsize)

            # Customize plot
            fig.suptitle(dataset['long_name'])
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

    def _plot_hovmoeller_time_vs_lat_or_lon_without_ref(self, plot_func,
                                                        dataset):
        """Plot time vs zonal or meridional Hovmoeller without reference."""
        plot_type = 'hovmoeller_time_vs_lat_or_lon'
        logger.info("Plotting Hovmoeller plots without reference dataset"
                    " for '%s'", self._get_label(dataset))

        # Make sure that the data has the correct dimensions
        cube = dataset['cube']
        dim_coords_dat = self._check_cube_dimensions(cube, plot_type)

        # Create plot with desired settings
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg['figure_kwargs'])
            axes = fig.add_subplot()
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs['axes'] = axes

            # Make sure time is on y-axis
            plot_kwargs['coords'] = list(reversed(dim_coords_dat))
            plot_hovmoeller = plot_func(cube, **plot_kwargs)

            # Setup colorbar
            fontsize = self.plots[plot_type]['fontsize']
            colorbar = fig.colorbar(plot_hovmoeller, ax=axes,
                                    **self._get_cbar_kwargs(plot_type))
            colorbar.set_label(self._get_cbar_label(plot_type, dataset),
                               fontsize=fontsize)
            colorbar.ax.tick_params(labelsize=fontsize)

            # Customize plot
            axes.set_title(self._get_label(dataset))
            fig.suptitle(dataset['long_name'])
            if 'latitude' in dim_coords_dat:
                axes.set_xlabel('latitude [°N]')
            elif 'longitude' in dim_coords_dat:
                axes.set_xlabel('longitude [°E]')
            axes.set_ylabel('time')
            if self.plots[plot_type]['time_format'] is not None:
                axes.get_yaxis().set_major_formatter(mdates.DateFormatter(
                    self.plots[plot_type]['time_format'])
                )
            if self.plots[plot_type]['show_y_minor_ticks']:
                axes.get_yaxis().set_minor_locator(AutoMinorLocator())
            if self.plots[plot_type]['show_x_minor_ticks']:
                axes.get_xaxis().set_minor_locator(AutoMinorLocator())
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

    @staticmethod
    def _check_cube_dimensions(cube, plot_type):
        """Check that cube has correct dimensional variables."""
        expected_dimensions_dict = {
            'annual_cycle': (['month_number'],),
            'map': (['latitude', 'longitude'],),
            'zonal_mean_profile': (['latitude', 'air_pressure'],
                                   ['latitude', 'altitude']),
            'timeseries': (['time'],),
            '1d_profile': (['air_pressure'],
                           ['altitude']),
            'variable_vs_lat': (['latitude'],),
            'hovmoeller_z_vs_time': (['time', 'air_pressure'],
                                     ['time', 'altitude']),
            'hovmoeller_time_vs_lat_or_lon': (['time', 'latitude'],
                                              ['time', 'longitude']),
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
                multi_dataset_facets[key] = datasets[0].get(key)
            else:
                multi_dataset_facets[key] = f'ambiguous_{key}'
        return multi_dataset_facets

    def _get_reference_dataset(self, datasets):
        """Extract reference dataset."""
        variable = datasets[0][self.cfg['group_variables_by']]
        ref_datasets = [d for d in datasets if
                        d.get('reference_for_monitor_diags', False)]
        if len(ref_datasets) > 1:
            raise ValueError(
                f"Expected at most 1 reference dataset (with "
                f"'reference_for_monitor_diags: true' for variable "
                f"'{variable}', got {len(ref_datasets):d}")
        if ref_datasets:
            return ref_datasets[0]
        return None

    def create_timeseries_plot(self, datasets):
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
        cubes = {}
        for dataset in datasets:
            ancestors.append(dataset['filename'])
            cube = dataset['cube']
            cubes[self._get_label(dataset)] = cube
            self._check_cube_dimensions(cube, plot_type)

            # Plot original time series
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs['axes'] = axes
            iris.plot.plot(cube, **plot_kwargs)

            # Plot annual means if desired
            annual_mean_kwargs = self.plots[plot_type]['annual_mean_kwargs']
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
        axes.set_title(multi_dataset_facets['long_name'])
        axes.set_xlabel('time')
        # apply time formatting
        if self.plots[plot_type]['time_format'] is not None:
            axes.get_xaxis().set_major_formatter(
                mdates.DateFormatter(self.plots[plot_type]['time_format']))
        axes.set_ylabel(
            f"{multi_dataset_facets[self.cfg['group_variables_by']]} "
            f"[{multi_dataset_facets['units']}]"
        )
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
            n: datasets[0][n] for n in ('short_name', 'long_name', 'units')
        }
        io.save_1d_data(cubes, netcdf_path, 'time', var_attrs)

        # Provenance tracking
        caption = (f"Time series of {multi_dataset_facets['long_name']} for "
                   f"various datasets.")
        provenance_record = {
            'ancestors': ancestors,
            'authors': ['schlund_manuel'],
            'caption': caption,
            'plot_types': ['line'],
            'long_names': [var_attrs['long_name']],
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)
            provenance_logger.log(netcdf_path, provenance_record)

    def create_annual_cycle_plot(self, datasets):
        """Create annual cycle plot."""
        plot_type = 'annual_cycle'
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        logger.info("Plotting %s", plot_type)
        fig = plt.figure(**self.cfg['figure_kwargs'])
        axes = fig.add_subplot()

        # Plot all datasets in one single figure
        ancestors = []
        cubes = {}
        for dataset in datasets:
            ancestors.append(dataset['filename'])
            cube = dataset['cube']
            cubes[self._get_label(dataset)] = cube
            self._check_cube_dimensions(cube, plot_type)

            # Plot annual cycle
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs['axes'] = axes
            iris.plot.plot(cube, **plot_kwargs)

        # Default plot appearance
        multi_dataset_facets = self._get_multi_dataset_facets(datasets)
        axes.set_title(multi_dataset_facets['long_name'])
        axes.set_xlabel('Month')
        axes.set_ylabel(
            f"{multi_dataset_facets[self.cfg['group_variables_by']]} "
            f"[{multi_dataset_facets['units']}]"
        )
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
            n: datasets[0][n] for n in ('short_name', 'long_name', 'units')
        }
        io.save_1d_data(cubes, netcdf_path, 'month_number', var_attrs)

        # Provenance tracking
        caption = (f"Annual cycle of {multi_dataset_facets['long_name']} for "
                   f"various datasets.")
        provenance_record = {
            'ancestors': ancestors,
            'authors': ['schlund_manuel'],
            'caption': caption,
            'plot_types': ['seas'],
            'long_names': [var_attrs['long_name']],
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)
            provenance_logger.log(netcdf_path, provenance_record)

    def create_map_plot(self, datasets):
        """Create map plot."""
        plot_type = 'map'
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        # Get reference dataset if possible
        ref_dataset = self._get_reference_dataset(datasets)
        if ref_dataset is None:
            logger.info("Plotting %s without reference dataset", plot_type)
        else:
            logger.info("Plotting %s with reference dataset '%s'", plot_type,
                        self._get_label(ref_dataset))

        # Get plot function
        plot_func = self._get_plot_func(plot_type)

        # Create a single plot for each dataset (incl. reference dataset if
        # given)
        for dataset in datasets:
            if dataset == ref_dataset:
                continue
            ancestors = [dataset['filename']]
            if ref_dataset is None:
                (plot_path, netcdf_paths) = (
                    self._plot_map_without_ref(plot_func, dataset)
                )
                caption = (
                    f"Map plot of {dataset['long_name']} of dataset "
                    f"{dataset['dataset']} (project {dataset['project']})."
                )
            else:
                (plot_path, netcdf_paths) = (
                    self._plot_map_with_ref(plot_func, dataset, ref_dataset)
                )
                caption = (
                    f"Map plot of {dataset['long_name']} of dataset "
                    f"{dataset['dataset']} (project {dataset['project']}) "
                    f"including bias relative to {ref_dataset['dataset']} "
                    f"(project {ref_dataset['project']})."
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
                'authors': ['schlund_manuel'],
                'caption': caption,
                'plot_types': ['map'],
                'long_names': [dataset['long_name']],
            }
            with ProvenanceLogger(self.cfg) as provenance_logger:
                provenance_logger.log(plot_path, provenance_record)
                for netcdf_path in netcdf_paths:
                    provenance_logger.log(netcdf_path, provenance_record)

    def create_zonal_mean_profile_plot(self, datasets):
        """Create zonal mean profile plot."""
        plot_type = 'zonal_mean_profile'
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        # Get reference dataset if possible
        ref_dataset = self._get_reference_dataset(datasets)
        if ref_dataset is None:
            logger.info("Plotting %s without reference dataset", plot_type)
        else:
            logger.info("Plotting %s with reference dataset '%s'", plot_type,
                        self._get_label(ref_dataset))

        # Get plot function
        plot_func = self._get_plot_func(plot_type)

        # Create a single plot for each dataset (incl. reference dataset if
        # given)
        for dataset in datasets:
            if dataset == ref_dataset:
                continue
            ancestors = [dataset['filename']]
            if ref_dataset is None:
                (plot_path, netcdf_paths) = (
                    self._plot_zonal_mean_profile_without_ref(plot_func,
                                                              dataset)
                )
                caption = (
                    f"Zonal mean profile of {dataset['long_name']} of dataset "
                    f"{dataset['dataset']} (project {dataset['project']})."
                )
            else:
                (plot_path, netcdf_paths) = (
                    self._plot_zonal_mean_profile_with_ref(plot_func, dataset,
                                                           ref_dataset)
                )
                caption = (
                    f"Zonal mean profile of {dataset['long_name']} of dataset "
                    f"{dataset['dataset']} (project {dataset['project']}) "
                    f"including bias relative to {ref_dataset['dataset']} "
                    f"(project {ref_dataset['project']})."
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
                'authors': ['schlund_manuel'],
                'caption': caption,
                'plot_types': ['vert'],
                'long_names': [dataset['long_name']],
            }
            with ProvenanceLogger(self.cfg) as provenance_logger:
                provenance_logger.log(plot_path, provenance_record)
                for netcdf_path in netcdf_paths:
                    provenance_logger.log(netcdf_path, provenance_record)

    def create_1d_profile_plot(self, datasets):
        """Create 1D profile plot."""
        plot_type = '1d_profile'
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        logger.info("Plotting %s", plot_type)
        fig = plt.figure(**self.cfg['figure_kwargs'])
        axes = fig.add_subplot()

        multi_dataset_facets = self._get_multi_dataset_facets(datasets)

        # Plot all datasets in one single figure
        ancestors = []
        cubes = {}
        for dataset in datasets:
            ancestors.append(dataset['filename'])
            cube = dataset['cube']
            cubes[self._get_label(dataset)] = cube
            self._check_cube_dimensions(cube, plot_type)

            # Plot 1D profile
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs['axes'] = axes

            iris.plot.plot(cube, **plot_kwargs)

        # Default plot appearance
        axes.set_title(multi_dataset_facets['long_name'])
        axes.set_xlabel(
            f"{multi_dataset_facets[self.cfg['group_variables_by']]} "
            f"[{multi_dataset_facets['units']}]"
        )
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
            n: datasets[0][n] for n in ('short_name', 'long_name', 'units')
        }
        io.save_1d_data(cubes, netcdf_path, z_coord.standard_name, var_attrs)

        # Provenance tracking
        caption = ("Vertical one-dimensional profile of "
                   f"{multi_dataset_facets['long_name']}"
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

    def create_variable_vs_lat_plot(self, datasets):
        """Create Variable as a function of latitude."""
        plot_type = 'variable_vs_lat'
        if plot_type not in self.plots:
            return
        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")
        logger.info("Plotting %s", plot_type)
        fig = plt.figure(**self.cfg['figure_kwargs'])
        axes = fig.add_subplot()

        # Plot all datasets in one single figure
        ancestors = []
        cubes = {}
        for dataset in datasets:
            ancestors.append(dataset['filename'])
            cube = dataset['cube']
            cubes[self._get_label(dataset)] = cube
            self._check_cube_dimensions(cube, plot_type)

            # Plot data
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs['axes'] = axes
            iris.plot.plot(cube, **plot_kwargs)

        # Default plot appearance
        multi_dataset_facets = self._get_multi_dataset_facets(datasets)
        axes.set_title(multi_dataset_facets['long_name'])
        axes.set_xlabel('latitude [°N]')
        axes.set_ylabel(
            f"{multi_dataset_facets[self.cfg['group_variables_by']]} "
            f"[{multi_dataset_facets['units']}]"
        )
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
            n: datasets[0][n] for n in ('short_name', 'long_name', 'units')
        }
        io.save_1d_data(cubes, netcdf_path, 'latitude', var_attrs)

        # Provenance tracking
        caption = (f"{multi_dataset_facets['long_name']} vs. latitude for "
                   f"various datasets.")
        provenance_record = {
            'ancestors': ancestors,
            'authors': ['sarauer_ellen'],
            'caption': caption,
            'plot_types': ['line'],
            'long_names': [var_attrs['long_name']],
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)
            provenance_logger.log(netcdf_path, provenance_record)

    def create_hovmoeller_z_vs_time_plot(self, datasets):
        """Create Hovmoeller Z vs. time plot."""
        plot_type = 'hovmoeller_z_vs_time'
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        # Get reference dataset if possible
        ref_dataset = self._get_reference_dataset(datasets)
        if ref_dataset is None:
            logger.info("Plotting %s without reference dataset", plot_type)
        else:
            logger.info("Plotting %s with reference dataset '%s'", plot_type,
                        self._get_label(ref_dataset))

        # Get plot function
        plot_func = self._get_plot_func(plot_type)

        # Create a single plot for each dataset (incl. reference dataset if
        # given)
        for dataset in datasets:
            if dataset == ref_dataset:
                continue
            ancestors = [dataset['filename']]
            if ref_dataset is None:
                (plot_path,
                 netcdf_paths) = (self._plot_hovmoeller_z_vs_time_without_ref(
                     plot_func, dataset))
                caption = (
                    f"Hovmoeller Z vs. time plot of {dataset['long_name']} "
                    f"of dataset "
                    f"{dataset['dataset']} (project {dataset['project']})."
                )
            else:
                (plot_path,
                 netcdf_paths) = (self._plot_hovmoeller_z_vs_time_with_ref(
                     plot_func, dataset, ref_dataset))
                caption = (
                    f"Hovmoeller Z vs. time plot of {dataset['long_name']} "
                    f"of dataset "
                    f"{dataset['dataset']} (project {dataset['project']}) "
                    f"including bias relative to {ref_dataset['dataset']} "
                    f"(project {ref_dataset['project']})."
                )
                ancestors.append(ref_dataset['filename'])

            # If statistics are shown add a brief description to the caption
            if self.plots[plot_type]['show_stats']:
                caption += (
                    " The number in the top left corner corresponds to the "
                    "spatiotemporal mean.")

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
                'authors': ['kuehbacher_birgit', 'heuer_helge'],
                'caption': caption,
                'plot_types': ['vert'],
                'long_names': [dataset['long_name']],
            }
            with ProvenanceLogger(self.cfg) as provenance_logger:
                provenance_logger.log(plot_path, provenance_record)
                for netcdf_path in netcdf_paths:
                    provenance_logger.log(netcdf_path, provenance_record)

    def create_hovmoeller_time_vs_lat_or_lon_plot(self, datasets):
        """Create the Hovmoeller plot with time vs latitude or longitude."""
        plot_type = 'hovmoeller_time_vs_lat_or_lon'
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        # Get reference dataset if possible
        ref_dataset = self._get_reference_dataset(datasets)
        if ref_dataset is None:
            logger.info("Plotting %s without reference dataset", plot_type)
        else:
            logger.info("Plotting %s with reference dataset '%s'", plot_type,
                        self._get_label(ref_dataset))

        # Get plot function
        plot_func = self._get_plot_func(plot_type)

        # Create a single plot for each dataset (incl. reference dataset if
        # given)
        for dataset in datasets:
            if dataset == ref_dataset:
                continue
            ancestors = [dataset['filename']]
            if ref_dataset is None:
                (plot_path, netcdf_paths) = (
                    self._plot_hovmoeller_time_vs_lat_or_lon_without_ref(
                        plot_func,
                        dataset)
                )
                caption = (
                    f"Hovmoeller plot of {dataset['long_name']} of dataset "
                    f"{dataset['dataset']} (project {dataset['project']})."
                )
            else:
                (plot_path, netcdf_paths) = (
                    self._plot_hovmoeller_time_vs_lat_or_lon_with_ref(
                        plot_func, dataset, ref_dataset)
                )
                caption = (
                    f"Hovmoeller plot of {dataset['long_name']} of dataset "
                    f"{dataset['dataset']} (project {dataset['project']}) "
                    f"including bias relative to {ref_dataset['dataset']} "
                    f"(project {ref_dataset['project']})."
                )
                ancestors.append(ref_dataset['filename'])

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
                'authors': ['schlund_manuel', 'kraft_jeremy', 'ruhe_lukas'],
                'caption': caption,
                'plot_types': ['zonal'],
                'long_names': [dataset['long_name']],
            }
            with ProvenanceLogger(self.cfg) as provenance_logger:
                provenance_logger.log(plot_path, provenance_record)
                for netcdf_path in netcdf_paths:
                    provenance_logger.log(netcdf_path, provenance_record)

    def compute(self):
        """Plot preprocessed data."""
        for (var_key, datasets) in self.grouped_input_data.items():
            logger.info("Processing variable %s", var_key)
            self.create_timeseries_plot(datasets)
            self.create_annual_cycle_plot(datasets)
            self.create_map_plot(datasets)
            self.create_zonal_mean_profile_plot(datasets)
            self.create_1d_profile_plot(datasets)
            self.create_variable_vs_lat_plot(datasets)
            self.create_hovmoeller_z_vs_time_plot(datasets)
            self.create_hovmoeller_time_vs_lat_or_lon_plot(datasets)


def main():
    """Run diagnostic."""
    with run_diagnostic() as config:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore',
                message="Using DEFAULT_SPHERICAL_EARTH_RADIUS",
                category=UserWarning,
                module='iris',
            )
            MultiDatasets(config).compute()


if __name__ == '__main__':
    main()
