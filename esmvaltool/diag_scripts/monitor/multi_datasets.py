#!/usr/bin/env python
"""Monitoring diagnostic to show multiple datasets in one plot (incl. biases).

Description
-----------
This diagnostic can be used to visualize multiple datasets in one plot.

For some plot types, a reference dataset can be defined. For this, use the
facet ``reference_for_monitor_diags: true`` in the definition of the dataset in
the recipe. Note that at most one reference dataset per variable is supported.

Please note that all benchmarking plot types (i.e. all plot types starting with
``benchmarking_``) require exactly one dataset (the dataset to be benchmarked)
to have the facet ``benchmark_dataset: true`` in the dataset entry of the
recipe.  For benchmarking line plots (i.e. ``benchmarking_annual_cycle``,
``benchmarking_diurnal_cycle``, ``benchmarking_timeseries``), it is recommended
to specify a particular line color and line style in the ``scripts`` section of
the recipe for the dataset to be benchmarked (``benchmark_dataset: true``) so
that this dataset is easy to identify in the plot.

Currently supported plot types (use the option ``plots`` to specify them):
    - Time series (plot type ``timeseries``): for each variable separately, all
      datasets are plotted in one single figure. Input data needs to be 1D with
      single dimension `time`.
    - Annual cycle (plot type ``annual_cycle``): for each variable separately,
      all datasets are plotted in one single figure. Input data needs to be 1D
      with single dimension `month_number`.
    - Diurnal cycle (plot type ``diurnal_cycle``): for each variable
      separately, all datasets are plotted in one single figure. Input data
      needs to be 1D with single dimension `hour`.
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
    - Hovmoeller annual cycle vs. latitude or longitude (plot type
      ``hovmoeller_anncyc_vs_lat_or_lon``): for each variable and dataset, an
      individual figure is plotted. If a reference dataset is defined, also
      include this dataset and a bias plot into the figure. Note that if a
      reference dataset is defined, all input datasets need to be given on the
      same temporal and horizontal grid (you can use the preprocessors
      :func:`esmvalcore.preprocessor.regrid_time` and
      :func:`esmvalcore.preprocessor.regrid` for this). Input data
      needs to be 2D with dimensions `month_number`, `latitude`/`longitude`.
    - Benchmarking plot annual cycles (``benchmarking_annual_cycle``):
      Same as plot type ``annual_cycle`` but including the range of metric
      results from an ensemble of models as shading.
    - Benchmarking box plots (``benchmarking_boxplot``):
      Box plots showing the metric results for given variables from a given
      model and the range from the first quartile to the third quartile, the
      median, and minimum and maximum values (excluding the outliers) from
      an ensemble of models for comparison.
    - Benchmarking plot diurnal cycles (``benchmarking_diurnal_cycle``):
      Same as plot type ``diurnal_cycle`` but including range of metric results
      from an ensemble of models as shading.
    - Benchmarking map plots (``benchmarking_map``):
      Same as plot type ``map`` but with stippled areas masking grid cells
      where the selected metric is smaller than the 90% percentile of
      corresponding values from an ensemble of models used for comparison.
    - Benchmarking plot time series (``benchmarking_timeseries``):
      Same as plot type ``timeseries`` but including the range of metric
      results from an ensemble of models as shading.
    - Benchmarking plot zonal mean profiles (plot type ``benchmarking_zonal``):
      Same as plot type ``zonal_mean_profile`` but with stippled areas masking
      grid cells where the selected metric is smaller than the 90% percentile
      of corresponding values from an ensemble of models used for comparison.

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
matplotlib_rc_params: dict, optional (default: {})
    Optional :class:`matplotlib.RcParams` used to customize matplotlib plots.
    Options given here will be passed to :func:`matplotlib.rc_context` and used
    for all plots produced with this diagnostic. Note: fontsizes specified here
    might be overwritten by the plot-type-specific option ``fontsize`` (see
    below).
plots: dict, optional
    Plot types plotted by this diagnostic (see list above). Dictionary keys
    must be ``timeseries``, ``annual_cycle``, ``map``, ``zonal_mean_profile``,
    ``1d_profile``, ``variable_vs_lat``, ``hovmoeller_z_vs_time``,
    ``hovmoeller_time_vs_lat_or_lon``, ``hovmoeller_anncyc_vs_lat_or_lon``.
    Dictionary values are dictionaries used
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
    replaced with the corresponding tags). ``{plot_dir}`` is replaced with the
    default ESMValTool plot directory (i.e.,
    ``output_dir/plots/diagnostic_name/script_name/``, see
    :ref:`esmvalcore:outputdata`).
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
hlines: list of dict, optional
    Horizontal lines to show in plot. Each list element corresponds to one
    line, and each list element should contain a dictionary with keywords
    arguments passed to :meth:`matplotlib.axes.Axes.axhline`. Example: ``[{y:
    0}, {y: 1, color: 'red'}]``.
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
    argument(s) for these functions (if values are dictionaries, these are
    interpreted as keyword arguments; otherwise a single argument is assumed).
    String arguments can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``. Examples: ``title: 'Awesome Plot of {long_name}'``, ``xlabel:
    '{short_name}'``, ``xlim: [0, 5]``.
time_format: str, optional (default: None)
    :func:`~datetime.datetime.strftime` format string that is used to format
    the time axis using :class:`matplotlib.dates.DateFormatter`. If ``None``,
    use the default formatting imposed by the iris plotting function.

Configuration options for plot type ``annual_cycle`` and ``diurnal_cycle``
--------------------------------------------------------------------------
gridline_kwargs: dict, optional
    Optional keyword arguments for grid lines. By default, ``color: lightgrey,
    alpha: 0.5`` are used. Use ``gridline_kwargs: false`` to not show grid
    lines.
hlines: list of dict, optional
    Horizontal lines to show in plot. Each list element corresponds to one
    line, and each list element should contain a dictionary with keywords
    arguments passed to :meth:`matplotlib.axes.Axes.axhline`. Example: ``[{y:
    0}, {y: 1, color: 'red'}]``.
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
    argument(s) for these functions (if values are dictionaries, these are
    interpreted as keyword arguments; otherwise a single argument is assumed).
    String arguments can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``. Examples: ``title: 'Awesome Plot of {long_name}'``, ``xlabel:
    '{short_name}'``, ``xlim: [0, 5]``.

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
fontsize: int, optional (default: None)
    Fontsize used for ticks, labels and titles. For the latter, use the given
    fontsize plus 2. Does not affect suptitles. If not given, use default
    matplotlib values. For a more fine-grained definition of fontsizes, use the
    option ``matplotlib_rc_params`` (see above).
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
    argument(s) for these functions (if values are dictionaries, these are
    interpreted as keyword arguments; otherwise a single argument is assumed).
    String arguments can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``. Examples: ``title: 'Awesome Plot of {long_name}'``, ``xlabel:
    '{short_name}'``, ``xlim: [0, 5]``.
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
fontsize: int, optional (default: None)
    Fontsize used for ticks, labels and titles. For the latter, use the given
    fontsize plus 2. Does not affect suptitles. If not given, use default
    matplotlib values. For a more fine-grained definition of fontsizes, use the
    option ``matplotlib_rc_params`` (see above).
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
    argument(s) for these functions (if values are dictionaries, these are
    interpreted as keyword arguments; otherwise a single argument is assumed).
    String arguments can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``. Examples: ``title: 'Awesome Plot of {long_name}'``, ``xlabel:
    '{short_name}'``, ``xlim: [0, 5]``.
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
hlines: list of dict, optional
    Horizontal lines to show in plot. Each list element corresponds to one
    line, and each list element should contain a dictionary with keywords
    arguments passed to :meth:`matplotlib.axes.Axes.axhline`. Example: ``[{y:
    0}, {y: 1, color: 'red'}]``.
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
    argument(s) for these functions (if values are dictionaries, these are
    interpreted as keyword arguments; otherwise a single argument is assumed).
    String arguments can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``. Examples: ``title: 'Awesome Plot of {long_name}'``, ``xlabel:
    '{short_name}'``, ``xlim: [0, 5]``.
show_y_minor_ticklabels: bool, optional (default: False)
    Show tick labels for the minor ticks on the Y axis.

Configuration options for plot type ``variable_vs_lat``
-------------------------------------------------------
gridline_kwargs: dict, optional
    Optional keyword arguments for grid lines. By default, ``color: lightgrey,
    alpha: 0.5`` are used. Use ``gridline_kwargs: false`` to not show grid
    lines.
hlines: list of dict, optional
    Horizontal lines to show in plot. Each list element corresponds to one
    line, and each list element should contain a dictionary with keywords
    arguments passed to :meth:`matplotlib.axes.Axes.axhline`. Example: ``[{y:
    0}, {y: 1, color: 'red'}]``.
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
    argument(s) for these functions (if values are dictionaries, these are
    interpreted as keyword arguments; otherwise a single argument is assumed).
    String arguments can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``. Examples: ``title: 'Awesome Plot of {long_name}'``, ``xlabel:
    '{short_name}'``, ``xlim: [0, 5]``.

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
fontsize: int, optional (default: None)
    Fontsize used for ticks, labels and titles. For the latter, use the given
    fontsize plus 2. Does not affect suptitles. If not given, use default
    matplotlib values. For a more fine-grained definition of fontsizes, use the
    option ``matplotlib_rc_params`` (see above).
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
    argument(s) for these functions (if values are dictionaries, these are
    interpreted as keyword arguments; otherwise a single argument is assumed).
    String arguments can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``. Examples: ``title: 'Awesome Plot of {long_name}'``, ``xlabel:
    '{short_name}'``, ``xlim: [0, 5]``.
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
fontsize: int, optional (default: None)
    Fontsize used for ticks, labels and titles. For the latter, use the given
    fontsize plus 2. Does not affect suptitles. If not given, use default
    matplotlib values. For a more fine-grained definition of fontsizes, use the
    option ``matplotlib_rc_params`` (see above).
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
    argument(s) for these functions (if values are dictionaries, these are
    interpreted as keyword arguments; otherwise a single argument is assumed).
    String arguments can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``. Examples: ``title: 'Awesome Plot of {long_name}'``, ``xlabel:
    '{short_name}'``, ``xlim: [0, 5]``.
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
time_on: str, optional (default: y-axis)
    Optional switch to change the orientation of the plot so that time is on
    the x-axis ``time_on: x-axis``. Default orientation is time on y-axis and
    lat/lon on x-axis.

Configuration options for plot type ``hovmoeller_anncyc_vs_lat_or_lon``
-----------------------------------------------------------------------
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
fontsize: int, optional (default: None)
    Fontsize used for ticks, labels and titles. For the latter, use the given
    fontsize plus 2. Does not affect suptitles. If not given, use default
    matplotlib values. For a more fine-grained definition of fontsizes, use the
    option ``matplotlib_rc_params`` (see above).
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
    argument(s) for these functions (if values are dictionaries, these are
    interpreted as keyword arguments; otherwise a single argument is assumed).
    String arguments can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``. Examples: ``title: 'Awesome Plot of {long_name}'``, ``xlabel:
    '{short_name}'``, ``xlim: [0, 5]``.
rasterize: bool, optional (default: True)
    If ``True``, use rasterization_ for profile plots to produce smaller files.
    This is only relevant for vector graphics (e.g., ``output_file_type:
    pdf,svg,ps``).
show_y_minor_ticks: bool, optional (default: True)
    Show minor ticks for time on the Y axis.
show_x_minor_ticks: bool, optional (default: True)
    Show minor ticks for latitude or longitude on the X axis.
time_on: str, optional (default: y-axis)
    Optional switch to change the orientation of the plot so that time is on
    the x-axis ``time_on: x-axis``. Default orientation is time on y-axis and
    lat/lon on x-axis.

Configuration options for plot type ``benchmarking_annual_cycle``
-----------------------------------------------------------------
Same as for plot type ``annual_cycle``.

Configuration options for plot type ``benchmarking_boxplot``
------------------------------------------------------------
fontsize: int, optional (default: None)
    Fontsize used for ticks, labels and titles. For the latter, use the given
    fontsize plus 2. Does not affect suptitles. If not given, use default
    matplotlib values. For a more fine-grained definition of fontsizes, use the
    option ``matplotlib_rc_params`` (see above).
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
pyplot_kwargs: dict, optional
    Optional calls to functions of :mod:`matplotlib.pyplot`. Dictionary keys
    are functions of :mod:`matplotlib.pyplot`. Dictionary values are used as
    argument(s) for these functions (if values are dictionaries, these are
    interpreted as keyword arguments; otherwise a single argument is assumed).
    String arguments can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``. Examples: ``title: 'Awesome Plot of {long_name}'``, ``xlabel:
    '{short_name}'``, ``xlim: [0, 5]``.
var_order: list of str, optional
    Optional list of strings containing variable names to define the order of
    the variables plotted.

Configuration options for plot type ``benchmarking_diurnal_cycle``
------------------------------------------------------------------
Same as for plot type ``diurnal_cycle``.

Configuration options for plot type ``benchmarking_map``
--------------------------------------------------------
cbar_label: str, optional (default: '{short_name} [{units}]')
    Colorbar label. Can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``.
cbar_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.colorbar`. By
    default, uses ``orientation: horizontal, aspect: 30``.
fontsize: int, optional (default: None)
    Fontsize used for ticks, labels and titles. For the latter, use the given
    fontsize plus 2. Does not affect suptitles. If not given, use default
    matplotlib values. For a more fine-grained definition of fontsizes, use the
    option ``matplotlib_rc_params`` (see above).
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
    argument(s) for these functions (if values are dictionaries, these are
    interpreted as keyword arguments; otherwise a single argument is assumed).
    String arguments can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``. Examples: ``title: 'Awesome Plot of {long_name}'``, ``xlabel:
    '{short_name}'``, ``xlim: [0, 5]``.
rasterize: bool, optional (default: True)
    If ``True``, use rasterization_ for map plots to produce smaller files.
    This is only relevant for vector graphics (e.g., ``output_file_type:
    pdf,svg,ps``).

Configuration options for plot type ``benchmarking_timeseries``
---------------------------------------------------------------
Same as for plot type ``timeseries``.

Configuration options for plot type ``benchmarking_zonal``
----------------------------------------------------------
cbar_label: str, optional (default: '{short_name} [{units}]')
    Colorbar label. Can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``.
cbar_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.colorbar`. By
    default, uses ``orientation: vertical``.
fontsize: int, optional (default: None)
    Fontsize used for ticks, labels and titles. For the latter, use the given
    fontsize plus 2. Does not affect suptitles. If not given, use default
    matplotlib values. For a more fine-grained definition of fontsizes, use the
    option ``matplotlib_rc_params`` (see above).
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
pyplot_kwargs: dict, optional
    Optional calls to functions of :mod:`matplotlib.pyplot`. Dictionary keys
    are functions of :mod:`matplotlib.pyplot`. Dictionary values are used as
    argument(s) for these functions (if values are dictionaries, these are
    interpreted as keyword arguments; otherwise a single argument is assumed).
    String arguments can include facets in curly brackets which will be derived
    from the corresponding dataset, e.g., ``{project}``, ``{short_name}``,
    ``{exp}``. Examples: ``title: 'Awesome Plot of {long_name}'``, ``xlabel:
    '{short_name}'``, ``xlim: [0, 5]``.
rasterize: bool, optional (default: True)
    If ``True``, use rasterization_ for profile plots to produce smaller files.
    This is only relevant for vector graphics (e.g., ``output_file_type:
    pdf,svg,ps``).
show_y_minor_ticklabels: bool, optional (default: False)
    Show tick labels for the minor ticks on the Y axis.


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
import dask.array as da
import iris
import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from iris.analysis.cartography import area_weights
from iris.coord_categorisation import add_year
from iris.coords import AuxCoord
from iris.exceptions import ConstraintMismatchError
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

        if "profile" in self.plots:
            logger.warning(
                "The plot_type ``profile`` for zonal mean profiles"
                " has been deprecated in ESMValTool version 2.9.0"
                " and is scheduled for removal in version 2.11.0."
                " Please use plot type ``zonal_mean_profile``"
                " instead. This is an exact replacement."
            )
            if "zonal_mean_profile" in self.plots:
                raise ValueError(
                    "Both ``profile`` and ``zonal_mean_profile`` is used."
                    " Please use ``zonal_mean_profile`` only."
                )
            self.plots["zonal_mean_profile"] = self.plots.pop("profile")

        # Check given plot types and set default settings for them
        self.supported_plot_types = [
            "timeseries",
            "annual_cycle",
            "diurnal_cycle",
            "map",
            "zonal_mean_profile",
            "1d_profile",
            "variable_vs_lat",
            "hovmoeller_z_vs_time",
            "hovmoeller_time_vs_lat_or_lon",
            "hovmoeller_anncyc_vs_lat_or_lon",
            "benchmarking_annual_cycle",
            "benchmarking_boxplot",
            "benchmarking_diurnal_cycle",
            "benchmarking_map",
            "benchmarking_timeseries",
            "benchmarking_zonal",
        ]
        for plot_type, plot_options in self.plots.items():
            if plot_type not in self.supported_plot_types:
                raise ValueError(
                    f"Got unexpected plot type '{plot_type}' for option "
                    f"'plots', expected one of {self.supported_plot_types}"
                )
            if plot_options is None:
                self.plots[plot_type] = {}

            # Default options for the different plot types
            if plot_type == "timeseries":
                self.plots[plot_type].setdefault("annual_mean_kwargs", {})
                self.plots[plot_type].setdefault("gridline_kwargs", {})
                self.plots[plot_type].setdefault("hlines", [])
                self.plots[plot_type].setdefault("legend_kwargs", {})
                self.plots[plot_type].setdefault("plot_kwargs", {})
                self.plots[plot_type].setdefault("pyplot_kwargs", {})
                self.plots[plot_type].setdefault("time_format", None)

            elif plot_type == "benchmarking_timeseries":
                self.plots[plot_type].setdefault("annual_mean_kwargs", {})
                self.plots[plot_type].setdefault("gridline_kwargs", {})
                self.plots[plot_type].setdefault("hlines", [])
                self.plots[plot_type].setdefault("legend_kwargs", {})
                self.plots[plot_type].setdefault("plot_kwargs", {})
                self.plots[plot_type].setdefault("pyplot_kwargs", {})
                self.plots[plot_type].setdefault("time_format", None)

            elif plot_type == "annual_cycle":
                self.plots[plot_type].setdefault("gridline_kwargs", {})
                self.plots[plot_type].setdefault("hlines", [])
                self.plots[plot_type].setdefault("legend_kwargs", {})
                self.plots[plot_type].setdefault("plot_kwargs", {})
                self.plots[plot_type].setdefault("pyplot_kwargs", {})

            elif plot_type == "benchmarking_annual_cycle":
                self.plots[plot_type].setdefault("gridline_kwargs", {})
                self.plots[plot_type].setdefault("hlines", [])
                self.plots[plot_type].setdefault("legend_kwargs", {})
                self.plots[plot_type].setdefault("plot_kwargs", {})
                self.plots[plot_type].setdefault("pyplot_kwargs", {})

            elif plot_type == "diurnal_cycle":
                self.plots[plot_type].setdefault("gridline_kwargs", {})
                self.plots[plot_type].setdefault("hlines", [])
                self.plots[plot_type].setdefault("legend_kwargs", {})
                self.plots[plot_type].setdefault("plot_kwargs", {})
                self.plots[plot_type].setdefault("pyplot_kwargs", {})

            elif plot_type == "benchmarking_diurnal_cycle":
                self.plots[plot_type].setdefault("gridline_kwargs", {})
                self.plots[plot_type].setdefault("hlines", [])
                self.plots[plot_type].setdefault("legend_kwargs", {})
                self.plots[plot_type].setdefault("plot_kwargs", {})
                self.plots[plot_type].setdefault("pyplot_kwargs", {})

            elif plot_type == "benchmarking_boxplot":
                self.plots[plot_type].setdefault("fontsize", None)
                self.plots[plot_type].setdefault("plot_kwargs", {})
                self.plots[plot_type].setdefault("pyplot_kwargs", {})
                self.plots[plot_type].setdefault("var_order", None)

            elif plot_type == "map":
                self.plots[plot_type].setdefault(
                    "cbar_label", "{short_name} [{units}]"
                )
                self.plots[plot_type].setdefault(
                    "cbar_label_bias", "Δ{short_name} [{units}]"
                )
                self.plots[plot_type].setdefault(
                    "cbar_kwargs", {"orientation": "horizontal", "aspect": 30}
                )
                self.plots[plot_type].setdefault("cbar_kwargs_bias", {})
                self.plots[plot_type].setdefault("common_cbar", False)
                self.plots[plot_type].setdefault("fontsize", None)
                self.plots[plot_type].setdefault("gridline_kwargs", {})
                self.plots[plot_type].setdefault("plot_func", "contourf")
                self.plots[plot_type].setdefault("plot_kwargs", {})
                self.plots[plot_type].setdefault("plot_kwargs_bias", {})
                self.plots[plot_type]["plot_kwargs_bias"].setdefault(
                    "cmap", "bwr"
                )
                self.plots[plot_type]["plot_kwargs_bias"].setdefault(
                    "norm", "centered"
                )
                if "projection" not in self.plots[plot_type]:
                    self.plots[plot_type].setdefault("projection", "Robinson")
                    self.plots[plot_type].setdefault(
                        "projection_kwargs", {"central_longitude": 10}
                    )
                else:
                    self.plots[plot_type].setdefault("projection_kwargs", {})
                self.plots[plot_type].setdefault("pyplot_kwargs", {})
                self.plots[plot_type].setdefault("rasterize", True)
                self.plots[plot_type].setdefault("show_stats", True)
                self.plots[plot_type].setdefault("x_pos_stats_avg", 0.0)
                self.plots[plot_type].setdefault("x_pos_stats_bias", 0.92)

            elif plot_type == "benchmarking_map":
                self.plots[plot_type].setdefault(
                    "cbar_label", "{short_name} [{units}]"
                )
                self.plots[plot_type].setdefault(
                    "cbar_kwargs", {"orientation": "horizontal", "aspect": 30}
                )
                self.plots[plot_type].setdefault("fontsize", None)
                self.plots[plot_type].setdefault("plot_func", "contourf")
                self.plots[plot_type].setdefault("plot_kwargs", {})
                if "projection" not in self.plots[plot_type]:
                    self.plots[plot_type].setdefault("projection", "Robinson")
                    self.plots[plot_type].setdefault(
                        "projection_kwargs", {"central_longitude": 10}
                    )
                else:
                    self.plots[plot_type].setdefault("projection_kwargs", {})
                self.plots[plot_type].setdefault("pyplot_kwargs", {})
                self.plots[plot_type].setdefault("rasterize", True)

            elif plot_type == "zonal_mean_profile":
                self.plots[plot_type].setdefault(
                    "cbar_label", "{short_name} [{units}]"
                )
                self.plots[plot_type].setdefault(
                    "cbar_label_bias", "Δ{short_name} [{units}]"
                )
                self.plots[plot_type].setdefault(
                    "cbar_kwargs", {"orientation": "vertical"}
                )
                self.plots[plot_type].setdefault("cbar_kwargs_bias", {})
                self.plots[plot_type].setdefault("common_cbar", False)
                self.plots[plot_type].setdefault("fontsize", None)
                self.plots[plot_type].setdefault("log_y", True)
                self.plots[plot_type].setdefault("plot_func", "contourf")
                self.plots[plot_type].setdefault("plot_kwargs", {})
                self.plots[plot_type].setdefault("plot_kwargs_bias", {})
                self.plots[plot_type]["plot_kwargs_bias"].setdefault(
                    "cmap", "bwr"
                )
                self.plots[plot_type]["plot_kwargs_bias"].setdefault(
                    "norm", "centered"
                )
                self.plots[plot_type].setdefault("pyplot_kwargs", {})
                self.plots[plot_type].setdefault("rasterize", True)
                self.plots[plot_type].setdefault("show_stats", True)
                self.plots[plot_type].setdefault(
                    "show_y_minor_ticklabels", False
                )
                self.plots[plot_type].setdefault("x_pos_stats_avg", 0.01)
                self.plots[plot_type].setdefault("x_pos_stats_bias", 0.7)

            elif plot_type == "benchmarking_zonal":
                self.plots[plot_type].setdefault(
                    "cbar_label", "{short_name} [{units}]"
                )
                self.plots[plot_type].setdefault(
                    "cbar_kwargs", {"orientation": "vertical"}
                )
                self.plots[plot_type].setdefault("fontsize", None)
                self.plots[plot_type].setdefault("log_y", True)
                self.plots[plot_type].setdefault("plot_func", "contourf")
                self.plots[plot_type].setdefault("plot_kwargs", {})
                self.plots[plot_type].setdefault("pyplot_kwargs", {})
                self.plots[plot_type].setdefault("rasterize", True)
                self.plots[plot_type].setdefault(
                    "show_y_minor_ticklabels", False
                )

            elif plot_type == "1d_profile":
                self.plots[plot_type].setdefault("aspect_ratio", 1.5)
                self.plots[plot_type].setdefault("gridline_kwargs", {})
                self.plots[plot_type].setdefault("hlines", [])
                self.plots[plot_type].setdefault("legend_kwargs", {})
                self.plots[plot_type].setdefault("log_x", False)
                self.plots[plot_type].setdefault("log_y", True)
                self.plots[plot_type].setdefault("plot_kwargs", {})
                self.plots[plot_type].setdefault("pyplot_kwargs", {})
                self.plots[plot_type].setdefault(
                    "show_y_minor_ticklabels", False
                )

            elif plot_type == "variable_vs_lat":
                self.plots[plot_type].setdefault("gridline_kwargs", {})
                self.plots[plot_type].setdefault("hlines", [])
                self.plots[plot_type].setdefault("legend_kwargs", {})
                self.plots[plot_type].setdefault("plot_kwargs", {})
                self.plots[plot_type].setdefault("pyplot_kwargs", {})

            elif plot_type == "hovmoeller_z_vs_time":
                self.plots[plot_type].setdefault(
                    "cbar_label", "{short_name} [{units}]"
                )
                self.plots[plot_type].setdefault(
                    "cbar_label_bias", "Δ{short_name} [{units}]"
                )
                self.plots[plot_type].setdefault(
                    "cbar_kwargs", {"orientation": "vertical"}
                )
                self.plots[plot_type].setdefault("cbar_kwargs_bias", {})
                self.plots[plot_type].setdefault("common_cbar", False)
                self.plots[plot_type].setdefault("fontsize", None)
                self.plots[plot_type].setdefault("log_y", True)
                self.plots[plot_type].setdefault("plot_func", "contourf")
                self.plots[plot_type].setdefault("plot_kwargs", {})
                self.plots[plot_type].setdefault("plot_kwargs_bias", {})
                self.plots[plot_type]["plot_kwargs_bias"].setdefault(
                    "cmap", "bwr"
                )
                self.plots[plot_type]["plot_kwargs_bias"].setdefault(
                    "norm", "centered"
                )
                self.plots[plot_type].setdefault("pyplot_kwargs", {})
                self.plots[plot_type].setdefault("rasterize", True)
                self.plots[plot_type].setdefault("show_stats", True)
                self.plots[plot_type].setdefault(
                    "show_y_minor_ticklabels", False
                )
                self.plots[plot_type].setdefault("time_format", None)
                self.plots[plot_type].setdefault("x_pos_stats_avg", 0.01)
                self.plots[plot_type].setdefault("x_pos_stats_bias", 0.7)

            elif plot_type == "hovmoeller_time_vs_lat_or_lon":
                self.plots[plot_type].setdefault(
                    "cbar_label", "{short_name} [{units}]"
                )
                self.plots[plot_type].setdefault(
                    "cbar_label_bias", "Δ{short_name} [{units}]"
                )
                self.plots[plot_type].setdefault(
                    "cbar_kwargs", {"orientation": "vertical"}
                )
                self.plots[plot_type].setdefault("cbar_kwargs_bias", {})
                self.plots[plot_type].setdefault("common_cbar", False)
                self.plots[plot_type].setdefault("fontsize", None)
                self.plots[plot_type].setdefault("plot_func", "contourf")
                self.plots[plot_type].setdefault("plot_kwargs", {})
                self.plots[plot_type].setdefault("plot_kwargs_bias", {})
                self.plots[plot_type]["plot_kwargs_bias"].setdefault(
                    "cmap", "bwr"
                )
                self.plots[plot_type]["plot_kwargs_bias"].setdefault(
                    "norm", "centered"
                )
                self.plots[plot_type].setdefault("pyplot_kwargs", {})
                self.plots[plot_type].setdefault("rasterize", True)
                self.plots[plot_type].setdefault("show_y_minor_ticks", True)
                self.plots[plot_type].setdefault("show_x_minor_ticks", True)
                self.plots[plot_type].setdefault("time_format", None)
                self.plots[plot_type].setdefault("time_on", "y-axis")

            elif plot_type == "hovmoeller_anncyc_vs_lat_or_lon":
                self.plots[plot_type].setdefault(
                    "cbar_label", "{short_name} [{units}]"
                )
                self.plots[plot_type].setdefault(
                    "cbar_label_bias", "Δ{short_name} [{units}]"
                )
                self.plots[plot_type].setdefault(
                    "cbar_kwargs", {"orientation": "vertical"}
                )
                self.plots[plot_type].setdefault("cbar_kwargs_bias", {})
                self.plots[plot_type].setdefault("common_cbar", False)
                self.plots[plot_type].setdefault("fontsize", None)
                self.plots[plot_type].setdefault("plot_func", "contourf")
                self.plots[plot_type].setdefault("plot_kwargs", {})
                self.plots[plot_type].setdefault("plot_kwargs_bias", {})
                self.plots[plot_type]["plot_kwargs_bias"].setdefault(
                    "cmap", "bwr"
                )
                self.plots[plot_type]["plot_kwargs_bias"].setdefault(
                    "norm", "centered"
                )
                self.plots[plot_type].setdefault("pyplot_kwargs", {})
                self.plots[plot_type].setdefault("rasterize", True)
                self.plots[plot_type].setdefault("show_y_minor_ticks", True)
                self.plots[plot_type].setdefault("show_x_minor_ticks", True)
                self.plots[plot_type].setdefault("time_on", "y-axis")

        # Check that facet_used_for_labels is present for every dataset
        for dataset in self.input_data:
            if self.cfg["facet_used_for_labels"] not in dataset:
                raise ValueError(
                    f"facet_used_for_labels "
                    f"'{self.cfg['facet_used_for_labels']}' not present for "
                    f"the following dataset:\n{pformat(dataset)}"
                )

        # Load seaborn settings
        sns.set_theme(**self.cfg["seaborn_settings"])

    def _add_colorbar(
        self,
        plot_type,
        plot_left,
        plot_right,
        axes_left,
        axes_right,
        dataset_left,
        dataset_right,
    ):
        """Add colorbar(s) for plots."""
        fontsize = (
            self.plots[plot_type]["fontsize"] or mpl.rcParams["axes.labelsize"]
        )
        cbar_kwargs = self._get_cbar_kwargs(plot_type)
        cbar_label_left = self._get_cbar_label(plot_type, dataset_left)
        cbar_label_right = self._get_cbar_label(plot_type, dataset_right)

        # Create one common colorbar for the top panels
        # Note: Increase aspect ratio for nicer looks
        if self.plots[plot_type]["common_cbar"]:
            if "aspect" in cbar_kwargs:
                cbar_kwargs["aspect"] += 20.0
            cbar = plt.colorbar(
                plot_left, ax=[axes_left, axes_right], **cbar_kwargs
            )
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

    def _add_stats(
        self, plot_type, axes, dim_coords, dataset, ref_dataset=None
    ):
        """Add text to plot that describes basic statistics."""
        if not self.plots[plot_type]["show_stats"]:
            return

        # Extract cube(s)
        cube = dataset["cube"]
        if ref_dataset is None:
            ref_cube = None
            label = self._get_label(dataset)
        else:
            ref_cube = ref_dataset["cube"]
            label = (
                f"{self._get_label(dataset)} vs. "
                f"{self._get_label(ref_dataset)}"
            )

        # Different options for the different plots types
        fontsize = 6.0
        y_pos = 0.95
        if all(
            [
                "x_pos_stats_avg" in self.plots[plot_type],
                "x_pos_stats_bias" in self.plots[plot_type],
            ]
        ):
            x_pos_bias = self.plots[plot_type]["x_pos_stats_bias"]
            x_pos = self.plots[plot_type]["x_pos_stats_avg"]
        else:
            raise NotImplementedError(f"plot_type '{plot_type}' not supported")

        # For zonal_mean_profile plots add scalar longitude coordinate
        # (necessary for calculation of area weights). The exact values for the
        # points/bounds of this coordinate do not matter since they don't
        # change the weights.
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

        # Mean
        weights = area_weights(cube)
        if ref_cube is None:
            mean = cube.collapsed(
                dim_coords, iris.analysis.MEAN, weights=weights
            )
            logger.info(
                "Area-weighted mean of %s for %s = %f%s",
                dataset["short_name"],
                label,
                mean.data,
                dataset["units"],
            )
        else:
            mean = (cube - ref_cube).collapsed(
                dim_coords, iris.analysis.MEAN, weights=weights
            )
            logger.info(
                "Area-weighted bias of %s for %s = %f%s",
                dataset["short_name"],
                label,
                mean.data,
                dataset["units"],
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
        rmse = (cube - ref_cube).collapsed(
            dim_coords, iris.analysis.RMS, weights=weights
        )
        if np.abs(rmse.data) >= 0.1:
            rmse_val = f"{rmse.data:.2f} {cube.units}"
        else:
            rmse_val = f"{rmse.data:.2e} {cube.units}"
        axes.text(
            x_pos_bias,
            y_pos,
            f"RMSE={rmse_val}",
            fontsize=fontsize,
            transform=axes.transAxes,
        )
        logger.info(
            "Area-weighted RMSE of %s for %s = %f%s",
            dataset["short_name"],
            label,
            rmse.data,
            dataset["units"],
        )

        # Weighted R2
        mask = np.ma.getmaskarray(cube.data).ravel()
        mask |= np.ma.getmaskarray(ref_cube.data).ravel()
        cube_data = cube.data.ravel()[~mask]
        ref_cube_data = ref_cube.data.ravel()[~mask]
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
            dataset["short_name"],
            label,
            r2_val,
        )

    def _get_custom_mpl_rc_params(self, plot_type):
        """Get custom matplotlib rcParams."""
        custom_rc_params = {}
        fontsize = self.plots[plot_type]["fontsize"]
        if fontsize is not None:
            custom_rc_params.update(
                {
                    "axes.titlesize": fontsize + 2.0,
                    "axes.labelsize": fontsize,
                    "xtick.labelsize": fontsize,
                    "ytick.labelsize": fontsize,
                }
            )
        return custom_rc_params

    def _get_label(self, dataset):
        """Get label of dataset."""
        return dataset[self.cfg["facet_used_for_labels"]]

    def _get_cbar_kwargs(self, plot_type, bias=False):
        """Get colorbar kwargs."""
        cbar_kwargs = deepcopy(self.plots[plot_type]["cbar_kwargs"])
        if bias:
            cbar_kwargs.update(self.plots[plot_type]["cbar_kwargs_bias"])
        return deepcopy(cbar_kwargs)

    def _get_cbar_label(self, plot_type, dataset, bias=False):
        """Get colorbar label."""
        if bias:
            cbar_label = self.plots[plot_type]["cbar_label_bias"]
            descr = f"cbar_label_bias of {plot_type} '{cbar_label}'"
        else:
            cbar_label = self.plots[plot_type]["cbar_label"]
            descr = f"cbar_label of {plot_type} '{cbar_label}'"
        cbar_label = self._fill_facet_placeholders(cbar_label, dataset, descr)
        return cbar_label

    def _get_gridline_kwargs(self, plot_type):
        """Get gridline kwargs."""
        gridline_kwargs = self.plots[plot_type]["gridline_kwargs"]
        return deepcopy(gridline_kwargs)

    def _get_map_projection(self):
        """Get projection used for map plots."""
        plot_type = "map"
        projection = self.plots[plot_type]["projection"]
        projection_kwargs = self.plots[plot_type]["projection_kwargs"]

        # Check if desired projection is valid
        if not hasattr(ccrs, projection):
            raise AttributeError(
                f"Got invalid projection '{projection}' for plotting "
                f"{plot_type}, expected class of cartopy.crs"
            )

        return getattr(ccrs, projection)(**projection_kwargs)

    def _get_benchmarking_projection(self):
        """Get projection used for benchmarking map plots."""
        plot_type = "benchmarking_map"
        projection = self.plots[plot_type]["projection"]
        projection_kwargs = self.plots[plot_type]["projection_kwargs"]

        # Check if desired projection is valid
        if not hasattr(ccrs, projection):
            raise AttributeError(
                f"Got invalid projection '{projection}' for plotting "
                f"{plot_type}, expected class of cartopy.crs"
            )

        return getattr(ccrs, projection)(**projection_kwargs)

    def _get_plot_func(self, plot_type):
        """Get plot function."""
        plot_func = self.plots[plot_type]["plot_func"]
        if not hasattr(iris.plot, plot_func):
            raise AttributeError(
                f"Got invalid plot function '{plot_func}' for plotting "
                f"{plot_type}, expected function of iris.plot"
            )
        logger.info(
            "Creating %s plots using function '%s'", plot_type, plot_func
        )
        return getattr(iris.plot, plot_func)

    def _get_plot_kwargs(self, plot_type, dataset, bias=False):
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

        # Default settings for different plot types
        if plot_type in (
            "timeseries",
            "annual_cycle",
            "benchmarking_annual_cycle",
            "1d_profile",
            "diurnal_cycle",
            "benchmarking_diurnal_cycle",
            "variable_vs_lat",
            "benchmarking_timeseries",
        ):
            plot_kwargs.setdefault("label", label)

        if plot_kwargs.get("norm") == "centered":
            norm = CenteredNorm(
                vcenter=plot_kwargs.pop("vcenter", 0.0),
                halfrange=plot_kwargs.pop("halfrange", None),
            )
            plot_kwargs["norm"] = norm

        return deepcopy(plot_kwargs)

    def _load_and_preprocess_data(self):
        """Load and preprocess data."""
        input_data = list(self.cfg["input_data"].values())

        for dataset in input_data:
            filename = dataset["filename"]
            logger.info("Loading %s", filename)
            cubes = iris.load(filename)
            if len(cubes) == 1:
                cube = cubes[0]
            else:
                var_name = dataset["short_name"]
                try:
                    cube = cubes.extract_cube(
                        iris.NameConstraint(var_name=var_name)
                    )
                except ConstraintMismatchError as exc:
                    var_names = [c.var_name for c in cubes]
                    raise ValueError(
                        f"Cannot load data: multiple variables ({var_names}) "
                        f"are available in file {filename}, but not the "
                        f"requested '{var_name}'"
                    ) from exc

            # Fix time coordinate if present
            if cube.coords("time", dim_coords=True):
                ih.unify_time_coord(cube)

            # Fix Z-coordinate if present
            if cube.coords("air_pressure", dim_coords=True):
                z_coord = cube.coord("air_pressure", dim_coords=True)
                z_coord.attributes["positive"] = "down"
                z_coord.convert_units("hPa")
            elif cube.coords("altitude", dim_coords=True):
                z_coord = cube.coord("altitude")
                z_coord.attributes["positive"] = "up"

            dataset["cube"] = cube

        return input_data

    def _plot_map_with_ref(self, plot_func, dataset, ref_dataset):
        """Plot map plot for single dataset with a reference dataset."""
        plot_type = "map"
        logger.info(
            "Plotting map with reference dataset '%s' for '%s'",
            self._get_label(ref_dataset),
            self._get_label(dataset),
        )

        # Make sure that the data has the correct dimensions
        cube = dataset["cube"]
        ref_cube = ref_dataset["cube"]
        dim_coords_dat = self._check_cube_dimensions(cube, plot_type)
        dim_coords_ref = self._check_cube_dimensions(ref_cube, plot_type)

        # Create single figure with multiple axes
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg["figure_kwargs"])
            gridspec = GridSpec(
                5, 4, figure=fig, height_ratios=[1.0, 1.0, 0.4, 1.0, 1.0]
            )

            # Options used for all subplots
            projection = self._get_map_projection()
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            gridline_kwargs = self._get_gridline_kwargs(plot_type)
            fontsize = (
                self.plots[plot_type]["fontsize"]
                or mpl.rcParams["axes.labelsize"]
            )

            # Plot dataset (top left)
            axes_data = fig.add_subplot(
                gridspec[0:2, 0:2], projection=projection
            )
            plot_kwargs["axes"] = axes_data
            if plot_func is iris.plot.contourf:
                # see https://github.com/SciTools/cartopy/issues/2457
                # and https://github.com/SciTools/cartopy/issues/2468
                plot_kwargs["transform_first"] = True
                npx = da if cube.has_lazy_data() else np
                cube_to_plot = cube.copy(
                    npx.ma.filled(cube.core_data(), np.nan)
                )
            else:
                cube_to_plot = cube
            plot_data = plot_func(cube_to_plot, **plot_kwargs)
            axes_data.coastlines()
            if gridline_kwargs is not False:
                axes_data.gridlines(**gridline_kwargs)
            axes_data.set_title(self._get_label(dataset), pad=3.0)
            self._add_stats(plot_type, axes_data, dim_coords_dat, dataset)
            self._process_pyplot_kwargs(plot_type, dataset)

            # Plot reference dataset (top right)
            # Note: make sure to use the same vmin and vmax than the top left
            # plot if a common cpltolorbar is desired
            axes_ref = fig.add_subplot(
                gridspec[0:2, 2:4], projection=projection
            )
            plot_kwargs["axes"] = axes_ref
            if self.plots[plot_type]["common_cbar"]:
                plot_kwargs.setdefault("vmin", plot_data.get_clim()[0])
                plot_kwargs.setdefault("vmax", plot_data.get_clim()[1])
            if plot_func is iris.plot.contourf:
                # see https://github.com/SciTools/cartopy/issues/2457
                # and https://github.com/SciTools/cartopy/issues/2468
                plot_kwargs["transform_first"] = True
                npx = da if ref_cube.has_lazy_data() else np
                ref_cube_to_plot = ref_cube.copy(
                    npx.ma.filled(ref_cube.core_data(), np.nan)
                )
            else:
                ref_cube_to_plot = ref_cube
            plot_ref = plot_func(ref_cube_to_plot, **plot_kwargs)
            axes_ref.coastlines()
            if gridline_kwargs is not False:
                axes_ref.gridlines(**gridline_kwargs)
            axes_ref.set_title(self._get_label(ref_dataset), pad=3.0)
            self._add_stats(plot_type, axes_ref, dim_coords_ref, ref_dataset)
            self._process_pyplot_kwargs(plot_type, ref_dataset)

            # Add colorbar(s)
            self._add_colorbar(
                plot_type,
                plot_data,
                plot_ref,
                axes_data,
                axes_ref,
                dataset,
                ref_dataset,
            )

            # Plot bias (bottom center)
            bias_cube = cube - ref_cube
            axes_bias = fig.add_subplot(
                gridspec[3:5, 1:3], projection=projection
            )
            plot_kwargs_bias = self._get_plot_kwargs(
                plot_type, dataset, bias=True
            )
            plot_kwargs_bias["axes"] = axes_bias
            if plot_func is iris.plot.contourf:
                # see https://github.com/SciTools/cartopy/issues/2457
                # and https://github.com/SciTools/cartopy/issues/2468
                plot_kwargs_bias["transform_first"] = True
                npx = da if bias_cube.has_lazy_data() else np
                bias_cube_to_plot = bias_cube.copy(
                    npx.ma.filled(bias_cube.core_data(), np.nan)
                )
            else:
                bias_cube_to_plot = bias_cube
            plot_bias = plot_func(bias_cube_to_plot, **plot_kwargs_bias)
            axes_bias.coastlines()
            if gridline_kwargs is not False:
                axes_bias.gridlines(**gridline_kwargs)
            axes_bias.set_title(
                f"{self._get_label(dataset)} - {self._get_label(ref_dataset)}",
                pad=3.0,
            )
            cbar_kwargs_bias = self._get_cbar_kwargs(plot_type, bias=True)
            cbar_bias = fig.colorbar(
                plot_bias, ax=axes_bias, **cbar_kwargs_bias
            )
            cbar_bias.set_label(
                self._get_cbar_label(plot_type, dataset, bias=True),
                fontsize=fontsize,
            )
            cbar_bias.ax.tick_params(labelsize=fontsize)
            self._add_stats(
                plot_type, axes_bias, dim_coords_dat, dataset, ref_dataset
            )

            # Customize plot
            fig.suptitle(dataset["long_name"])
            self._process_pyplot_kwargs(plot_type, dataset)

            # Rasterization
            if self.plots[plot_type]["rasterize"]:
                self._set_rasterized([axes_data, axes_ref, axes_bias])

        # File paths
        plot_path = self.get_plot_path(plot_type, dataset)
        netcdf_path = get_diagnostic_filename(
            Path(plot_path).stem + "_{pos}", self.cfg
        )
        netcdf_paths = {
            netcdf_path.format(pos="top_left"): cube,
            netcdf_path.format(pos="top_right"): ref_cube,
            netcdf_path.format(pos="bottom"): bias_cube,
        }

        return (plot_path, netcdf_paths)

    def _plot_map_without_ref(self, plot_func, dataset):
        """Plot map plot for single dataset without a reference dataset."""
        plot_type = "map"
        logger.info(
            "Plotting map without reference dataset for '%s'",
            self._get_label(dataset),
        )

        # Make sure that the data has the correct dimensions
        cube = dataset["cube"]
        dim_coords_dat = self._check_cube_dimensions(cube, plot_type)

        # Create plot with desired settings
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg["figure_kwargs"])
            axes = fig.add_subplot(projection=self._get_map_projection())
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs["axes"] = axes
            if plot_func is iris.plot.contourf:
                # see https://github.com/SciTools/cartopy/issues/2457
                # and https://github.com/SciTools/cartopy/issues/2468
                plot_kwargs["transform_first"] = True
                npx = da if cube.has_lazy_data() else np
                cube_to_plot = cube.copy(
                    npx.ma.filled(cube.core_data(), np.nan)
                )
            else:
                cube_to_plot = cube
            plot_map = plot_func(cube_to_plot, **plot_kwargs)
            axes.coastlines()
            gridline_kwargs = self._get_gridline_kwargs(plot_type)
            if gridline_kwargs is not False:
                axes.gridlines(**gridline_kwargs)

            # Print statistics if desired
            self._add_stats(plot_type, axes, dim_coords_dat, dataset)

            # Setup colorbar
            fontsize = (
                self.plots[plot_type]["fontsize"]
                or mpl.rcParams["axes.labelsize"]
            )
            colorbar = fig.colorbar(
                plot_map, ax=axes, **self._get_cbar_kwargs(plot_type)
            )
            colorbar.set_label(
                self._get_cbar_label(plot_type, dataset), fontsize=fontsize
            )
            colorbar.ax.tick_params(labelsize=fontsize)

            # Customize plot
            axes.set_title(self._get_label(dataset))
            fig.suptitle(dataset["long_name"])
            self._process_pyplot_kwargs(plot_type, dataset)

            # Rasterization
            if self.plots[plot_type]["rasterize"]:
                self._set_rasterized([axes])

        # File paths
        plot_path = self.get_plot_path(plot_type, dataset)
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)

        return (plot_path, {netcdf_path: cube})

    def _plot_zonal_mean_profile_with_ref(
        self, plot_func, dataset, ref_dataset
    ):
        """Plot zonal mean profile for single dataset with reference."""
        plot_type = "zonal_mean_profile"
        logger.info(
            "Plotting zonal mean profile with reference dataset '%s' for '%s'",
            self._get_label(ref_dataset),
            self._get_label(dataset),
        )

        # Make sure that the data has the correct dimensions
        cube = dataset["cube"]
        ref_cube = ref_dataset["cube"]
        dim_coords_dat = self._check_cube_dimensions(cube, plot_type)
        dim_coords_ref = self._check_cube_dimensions(ref_cube, plot_type)

        # Create single figure with multiple axes
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg["figure_kwargs"])
            gridspec = GridSpec(
                5, 4, figure=fig, height_ratios=[1.0, 1.0, 0.4, 1.0, 1.0]
            )

            # Options used for all subplots
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            fontsize = (
                self.plots[plot_type]["fontsize"]
                or mpl.rcParams["axes.labelsize"]
            )

            # Plot dataset (top left)
            axes_data = fig.add_subplot(gridspec[0:2, 0:2])
            plot_kwargs["axes"] = axes_data
            plot_data = plot_func(cube, **plot_kwargs)
            axes_data.set_title(self._get_label(dataset), pad=3.0)
            z_coord = cube.coord(axis="Z")
            axes_data.set_ylabel(f"{z_coord.long_name} [{z_coord.units}]")
            if self.plots[plot_type]["log_y"]:
                axes_data.set_yscale("log")
                axes_data.get_yaxis().set_major_formatter(
                    FormatStrFormatter("%.1f")
                )
            if self.plots[plot_type]["show_y_minor_ticklabels"]:
                axes_data.get_yaxis().set_minor_formatter(
                    FormatStrFormatter("%.1f")
                )
            else:
                axes_data.get_yaxis().set_minor_formatter(NullFormatter())
            self._add_stats(plot_type, axes_data, dim_coords_dat, dataset)
            self._process_pyplot_kwargs(plot_type, dataset)

            # Plot reference dataset (top right)
            # Note: make sure to use the same vmin and vmax than the top left
            # plot if a common colorbar is desired
            axes_ref = fig.add_subplot(
                gridspec[0:2, 2:4], sharex=axes_data, sharey=axes_data
            )
            plot_kwargs["axes"] = axes_ref
            if self.plots[plot_type]["common_cbar"]:
                plot_kwargs.setdefault("vmin", plot_data.get_clim()[0])
                plot_kwargs.setdefault("vmax", plot_data.get_clim()[1])
            plot_ref = plot_func(ref_cube, **plot_kwargs)
            axes_ref.set_title(self._get_label(ref_dataset), pad=3.0)
            plt.setp(axes_ref.get_yticklabels(), visible=False)
            self._add_stats(plot_type, axes_ref, dim_coords_ref, ref_dataset)
            self._process_pyplot_kwargs(plot_type, ref_dataset)

            # Add colorbar(s)
            self._add_colorbar(
                plot_type,
                plot_data,
                plot_ref,
                axes_data,
                axes_ref,
                dataset,
                ref_dataset,
            )

            # Plot bias (bottom center)
            bias_cube = cube - ref_cube
            axes_bias = fig.add_subplot(
                gridspec[3:5, 1:3], sharex=axes_data, sharey=axes_data
            )
            plot_kwargs_bias = self._get_plot_kwargs(
                plot_type, dataset, bias=True
            )
            plot_kwargs_bias["axes"] = axes_bias
            plot_bias = plot_func(bias_cube, **plot_kwargs_bias)
            axes_bias.set_title(
                f"{self._get_label(dataset)} - {self._get_label(ref_dataset)}",
                pad=3.0,
            )
            axes_bias.set_xlabel("latitude [°N]")
            axes_bias.set_ylabel(f"{z_coord.long_name} [{z_coord.units}]")
            cbar_kwargs_bias = self._get_cbar_kwargs(plot_type, bias=True)
            cbar_bias = fig.colorbar(
                plot_bias, ax=axes_bias, **cbar_kwargs_bias
            )
            cbar_bias.set_label(
                self._get_cbar_label(plot_type, dataset, bias=True),
                fontsize=fontsize,
            )
            cbar_bias.ax.tick_params(labelsize=fontsize)
            self._add_stats(
                plot_type, axes_bias, dim_coords_dat, dataset, ref_dataset
            )

            # Customize plot
            fig.suptitle(dataset["long_name"])
            self._process_pyplot_kwargs(plot_type, dataset)

            # Rasterization
            if self.plots[plot_type]["rasterize"]:
                self._set_rasterized([axes_data, axes_ref, axes_bias])

        # File paths
        plot_path = self.get_plot_path(plot_type, dataset)
        netcdf_path = get_diagnostic_filename(
            Path(plot_path).stem + "_{pos}", self.cfg
        )
        netcdf_paths = {
            netcdf_path.format(pos="top_left"): cube,
            netcdf_path.format(pos="top_right"): ref_cube,
            netcdf_path.format(pos="bottom"): bias_cube,
        }

        return (plot_path, netcdf_paths)

    def _plot_zonal_mean_profile_without_ref(self, plot_func, dataset):
        """Plot zonal mean profile for single dataset without reference."""
        plot_type = "zonal_mean_profile"
        logger.info(
            "Plotting zonal mean profile without reference dataset for '%s'",
            self._get_label(dataset),
        )

        # Make sure that the data has the correct dimensions
        cube = dataset["cube"]
        dim_coords_dat = self._check_cube_dimensions(cube, plot_type)

        # Create plot with desired settings
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg["figure_kwargs"])
            axes = fig.add_subplot()
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs["axes"] = axes
            plot_zonal_mean_profile = plot_func(cube, **plot_kwargs)

            # Print statistics if desired
            self._add_stats(plot_type, axes, dim_coords_dat, dataset)

            # Setup colorbar
            fontsize = (
                self.plots[plot_type]["fontsize"]
                or mpl.rcParams["axes.labelsize"]
            )
            colorbar = fig.colorbar(
                plot_zonal_mean_profile,
                ax=axes,
                **self._get_cbar_kwargs(plot_type),
            )
            colorbar.set_label(
                self._get_cbar_label(plot_type, dataset), fontsize=fontsize
            )
            colorbar.ax.tick_params(labelsize=fontsize)

            # Customize plot
            axes.set_title(self._get_label(dataset))
            fig.suptitle(dataset["long_name"])
            axes.set_xlabel("latitude [°N]")
            z_coord = cube.coord(axis="Z")
            axes.set_ylabel(f"{z_coord.long_name} [{z_coord.units}]")
            if self.plots[plot_type]["log_y"]:
                axes.set_yscale("log")
                axes.get_yaxis().set_major_formatter(
                    FormatStrFormatter("%.1f")
                )
            if self.plots[plot_type]["show_y_minor_ticklabels"]:
                axes.get_yaxis().set_minor_formatter(
                    FormatStrFormatter("%.1f")
                )
            else:
                axes.get_yaxis().set_minor_formatter(NullFormatter())
            self._process_pyplot_kwargs(plot_type, dataset)

            # Rasterization
            if self.plots[plot_type]["rasterize"]:
                self._set_rasterized([axes])

        # File paths
        plot_path = self.get_plot_path(plot_type, dataset)
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)

        return (plot_path, {netcdf_path: cube})

    def _plot_hovmoeller_z_vs_time_without_ref(self, plot_func, dataset):
        """Plot Hovmoeller Z vs. time for single dataset without reference."""
        plot_type = "hovmoeller_z_vs_time"
        logger.info(
            "Plotting Hovmoeller Z vs. time without reference dataset"
            " for '%s'",
            self._get_label(dataset),
        )

        # Make sure that the data has the correct dimensions
        cube = dataset["cube"]
        dim_coords_dat = self._check_cube_dimensions(cube, plot_type)

        # Create plot with desired settings
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg["figure_kwargs"])
            axes = fig.add_subplot()
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs["axes"] = axes
            plot_hovmoeller = plot_func(cube, **plot_kwargs)

            # Print statistics if desired
            self._add_stats(plot_type, axes, dim_coords_dat, dataset)

            # Setup colorbar
            fontsize = (
                self.plots[plot_type]["fontsize"]
                or mpl.rcParams["axes.labelsize"]
            )
            colorbar = fig.colorbar(
                plot_hovmoeller, ax=axes, **self._get_cbar_kwargs(plot_type)
            )
            colorbar.set_label(
                self._get_cbar_label(plot_type, dataset), fontsize=fontsize
            )
            colorbar.ax.tick_params(labelsize=fontsize)

            # Customize plot
            axes.set_title(self._get_label(dataset))
            fig.suptitle(dataset["long_name"])
            z_coord = cube.coord(axis="Z")
            axes.set_ylabel(f"{z_coord.long_name} [{z_coord.units}]")
            if self.plots[plot_type]["log_y"]:
                axes.set_yscale("log")
                axes.get_yaxis().set_major_formatter(
                    FormatStrFormatter("%.1f")
                )
            if self.plots[plot_type]["show_y_minor_ticklabels"]:
                axes.get_yaxis().set_minor_formatter(
                    FormatStrFormatter("%.1f")
                )
            else:
                axes.get_yaxis().set_minor_formatter(NullFormatter())
            if self.plots[plot_type]["time_format"] is not None:
                axes.get_xaxis().set_major_formatter(
                    mdates.DateFormatter(self.plots[plot_type]["time_format"])
                )
            axes.set_xlabel("time")
            self._process_pyplot_kwargs(plot_type, dataset)

            # Rasterization
            if self.plots[plot_type]["rasterize"]:
                self._set_rasterized([axes])

        # File paths
        plot_path = self.get_plot_path(plot_type, dataset)
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)

        return (plot_path, {netcdf_path: cube})

    def _plot_hovmoeller_z_vs_time_with_ref(
        self, plot_func, dataset, ref_dataset
    ):
        """Plot Hovmoeller Z vs. time for single dataset with reference."""
        plot_type = "hovmoeller_z_vs_time"
        logger.info(
            "Plotting Hovmoeller z vs. time with reference dataset"
            " '%s' for '%s'",
            self._get_label(ref_dataset),
            self._get_label(dataset),
        )

        # Make sure that the data has the correct dimensions
        cube = dataset["cube"]
        ref_cube = ref_dataset["cube"]
        dim_coords_dat = self._check_cube_dimensions(cube, plot_type)
        dim_coords_ref = self._check_cube_dimensions(ref_cube, plot_type)

        # Create single figure with multiple axes
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg["figure_kwargs"])
            gridspec = GridSpec(
                5, 4, figure=fig, height_ratios=[1.0, 1.0, 0.4, 1.0, 1.0]
            )

            # Options used for all subplots
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            fontsize = (
                self.plots[plot_type]["fontsize"]
                or mpl.rcParams["axes.labelsize"]
            )

            # Plot dataset (top left)
            axes_data = fig.add_subplot(gridspec[0:2, 0:2])
            plot_kwargs["axes"] = axes_data
            plot_data = plot_func(cube, **plot_kwargs)
            axes_data.set_title(self._get_label(dataset), pad=3.0)
            z_coord = cube.coord(axis="Z")
            axes_data.set_ylabel(f"{z_coord.long_name} [{z_coord.units}]")
            if self.plots[plot_type]["log_y"]:
                axes_data.set_yscale("log")
                axes_data.get_yaxis().set_major_formatter(
                    FormatStrFormatter("%.1f")
                )
            if self.plots[plot_type]["show_y_minor_ticklabels"]:
                axes_data.get_yaxis().set_minor_formatter(
                    FormatStrFormatter("%.1f")
                )
            else:
                axes_data.get_yaxis().set_minor_formatter(NullFormatter())
            if self.plots[plot_type]["time_format"] is not None:
                axes_data.get_xaxis().set_major_formatter(
                    mdates.DateFormatter(self.plots[plot_type]["time_format"])
                )
            self._add_stats(plot_type, axes_data, dim_coords_dat, dataset)
            self._process_pyplot_kwargs(plot_type, dataset)

            # Plot reference dataset (top right)
            # Note: make sure to use the same vmin and vmax than the top left
            # plot if a common colorbar is desired
            axes_ref = fig.add_subplot(
                gridspec[0:2, 2:4], sharex=axes_data, sharey=axes_data
            )
            plot_kwargs["axes"] = axes_ref
            if self.plots[plot_type]["common_cbar"]:
                plot_kwargs.setdefault("vmin", plot_data.get_clim()[0])
                plot_kwargs.setdefault("vmax", plot_data.get_clim()[1])
            plot_ref = plot_func(ref_cube, **plot_kwargs)
            axes_ref.set_title(self._get_label(ref_dataset), pad=3.0)
            plt.setp(axes_ref.get_yticklabels(), visible=False)
            self._add_stats(plot_type, axes_ref, dim_coords_ref, ref_dataset)
            self._process_pyplot_kwargs(plot_type, ref_dataset)

            # Add colorbar(s)
            self._add_colorbar(
                plot_type,
                plot_data,
                plot_ref,
                axes_data,
                axes_ref,
                dataset,
                ref_dataset,
            )

            # Plot bias (bottom center)
            bias_cube = cube - ref_cube
            axes_bias = fig.add_subplot(
                gridspec[3:5, 1:3], sharex=axes_data, sharey=axes_data
            )
            plot_kwargs_bias = self._get_plot_kwargs(
                plot_type, dataset, bias=True
            )
            plot_kwargs_bias["axes"] = axes_bias
            plot_bias = plot_func(bias_cube, **plot_kwargs_bias)
            axes_bias.set_title(
                f"{self._get_label(dataset)} - {self._get_label(ref_dataset)}",
                pad=3.0,
            )
            axes_bias.set_xlabel("time")
            axes_bias.set_ylabel(f"{z_coord.long_name} [{z_coord.units}]")
            cbar_kwargs_bias = self._get_cbar_kwargs(plot_type, bias=True)
            cbar_bias = fig.colorbar(
                plot_bias, ax=axes_bias, **cbar_kwargs_bias
            )
            cbar_bias.set_label(
                self._get_cbar_label(plot_type, dataset, bias=True),
                fontsize=fontsize,
            )
            cbar_bias.ax.tick_params(labelsize=fontsize)
            self._add_stats(
                plot_type, axes_bias, dim_coords_dat, dataset, ref_dataset
            )

            # Customize plot
            fig.suptitle(dataset["long_name"])
            self._process_pyplot_kwargs(plot_type, dataset)

            # Rasterization
            if self.plots[plot_type]["rasterize"]:
                self._set_rasterized([axes_data, axes_ref, axes_bias])

        # File paths
        plot_path = self.get_plot_path(plot_type, dataset)
        netcdf_path = get_diagnostic_filename(
            Path(plot_path).stem + "_{pos}", self.cfg
        )
        netcdf_paths = {
            netcdf_path.format(pos="top_left"): cube,
            netcdf_path.format(pos="top_right"): ref_cube,
            netcdf_path.format(pos="bottom"): bias_cube,
        }

        return (plot_path, netcdf_paths)

    def _plot_hovmoeller_time_vs_lat_or_lon_with_ref(
        self, plot_func, dataset, ref_dataset
    ):
        """Plot the hovmoeller profile for single dataset with reference."""
        plot_type = "hovmoeller_time_vs_lat_or_lon"
        logger.info(
            "Plotting Hovmoeller plots with reference dataset '%s' for '%s'",
            self._get_label(ref_dataset),
            self._get_label(dataset),
        )

        # Make sure that the data has the correct dimensions
        cube = dataset["cube"]
        ref_cube = ref_dataset["cube"]
        dim_coords_dat = self._check_cube_dimensions(cube, plot_type)
        self._check_cube_dimensions(ref_cube, plot_type)
        if "latitude" in dim_coords_dat:
            non_time_label = "latitude [°N]"
        else:
            non_time_label = "longitude [°E]"

        # Create single figure with multiple axes
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg["figure_kwargs"])
            gridspec = GridSpec(
                5, 4, figure=fig, height_ratios=[1.0, 1.0, 0.4, 1.0, 1.0]
            )

            # Options used for all subplots
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            fontsize = (
                self.plots[plot_type]["fontsize"]
                or mpl.rcParams["axes.labelsize"]
            )

            # Plot dataset (top left)
            axes_data = fig.add_subplot(gridspec[0:2, 0:2])
            plot_kwargs["axes"] = axes_data
            if self.plots[plot_type]["time_on"] == "x-axis":
                plot_kwargs["coords"] = list(dim_coords_dat)
                x_label = "time"
                y_label = non_time_label
                time_axis = axes_data.get_xaxis()
            else:
                plot_kwargs["coords"] = list(reversed(dim_coords_dat))
                x_label = non_time_label
                y_label = "time"
                time_axis = axes_data.get_yaxis()
            plot_data = plot_func(cube, **plot_kwargs)
            axes_data.set_title(self._get_label(dataset), pad=3.0)
            axes_data.set_ylabel(y_label)
            if self.plots[plot_type]["time_format"] is not None:
                time_axis.set_major_formatter(
                    mdates.DateFormatter(self.plots[plot_type]["time_format"])
                )
            if self.plots[plot_type]["show_y_minor_ticks"]:
                axes_data.get_yaxis().set_minor_locator(AutoMinorLocator())
            if self.plots[plot_type]["show_x_minor_ticks"]:
                axes_data.get_xaxis().set_minor_locator(AutoMinorLocator())
            self._process_pyplot_kwargs(plot_type, dataset)

            # Plot reference dataset (top right)
            # Note: make sure to use the same vmin and vmax than the top left
            # plot if a common colorbar is desired
            axes_ref = fig.add_subplot(
                gridspec[0:2, 2:4], sharex=axes_data, sharey=axes_data
            )
            plot_kwargs["axes"] = axes_ref
            if self.plots[plot_type]["common_cbar"]:
                plot_kwargs.setdefault("vmin", plot_data.get_clim()[0])
                plot_kwargs.setdefault("vmax", plot_data.get_clim()[1])
            plot_ref = plot_func(ref_cube, **plot_kwargs)
            axes_ref.set_title(self._get_label(ref_dataset), pad=3.0)
            plt.setp(axes_ref.get_yticklabels(), visible=False)
            self._process_pyplot_kwargs(plot_type, ref_dataset)

            # Add colorbar(s)
            self._add_colorbar(
                plot_type,
                plot_data,
                plot_ref,
                axes_data,
                axes_ref,
                dataset,
                ref_dataset,
            )

            # Plot bias (bottom center)
            bias_cube = cube - ref_cube
            axes_bias = fig.add_subplot(
                gridspec[3:5, 1:3], sharex=axes_data, sharey=axes_data
            )
            plot_kwargs_bias = self._get_plot_kwargs(
                plot_type, dataset, bias=True
            )
            plot_kwargs_bias["axes"] = axes_bias
            plot_kwargs_bias["coords"] = plot_kwargs["coords"]
            plot_bias = plot_func(bias_cube, **plot_kwargs_bias)
            axes_bias.set_title(
                f"{self._get_label(dataset)} - {self._get_label(ref_dataset)}",
                pad=3.0,
            )
            axes_bias.set_xlabel(x_label)
            axes_bias.set_ylabel(y_label)
            cbar_kwargs_bias = self._get_cbar_kwargs(plot_type, bias=True)
            cbar_bias = fig.colorbar(
                plot_bias, ax=axes_bias, **cbar_kwargs_bias
            )
            cbar_bias.set_label(
                self._get_cbar_label(plot_type, dataset, bias=True),
                fontsize=fontsize,
            )
            cbar_bias.ax.tick_params(labelsize=fontsize)

            # Customize plot
            fig.suptitle(dataset["long_name"])
            self._process_pyplot_kwargs(plot_type, dataset)

            # Rasterization
            if self.plots[plot_type]["rasterize"]:
                self._set_rasterized([axes_data, axes_ref, axes_bias])

        # File paths
        plot_path = self.get_plot_path(plot_type, dataset)
        netcdf_path = get_diagnostic_filename(
            Path(plot_path).stem + "_{pos}", self.cfg
        )
        netcdf_paths = {
            netcdf_path.format(pos="top_left"): cube,
            netcdf_path.format(pos="top_right"): ref_cube,
            netcdf_path.format(pos="bottom"): bias_cube,
        }

        return (plot_path, netcdf_paths)

    def _plot_hovmoeller_time_vs_lat_or_lon_without_ref(
        self, plot_func, dataset
    ):
        """Plot time vs zonal or meridional Hovmoeller without reference."""
        plot_type = "hovmoeller_time_vs_lat_or_lon"
        logger.info(
            "Plotting Hovmoeller plots without reference dataset for '%s'",
            self._get_label(dataset),
        )

        # Make sure that the data has the correct dimensions
        cube = dataset["cube"]
        dim_coords_dat = self._check_cube_dimensions(cube, plot_type)
        if "latitude" in dim_coords_dat:
            non_time_label = "latitude [°N]"
        else:
            non_time_label = "longitude [°E]"

        # Create plot with desired settings
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg["figure_kwargs"])
            axes = fig.add_subplot()
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs["axes"] = axes

            # Put time on desired axis
            if self.plots[plot_type]["time_on"] == "x-axis":
                plot_kwargs["coords"] = list(dim_coords_dat)
                x_label = "time"
                y_label = non_time_label
                time_axis = axes.get_xaxis()
            else:
                plot_kwargs["coords"] = list(reversed(dim_coords_dat))
                x_label = non_time_label
                y_label = "time"
                time_axis = axes.get_yaxis()
            plot_hovmoeller = plot_func(cube, **plot_kwargs)

            # Setup colorbar
            fontsize = (
                self.plots[plot_type]["fontsize"]
                or mpl.rcParams["axes.labelsize"]
            )
            colorbar = fig.colorbar(
                plot_hovmoeller, ax=axes, **self._get_cbar_kwargs(plot_type)
            )
            colorbar.set_label(
                self._get_cbar_label(plot_type, dataset), fontsize=fontsize
            )
            colorbar.ax.tick_params(labelsize=fontsize)

            # Customize plot
            axes.set_title(self._get_label(dataset))
            fig.suptitle(dataset["long_name"])
            axes.set_xlabel(x_label)
            axes.set_ylabel(y_label)
            if self.plots[plot_type]["time_format"] is not None:
                time_axis.set_major_formatter(
                    mdates.DateFormatter(self.plots[plot_type]["time_format"])
                )
            if self.plots[plot_type]["show_y_minor_ticks"]:
                axes.get_yaxis().set_minor_locator(AutoMinorLocator())
            if self.plots[plot_type]["show_x_minor_ticks"]:
                axes.get_xaxis().set_minor_locator(AutoMinorLocator())
            self._process_pyplot_kwargs(plot_type, dataset)

            # Rasterization
            if self.plots[plot_type]["rasterize"]:
                self._set_rasterized([axes])

        # File paths
        plot_path = self.get_plot_path(plot_type, dataset)
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)
        return (plot_path, {netcdf_path: cube})

    def _plot_hovmoeller_anncyc_vs_lat_or_lon_with_ref(
        self, plot_func, dataset, ref_dataset
    ):
        """Plot the hovmoeller profile for single dataset with reference."""
        plot_type = "hovmoeller_anncyc_vs_lat_or_lon"
        logger.info(
            "Plotting Hovmoeller plots with reference dataset '%s' for '%s'",
            self._get_label(ref_dataset),
            self._get_label(dataset),
        )

        # Make sure that the data has the correct dimensions
        cube = dataset["cube"]
        ref_cube = ref_dataset["cube"]
        dim_coords_dat = self._check_cube_dimensions(cube, plot_type)
        self._check_cube_dimensions(ref_cube, plot_type)
        if "latitude" in dim_coords_dat:
            non_time_label = "latitude [°N]"
        else:
            non_time_label = "longitude [°E]"

        # Create single figure with multiple axes
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg["figure_kwargs"])
            gridspec = GridSpec(
                5, 4, figure=fig, height_ratios=[1.0, 1.0, 0.4, 1.0, 1.0]
            )

            # Options used for all subplots
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            fontsize = (
                self.plots[plot_type]["fontsize"]
                or mpl.rcParams["axes.labelsize"]
            )

            # Plot dataset (top left)
            axes_data = fig.add_subplot(gridspec[0:2, 0:2])
            plot_kwargs["axes"] = axes_data
            if self.plots[plot_type]["time_on"] == "x-axis":
                plot_kwargs["coords"] = list(dim_coords_dat)
                x_label = "month"
                y_label = non_time_label
            else:
                plot_kwargs["coords"] = list(reversed(dim_coords_dat))
                x_label = non_time_label
                y_label = "month"
            plot_data = plot_func(cube, **plot_kwargs)
            axes_data.set_title(self._get_label(dataset), pad=3.0)
            axes_data.set_ylabel(y_label)
            if self.plots[plot_type]["show_y_minor_ticks"]:
                axes_data.get_yaxis().set_minor_locator(AutoMinorLocator())
            if self.plots[plot_type]["show_x_minor_ticks"]:
                axes_data.get_xaxis().set_minor_locator(AutoMinorLocator())
            self._process_pyplot_kwargs(plot_type, dataset)

            # Plot reference dataset (top right)
            # Note: make sure to use the same vmin and vmax than the top left
            # plot if a common colorbar is desired
            axes_ref = fig.add_subplot(
                gridspec[0:2, 2:4], sharex=axes_data, sharey=axes_data
            )
            plot_kwargs["axes"] = axes_ref
            if self.plots[plot_type]["common_cbar"]:
                plot_kwargs.setdefault("vmin", plot_data.get_clim()[0])
                plot_kwargs.setdefault("vmax", plot_data.get_clim()[1])
            plot_ref = plot_func(ref_cube, **plot_kwargs)
            axes_ref.set_title(self._get_label(ref_dataset), pad=3.0)
            plt.setp(axes_ref.get_yticklabels(), visible=False)
            self._process_pyplot_kwargs(plot_type, ref_dataset)

            # Add colorbar(s)
            self._add_colorbar(
                plot_type,
                plot_data,
                plot_ref,
                axes_data,
                axes_ref,
                dataset,
                ref_dataset,
            )

            # Plot bias (bottom center)
            bias_cube = cube - ref_cube
            axes_bias = fig.add_subplot(
                gridspec[3:5, 1:3], sharex=axes_data, sharey=axes_data
            )
            plot_kwargs_bias = self._get_plot_kwargs(
                plot_type, dataset, bias=True
            )
            plot_kwargs_bias["axes"] = axes_bias
            plot_kwargs_bias["coords"] = plot_kwargs["coords"]
            plot_bias = plot_func(bias_cube, **plot_kwargs_bias)
            axes_bias.set_title(
                f"{self._get_label(dataset)} - {self._get_label(ref_dataset)}",
                pad=3.0,
            )
            axes_bias.set_xlabel(x_label)
            axes_bias.set_ylabel(y_label)
            cbar_kwargs_bias = self._get_cbar_kwargs(plot_type, bias=True)
            cbar_bias = fig.colorbar(
                plot_bias, ax=axes_bias, **cbar_kwargs_bias
            )
            cbar_bias.set_label(
                self._get_cbar_label(plot_type, dataset, bias=True),
                fontsize=fontsize,
            )
            cbar_bias.ax.tick_params(labelsize=fontsize)

            # Customize plot
            fig.suptitle(dataset["long_name"])
            self._process_pyplot_kwargs(plot_type, dataset)

            # Rasterization
            if self.plots[plot_type]["rasterize"]:
                self._set_rasterized([axes_data, axes_ref, axes_bias])

        # File paths
        plot_path = self.get_plot_path(plot_type, dataset)
        netcdf_path = get_diagnostic_filename(
            Path(plot_path).stem + "_{pos}", self.cfg
        )
        netcdf_paths = {
            netcdf_path.format(pos="top_left"): cube,
            netcdf_path.format(pos="top_right"): ref_cube,
            netcdf_path.format(pos="bottom"): bias_cube,
        }

        return (plot_path, netcdf_paths)

    def _plot_hovmoeller_anncyc_vs_lat_or_lon_without_ref(
        self, plot_func, dataset
    ):
        """Plot annual cylce vs zonal or meridional Hovmoeller without reference."""
        plot_type = "hovmoeller_anncyc_vs_lat_or_lon"
        logger.info(
            "Plotting Hovmoeller plots without reference dataset for '%s'",
            self._get_label(dataset),
        )

        # Make sure that the data has the correct dimensions
        cube = dataset["cube"]
        dim_coords_dat = self._check_cube_dimensions(cube, plot_type)
        if "latitude" in dim_coords_dat:
            non_time_label = "latitude [°N]"
        else:
            non_time_label = "longitude [°E]"

        # Create plot with desired settings
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg["figure_kwargs"])
            axes = fig.add_subplot()
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs["axes"] = axes

            # Put time on desired axis
            if self.plots[plot_type]["time_on"] == "x-axis":
                plot_kwargs["coords"] = list(dim_coords_dat)
                x_label = "month"
                y_label = non_time_label
            else:
                plot_kwargs["coords"] = list(reversed(dim_coords_dat))
                x_label = non_time_label
                y_label = "month"
            plot_hovmoeller = plot_func(cube, **plot_kwargs)

            # Setup colorbar
            fontsize = (
                self.plots[plot_type]["fontsize"]
                or mpl.rcParams["axes.labelsize"]
            )
            colorbar = fig.colorbar(
                plot_hovmoeller, ax=axes, **self._get_cbar_kwargs(plot_type)
            )
            colorbar.set_label(
                self._get_cbar_label(plot_type, dataset), fontsize=fontsize
            )
            colorbar.ax.tick_params(labelsize=fontsize)

            # Customize plot
            axes.set_title(self._get_label(dataset))
            fig.suptitle(dataset["long_name"])
            axes.set_xlabel(x_label)
            axes.set_ylabel(y_label)
            if self.plots[plot_type]["show_y_minor_ticks"]:
                axes.get_yaxis().set_minor_locator(AutoMinorLocator())
            if self.plots[plot_type]["show_x_minor_ticks"]:
                axes.get_xaxis().set_minor_locator(AutoMinorLocator())
            self._process_pyplot_kwargs(plot_type, dataset)

            # Rasterization
            if self.plots[plot_type]["rasterize"]:
                self._set_rasterized([axes])

        # File paths
        plot_path = self.get_plot_path(plot_type, dataset)
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)
        return (plot_path, {netcdf_path: cube})

    def _plot_benchmarking_map(
        self, plot_func, dataset, percentile_dataset, metric
    ):
        """Plot benchmarking map plot."""
        plot_type = "benchmarking_map"
        logger.info(
            "Plotting benchmarking map for '%s'", self._get_label(dataset)
        )

        # Make sure that the data has the correct dimensions
        cube = dataset["cube"]

        # Create plot with desired settings
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg["figure_kwargs"])
            axes = fig.add_subplot(
                projection=self._get_benchmarking_projection()
            )
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs["axes"] = axes
            plot_kwargs["extend"] = "both"

            # apply stippling (dots) to all grid cells that do not exceed
            # the upper percentile given by 'percentile_dataset[]'
            mask_cube = self._get_benchmark_mask(
                cube, percentile_dataset, metric
            )
            hatching_plot_kwargs = {
                "colors": "none",
                "levels": [0.5, 1.5],
                "hatches": ["......"],
            }

            if plot_func is iris.plot.contourf:
                # see https://github.com/SciTools/cartopy/issues/2457
                # and https://github.com/SciTools/cartopy/issues/2468
                plot_kwargs["transform_first"] = True
                hatching_plot_kwargs["transform_first"] = True
                npx = da if cube.has_lazy_data() else np
                cube_to_plot = cube.copy(
                    npx.ma.filled(cube.core_data(), np.nan)
                )
                mask_cube_to_plot = mask_cube.copy(
                    npx.ma.filled(mask_cube.core_data(), np.nan)
                )
            else:
                cube_to_plot = cube
                mask_cube_to_plot = mask_cube

            # Plot
            plot_map = plot_func(cube_to_plot, **plot_kwargs)
            hatching = plot_func(mask_cube_to_plot, **hatching_plot_kwargs)

            # set color for stippling to 'black' (default = 'white')
            hatching.set_edgecolor("black")
            hatching.set_linewidth(0.0)

            axes.coastlines()

            # Setup colorbar
            fontsize = (
                self.plots[plot_type]["fontsize"]
                or mpl.rcParams["axes.labelsize"]
            )
            colorbar = fig.colorbar(
                plot_map, ax=axes, **self._get_cbar_kwargs(plot_type)
            )
            colorbar.set_label(
                self._get_cbar_label(plot_type, dataset), fontsize=fontsize
            )
            colorbar.ax.tick_params(labelsize=fontsize)

            # Customize plot
            axes.set_title(self._get_label(dataset))
            fig.suptitle(dataset["long_name"])
            self._process_pyplot_kwargs(plot_type, dataset)

            # Rasterization
            if self.plots[plot_type]["rasterize"]:
                self._set_rasterized([axes])

        # File paths
        plot_path = self.get_plot_path(plot_type, dataset)
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)

        return (plot_path, {netcdf_path: cube})

    def _plot_benchmarking_boxplot(self, df, cubes, variables, datasets):
        """Plot benchmarking boxplot."""
        plot_type = "benchmarking_boxplot"
        logger.info(
            "Plotting benchmarking boxplot for '%s'",
            self._get_label(datasets[0]),
        )

        # Create plot with desired settings
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg["figure_kwargs"])
            metric = cubes[0].long_name.partition("of")[0]
            fig.suptitle(f"{metric}of {self._get_label(datasets[0])}")

            sns.set_style("darkgrid")

            for i, var in enumerate(variables):
                axes = plt.subplot(1, len(variables), i + 1)
                plot_kwargs = self._get_plot_kwargs(plot_type, datasets[i])
                plot_kwargs["axes"] = axes

                plot_boxplot = sns.boxplot(data=df[df["Variable"] == var])
                plot_boxplot.set(xticklabels=[])

                plt.scatter(
                    0,
                    cubes[i].data,
                    marker="x",
                    s=200,
                    linewidths=2,
                    color="red",
                    zorder=3,
                )

                plt.xlabel(var)
                if cubes[i].units != 1:
                    plt.ylabel(cubes[i].units)

                # Customize plot
                self._process_pyplot_kwargs(plot_type, datasets[i])

        # File paths
        datasets[0]["variable_group"] = datasets[0]["short_name"].partition(
            "_"
        )[0]
        plot_path = self.get_plot_path(plot_type, datasets[0])
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)

        return (plot_path, {netcdf_path: cubes[0]})

    def _plot_benchmarking_zonal(
        self, plot_func, dataset, percentile_dataset, metric
    ):
        """Plot benchmarking zonal mean profile."""
        plot_type = "benchmarking_zonal"
        logger.info(
            "Plotting benchmarking zonal mean profile for '%s'",
            self._get_label(dataset),
        )

        # Make sure that the data has the correct dimensions
        cube = dataset["cube"]

        # Create plot with desired settings
        with mpl.rc_context(self._get_custom_mpl_rc_params(plot_type)):
            fig = plt.figure(**self.cfg["figure_kwargs"])
            axes = fig.add_subplot()
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs["axes"] = axes
            plot_kwargs["extend"] = "both"
            plot_benchmarking_zonal = plot_func(cube, **plot_kwargs)

            # apply stippling (dots) to all grid cells that do not exceed
            # the upper percentile given by 'percentile_dataset[]'

            mask_cube = self._get_benchmark_mask(
                cube, percentile_dataset, metric
            )
            hatching = plot_func(
                mask_cube,
                colors="none",
                levels=[0.5, 1.5],
                hatches=["......"],
            )

            # set color for stippling to 'black' (default = 'white')
            hatching.set_edgecolor("black")
            hatching.set_linewidth(0.0)

            # Setup colorbar
            fontsize = (
                self.plots[plot_type]["fontsize"]
                or mpl.rcParams["axes.labelsize"]
            )
            colorbar = fig.colorbar(
                plot_benchmarking_zonal,
                ax=axes,
                **self._get_cbar_kwargs(plot_type),
            )
            colorbar.set_label(
                self._get_cbar_label(plot_type, dataset), fontsize=fontsize
            )
            colorbar.ax.tick_params(labelsize=fontsize)

            # Customize plot
            axes.set_title(self._get_label(dataset))
            fig.suptitle(dataset["long_name"])
            axes.set_xlabel("latitude [°N]")
            z_coord = cube.coord(axis="Z")
            axes.set_ylabel(f"{z_coord.long_name} [{z_coord.units}]")
            if self.plots[plot_type]["log_y"]:
                axes.set_yscale("log")
                axes.get_yaxis().set_major_formatter(
                    FormatStrFormatter("%.1f")
                )
            if self.plots[plot_type]["show_y_minor_ticklabels"]:
                axes.get_yaxis().set_minor_formatter(
                    FormatStrFormatter("%.1f")
                )
            else:
                axes.get_yaxis().set_minor_formatter(NullFormatter())
            self._process_pyplot_kwargs(plot_type, dataset)

            # Rasterization
            if self.plots[plot_type]["rasterize"]:
                self._set_rasterized([axes])

        # File paths
        plot_path = self.get_plot_path(plot_type, dataset)
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)

        return (plot_path, {netcdf_path: cube})

    def _process_pyplot_kwargs(self, plot_type, dataset):
        """Process functions for :mod:`matplotlib.pyplot`."""
        pyplot_kwargs = self.plots[plot_type]["pyplot_kwargs"]
        for func, arg in pyplot_kwargs.items():
            if isinstance(arg, str):
                arg = self._fill_facet_placeholders(
                    arg,
                    dataset,
                    f"pyplot_kwargs of {plot_type} '{func}: {arg}'",
                )
            if arg is None:
                getattr(plt, func)()
            elif isinstance(arg, dict):
                getattr(plt, func)(**arg)
            else:
                getattr(plt, func)(arg)

    @staticmethod
    def _check_cube_dimensions(cube, plot_type):
        """Check that cube has correct dimensional variables."""
        expected_dimensions_dict = {
            "annual_cycle": (["month_number"],),
            "benchmarking_boxplot": ([""],),
            "diurnal_cycle": (["hour"],),
            "map": (["latitude", "longitude"],),
            "benchmarking_annual_cycle": (["month_number"],),
            "benchmarking_diurnal_cycle": (["hour"],),
            "benchmarking_map": (["latitude", "longitude"],),
            "benchmarking_timeseries": (["time"],),
            "benchmarking_zonal": (
                ["latitude", "air_pressure"],
                ["latitude", "altitude"],
            ),
            "zonal_mean_profile": (
                ["latitude", "air_pressure"],
                ["latitude", "altitude"],
            ),
            "timeseries": (["time"],),
            "1d_profile": (["air_pressure"], ["altitude"]),
            "variable_vs_lat": (["latitude"],),
            "hovmoeller_z_vs_time": (
                ["time", "air_pressure"],
                ["time", "altitude"],
            ),
            "hovmoeller_time_vs_lat_or_lon": (
                ["time", "latitude"],
                ["time", "longitude"],
            ),
            "hovmoeller_anncyc_vs_lat_or_lon": (
                ["month_number", "latitude"],
                ["month_number", "longitude"],
            ),
        }
        if plot_type not in expected_dimensions_dict:
            raise NotImplementedError(f"plot_type '{plot_type}' not supported")
        expected_dimensions = expected_dimensions_dict[plot_type]
        for dims in expected_dimensions:
            cube_dims = [cube.coords(dim, dim_coords=True) for dim in dims]
            if all(cube_dims) and cube.ndim == len(dims):
                return dims
        expected_dims_str = " or ".join(
            [str(dims) for dims in expected_dimensions]
        )
        raise ValueError(
            f"Expected cube that exactly has the dimensional coordinates "
            f"{expected_dims_str}, got {cube.summary(shorten=True)}"
        )

    @staticmethod
    def _fill_facet_placeholders(string, dataset, description):
        """Fill facet placeholders."""
        try:
            string = string.format(**dataset)
        except KeyError as exc:
            raise ValueError(
                f"Not all necessary facets in {description} available for "
                f"dataset\n{pformat(dataset)}"
            ) from exc
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
                multi_dataset_facets[key] = f"ambiguous_{key}"
        return multi_dataset_facets

    def _get_reference_dataset(self, datasets):
        """Extract reference dataset."""
        variable = datasets[0][self.cfg["group_variables_by"]]
        ref_datasets = [
            d for d in datasets if d.get("reference_for_monitor_diags", False)
        ]
        if len(ref_datasets) > 1:
            raise ValueError(
                f"Expected at most 1 reference dataset (with "
                f"'reference_for_monitor_diags: true' for variable "
                f"'{variable}', got {len(ref_datasets):d}"
            )
        if ref_datasets:
            return ref_datasets[0]
        return None

    def _get_benchmarking_reference(self, datasets):
        """Extract reference dataset for calculation of benchmarking metric."""
        variable = datasets[0][self.cfg["group_variables_by"]]
        ref_datasets = [
            d for d in datasets if d.get("reference_for_metric", False)
        ]

        if len(ref_datasets) == 1:
            return ref_datasets[0]

        # try variable attribute "reference_dataset"
        for dataset in datasets:
            print(dataset.get("reference_dataset"))
            print(dataset.get("dataset"))
            if dataset.get("reference_dataset") == dataset.get("dataset"):
                ref_datasets = dataset
                break
        if len(ref_datasets) != 1:
            raise ValueError(
                f"Expected exactly 1 reference dataset for variable "
                f"'{variable}', got {len(ref_datasets)}"
            )
        return None

    def _get_benchmark_datasets(self, datasets):
        """Get dataset to be benchmarked."""
        variable = datasets[0][self.cfg["group_variables_by"]]
        benchmark_datasets = [
            d for d in datasets if d.get("benchmark_dataset", False)
        ]
        if len(benchmark_datasets) >= 1:
            return benchmark_datasets

        raise ValueError(
            f"Expected at least 1 benchmark dataset (with "
            f"'benchmark_dataset: true' for variable "
            f"'{variable}'), got {len(benchmark_datasets):d}"
        )

    def _get_benchmark_group(self, datasets):
        """Get datasets for benchmarking."""
        benchmark_datasets = [
            d
            for d in datasets
            if not (
                d.get("benchmark_dataset", False)
                or d.get("reference_for_metric", False)
            )
        ]
        return benchmark_datasets

    def _get_benchmark_mask(self, cube, percentile_dataset, metric):
        """Create mask for benchmarking cube depending on metric."""
        mask_cube = cube.copy()

        idx0 = 0  # index largest percentile
        idx1 = len(percentile_dataset) - 1  # index smallest percentile

        if metric == "bias":
            maxabs_perc = np.maximum(
                np.abs(percentile_dataset[idx0].data),
                np.abs(percentile_dataset[idx1].data),
            )
            mask = np.where(np.abs(cube.data) >= maxabs_perc, 0, 1)
        elif metric == "emd":
            mask = np.where(cube.data >= percentile_dataset[idx0].data, 0, 1)
        elif metric == "pearsonr":
            mask = np.where(cube.data <= percentile_dataset[idx0].data, 0, 1)
        elif metric == "rmse":
            mask = np.where(cube.data >= percentile_dataset[idx0].data, 0, 1)
        else:
            raise ValueError(
                f"Could not create benchmarking mask, unknown benchmarking "
                f"metric: '{metric}'"
            )

        mask_cube.data = mask
        return mask_cube

    def _get_benchmark_metric(self, datasets):
        """Get benchmarking metric."""
        short_name = datasets[0].get("short_name", "")
        if "rmse" in short_name:
            metric = "rmse"
        elif "pearsonr" in short_name:
            metric = "pearsonr"
        elif "emd" in short_name:
            metric = "emd"
        else:
            metric = "bias"  # default
            logger.info(
                "Could not determine metric from short_name, "
                "assuming benchmarking metric = %s",
                metric,
            )
        return metric

    def _get_benchmark_percentiles(self, datasets):
        """Get percentile datasets from multi-model statistics preprocessor."""
        variable = datasets[0][self.cfg["group_variables_by"]]
        percentiles = []
        for dataset in datasets:
            statistics = dataset.get("multi_model_statistics")
            if statistics:
                if "Percentile" in statistics:
                    percentiles.append(dataset)

        # *** sort percentiles by size ***

        # get percentiles as integers
        iperc = []
        for dataset in percentiles:
            stat = dataset.get("multi_model_statistics")
            perc = stat.replace("MultiModelPercentile", "")
            iperc.append(int(perc))

        idx = list(range(len(percentiles)))
        # sort list of percentile datasets by percentile with highest
        # percentile first (descending order)
        zipped_pairs = zip(iperc, idx)
        zval = [x for _, x in sorted(zipped_pairs, reverse=True)]
        perc_sorted = [percentiles[i] for i in zval]
        percentiles = perc_sorted

        # get number of percentiles expected depending on benchmarking metric

        metric = self._get_benchmark_metric(datasets)

        if metric == "bias":
            numperc = 2
        elif metric == "rmse":
            numperc = 1
        elif metric == "pearsonr":
            numperc = 1
        elif metric == "emd":
            numperc = 1
        else:
            raise ValueError(f"Unknown benchmarking metric: '{metric}'.")

        if len(percentiles) >= numperc:
            return percentiles

        raise ValueError(
            f"Expected at least '{numperc}' percentile datasets (created "
            f"'with multi-model statistics preprocessor for variable "
            f"'{variable}'), got {len(percentiles):d}"
        )

    def create_timeseries_plot(self, datasets):
        """Create time series plot."""
        plot_type = "timeseries"
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        logger.info("Plotting %s", plot_type)
        fig = plt.figure(**self.cfg["figure_kwargs"])
        axes = fig.add_subplot()

        # Plot all datasets in one single figure
        ancestors = []
        cubes = {}
        for dataset in datasets:
            ancestors.append(dataset["filename"])
            cube = dataset["cube"]
            cubes[self._get_label(dataset)] = cube
            self._check_cube_dimensions(cube, plot_type)

            # Plot original time series
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs["axes"] = axes
            iris.plot.plot(cube, **plot_kwargs)

            # Plot annual means if desired
            annual_mean_kwargs = self.plots[plot_type]["annual_mean_kwargs"]
            if annual_mean_kwargs is not False:
                logger.debug("Plotting annual means")
                if not cube.coords("year"):
                    add_year(cube, "time")
                annual_mean_cube = cube.aggregated_by(
                    "year", iris.analysis.MEAN
                )
                plot_kwargs.pop("label", None)
                plot_kwargs.update(annual_mean_kwargs)
                iris.plot.plot(annual_mean_cube, **plot_kwargs)

        # Plot horizontal lines
        for hline_kwargs in self.plots[plot_type]["hlines"]:
            axes.axhline(**hline_kwargs)

        # Default plot appearance
        multi_dataset_facets = self._get_multi_dataset_facets(datasets)
        axes.set_title(multi_dataset_facets["long_name"])
        axes.set_xlabel("time")
        # apply time formatting
        if self.plots[plot_type]["time_format"] is not None:
            axes.get_xaxis().set_major_formatter(
                mdates.DateFormatter(self.plots[plot_type]["time_format"])
            )
        axes.set_ylabel(
            f"{multi_dataset_facets[self.cfg['group_variables_by']]} "
            f"[{multi_dataset_facets['units']}]"
        )
        gridline_kwargs = self._get_gridline_kwargs(plot_type)
        if gridline_kwargs is not False:
            axes.grid(**gridline_kwargs)

        # Legend
        legend_kwargs = self.plots[plot_type]["legend_kwargs"]
        if legend_kwargs is not False:
            axes.legend(**legend_kwargs)

        # Customize plot appearance
        self._process_pyplot_kwargs(plot_type, multi_dataset_facets)

        # Save plot
        plot_path = self.get_plot_path(plot_type, multi_dataset_facets)
        fig.savefig(plot_path, **self.cfg["savefig_kwargs"])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save netCDF file
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)
        var_attrs = {
            n: datasets[0][n] for n in ("short_name", "long_name", "units")
        }
        io.save_1d_data(cubes, netcdf_path, "time", var_attrs)

        # Provenance tracking
        caption = (
            f"Time series of {multi_dataset_facets['long_name']} for "
            f"various datasets."
        )
        provenance_record = {
            "ancestors": ancestors,
            "authors": ["schlund_manuel"],
            "caption": caption,
            "plot_types": ["line"],
            "long_names": [var_attrs["long_name"]],
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)
            provenance_logger.log(netcdf_path, provenance_record)

    def create_benchmarking_timeseries(self, datasets):
        """Create time series benchmarking plot."""
        plot_type = "benchmarking_timeseries"
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        logger.info("Plotting %s", plot_type)

        # Get dataset to be benchmarked
        plot_datasets = self._get_benchmark_datasets(datasets)
        # Get percentiles from multi-model statistics
        percentile_dataset = self._get_benchmark_percentiles(datasets)

        fig = plt.figure(**self.cfg["figure_kwargs"])
        axes = fig.add_subplot()

        # load data

        percentile_data = []

        for dataset_to_load in percentile_dataset:
            filename = dataset_to_load["filename"]
            logger.info("Loading %s", filename)
            cube = iris.load_cube(filename)
            percentile_data.append(cube)

        # Plot all datasets in one single figure
        ancestors = []
        cubes = {}

        for dataset in plot_datasets:
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            iris.plot.plot(dataset["cube"], **plot_kwargs)

        yval2 = percentile_dataset[0]["cube"]
        if len(percentile_dataset) > 1:
            idx = len(percentile_dataset) - 1
            yval1 = percentile_dataset[idx]["cube"]
        else:
            yval1 = yval2.copy()
            ymin, __ = axes.get_ylim()
            yval1.data = np.full(len(yval1.data), ymin)

        dataset = plot_datasets[0]
        iris.plot.fill_between(
            dataset["cube"].coord("time"),
            yval1,
            yval2,
            facecolor="lightblue",
            edgecolor="lightblue",
            linewidth=3,
            zorder=1,
            alpha=0.8,
        )

        # Plot horizontal lines
        for hline_kwargs in self.plots[plot_type]["hlines"]:
            axes.axhline(**hline_kwargs)

        # Default plot appearance
        multi_dataset_facets = self._get_multi_dataset_facets(datasets)
        axes.set_title(multi_dataset_facets["long_name"])
        axes.set_xlabel("time")
        # apply time formatting
        if self.plots[plot_type]["time_format"] is not None:
            axes.get_xaxis().set_major_formatter(
                mdates.DateFormatter(self.plots[plot_type]["time_format"])
            )
        axes.set_ylabel(
            f"{multi_dataset_facets[self.cfg['group_variables_by']]} "
            f"[{multi_dataset_facets['units']}]"
        )
        gridline_kwargs = self._get_gridline_kwargs(plot_type)
        if gridline_kwargs is not False:
            axes.grid(**gridline_kwargs)

        # Legend
        legend_kwargs = self.plots[plot_type]["legend_kwargs"]
        if legend_kwargs is not False:
            axes.legend(**legend_kwargs)

        # Customize plot appearance
        self._process_pyplot_kwargs(plot_type, multi_dataset_facets)

        # Save plot
        plot_path = self.get_plot_path(plot_type, multi_dataset_facets)
        fig.savefig(plot_path, **self.cfg["savefig_kwargs"])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save netCDF file
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)
        var_attrs = {
            n: datasets[0][n] for n in ("short_name", "long_name", "units")
        }
        cubes[self._get_label(dataset)] = dataset["cube"]
        io.save_1d_data(cubes, netcdf_path, "time", var_attrs)

        # Provenance tracking
        caption = (
            f"Time series of {multi_dataset_facets['long_name']} for "
            f"various datasets."
        )
        provenance_record = {
            "ancestors": ancestors,
            "authors": ["schlund_manuel"],
            "caption": caption,
            "plot_types": ["line"],
            "long_names": [var_attrs["long_name"]],
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)
            provenance_logger.log(netcdf_path, provenance_record)

    def create_annual_cycle_plot(self, datasets):
        """Create annual cycle plot."""
        plot_type = "annual_cycle"
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        logger.info("Plotting %s", plot_type)
        fig = plt.figure(**self.cfg["figure_kwargs"])
        axes = fig.add_subplot()

        # Plot all datasets in one single figure
        ancestors = []
        cubes = {}
        for dataset in datasets:
            ancestors.append(dataset["filename"])
            cube = dataset["cube"]
            cubes[self._get_label(dataset)] = cube
            self._check_cube_dimensions(cube, plot_type)

            # Plot annual cycle
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs["axes"] = axes
            iris.plot.plot(cube, **plot_kwargs)

        # Plot horizontal lines
        for hline_kwargs in self.plots[plot_type]["hlines"]:
            axes.axhline(**hline_kwargs)

        # Default plot appearance
        multi_dataset_facets = self._get_multi_dataset_facets(datasets)
        axes.set_title(multi_dataset_facets["long_name"])
        axes.set_xlabel("month")
        axes.set_ylabel(
            f"{multi_dataset_facets[self.cfg['group_variables_by']]} "
            f"[{multi_dataset_facets['units']}]"
        )
        axes.set_xticks(range(1, 13), [str(m) for m in range(1, 13)])
        gridline_kwargs = self._get_gridline_kwargs(plot_type)
        if gridline_kwargs is not False:
            axes.grid(**gridline_kwargs)

        # Legend
        legend_kwargs = self.plots[plot_type]["legend_kwargs"]
        if legend_kwargs is not False:
            axes.legend(**legend_kwargs)

        # Customize plot appearance
        self._process_pyplot_kwargs(plot_type, multi_dataset_facets)

        # Save plot
        plot_path = self.get_plot_path(plot_type, multi_dataset_facets)
        fig.savefig(plot_path, **self.cfg["savefig_kwargs"])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save netCDF file
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)
        var_attrs = {
            n: datasets[0][n] for n in ("short_name", "long_name", "units")
        }
        io.save_1d_data(cubes, netcdf_path, "month_number", var_attrs)

        # Provenance tracking
        caption = (
            f"Annual cycle of {multi_dataset_facets['long_name']} for "
            f"various datasets."
        )
        provenance_record = {
            "ancestors": ancestors,
            "authors": ["schlund_manuel"],
            "caption": caption,
            "plot_types": ["seas"],
            "long_names": [var_attrs["long_name"]],
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)
            provenance_logger.log(netcdf_path, provenance_record)

    def create_benchmarking_annual(self, datasets):
        """Create benchmarking annual cycle plot."""
        plot_type = "benchmarking_annual_cycle"
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        logger.info("Plotting %s", plot_type)

        # Get dataset to be benchmarked
        plot_datasets = self._get_benchmark_datasets(datasets)
        # Get percentiles from multi-model statistics
        percentile_dataset = self._get_benchmark_percentiles(datasets)

        fig = plt.figure(**self.cfg["figure_kwargs"])
        axes = fig.add_subplot()

        # Plot all datasets in one single figure
        ancestors = []
        cubes = {}

        # Plot annual cycle(s)
        for dataset in plot_datasets:
            cube = dataset["cube"]
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs["axes"] = axes
            iris.plot.plot(cube, **plot_kwargs)

        yval2 = percentile_dataset[0]["cube"]
        if len(percentile_dataset) > 1:
            idx = len(percentile_dataset) - 1
            yval1 = percentile_dataset[idx]["cube"]
        else:
            yval1 = yval2.copy()
            ymin, __ = axes.get_ylim()
            yval1.data = np.full(len(yval1.data), ymin)

        iris.plot.fill_between(
            cube.coord("month_number"),
            yval1,
            yval2,
            facecolor="lightblue",
            linewidth=0,
            zorder=1,
            alpha=0.8,
        )

        # Plot horizontal lines
        for hline_kwargs in self.plots[plot_type]["hlines"]:
            axes.axhline(**hline_kwargs)

        # Default plot appearance
        multi_dataset_facets = self._get_multi_dataset_facets(datasets)
        axes.set_title(multi_dataset_facets["long_name"])
        axes.set_xlabel("month")
        axes.set_ylabel(
            f"{multi_dataset_facets[self.cfg['group_variables_by']]} "
            f"[{multi_dataset_facets['units']}]"
        )
        axes.set_xticks(range(1, 13), [str(m) for m in range(1, 13)])
        gridline_kwargs = self._get_gridline_kwargs(plot_type)
        if gridline_kwargs is not False:
            axes.grid(**gridline_kwargs)

        # Legend
        legend_kwargs = self.plots[plot_type]["legend_kwargs"]
        if legend_kwargs is not False:
            axes.legend(**legend_kwargs)

        # Customize plot appearance
        self._process_pyplot_kwargs(plot_type, multi_dataset_facets)

        # Save plot
        plot_path = self.get_plot_path(plot_type, multi_dataset_facets)
        fig.savefig(plot_path, **self.cfg["savefig_kwargs"])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save netCDF file
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)
        var_attrs = {
            n: datasets[0][n] for n in ("short_name", "long_name", "units")
        }
        dataset = plot_datasets[0]
        cubes[self._get_label(dataset)] = dataset["cube"]
        io.save_1d_data(cubes, netcdf_path, "month_number", var_attrs)

        # Provenance tracking
        caption = (
            f"Annual cycle of {multi_dataset_facets['long_name']} for "
            f"various datasets."
        )
        provenance_record = {
            "ancestors": ancestors,
            "authors": ["schlund_manuel"],
            "caption": caption,
            "plot_types": ["seas"],
            "long_names": [var_attrs["long_name"]],
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)
            provenance_logger.log(netcdf_path, provenance_record)

    def create_diurnal_cycle_plot(self, datasets):
        """Create diurnal cycle plot."""
        plot_type = "diurnal_cycle"
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        logger.info("Plotting %s", plot_type)
        fig = plt.figure(**self.cfg["figure_kwargs"])
        axes = fig.add_subplot()

        # Plot all datasets in one single figure
        ancestors = []
        cubes = {}
        for dataset in datasets:
            ancestors.append(dataset["filename"])
            cube = dataset["cube"]
            cubes[self._get_label(dataset)] = cube
            self._check_cube_dimensions(cube, plot_type)

            # Plot diurnal cycle
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs["axes"] = axes
            iris.plot.plot(cube, **plot_kwargs)

        # Plot horizontal lines
        for hline_kwargs in self.plots[plot_type]["hlines"]:
            axes.axhline(**hline_kwargs)

        # Default plot appearance
        multi_dataset_facets = self._get_multi_dataset_facets(datasets)
        axes.set_title(multi_dataset_facets["long_name"])
        axes.set_xlabel("Hour")
        axes.set_ylabel(
            f"{multi_dataset_facets[self.cfg['group_variables_by']]} "
            f"[{multi_dataset_facets['units']}]"
        )
        axes.set_xticks(range(0, 24), minor=True)
        axes.set_xticks(range(0, 24, 3), [str(m) for m in range(0, 24, 3)])
        gridline_kwargs = self._get_gridline_kwargs(plot_type)
        if gridline_kwargs is not False:
            axes.grid(**gridline_kwargs)

        # Legend
        legend_kwargs = self.plots[plot_type]["legend_kwargs"]
        if legend_kwargs is not False:
            axes.legend(**legend_kwargs)

        # Customize plot appearance
        self._process_pyplot_kwargs(plot_type, multi_dataset_facets)

        # Save plot
        plot_path = self.get_plot_path(plot_type, multi_dataset_facets)
        fig.savefig(plot_path, **self.cfg["savefig_kwargs"])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save netCDF file
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)
        var_attrs = {
            n: datasets[0][n] for n in ("short_name", "long_name", "units")
        }
        io.save_1d_data(cubes, netcdf_path, "hour", var_attrs)

        # Provenance tracking
        caption = (
            f"Diurnal cycle of {multi_dataset_facets['long_name']} for "
            f"various datasets."
        )
        provenance_record = {
            "ancestors": ancestors,
            "authors": ["schlund_manuel"],
            "caption": caption,
            "plot_types": ["seas"],
            "long_names": [var_attrs["long_name"]],
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)
            provenance_logger.log(netcdf_path, provenance_record)

    def create_benchmarking_diurnal(self, datasets):
        """Create benchmarking diurnal cycle plot."""
        plot_type = "benchmarking_diurnal_cycle"
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        logger.info("Plotting %s", plot_type)

        # Get dataset to be benchmarked
        plot_datasets = self._get_benchmark_datasets(datasets)
        # Get percentiles from multi-model statistics
        percentile_dataset = self._get_benchmark_percentiles(datasets)

        fig = plt.figure(**self.cfg["figure_kwargs"])
        axes = fig.add_subplot()

        # Plot all datasets in one single figure
        ancestors = []
        cubes = {}

        # Plot diurnal cycle(s)
        for dataset in plot_datasets:
            cube = dataset["cube"]
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs["axes"] = axes
            iris.plot.plot(cube, **plot_kwargs)

        yval2 = percentile_dataset[0]["cube"]
        if len(percentile_dataset) > 1:
            idx = len(percentile_dataset) - 1
            yval1 = percentile_dataset[idx]["cube"]
        else:
            yval1 = yval2.copy()
            ymin, __ = axes.get_ylim()
            yval1.data = np.full(len(yval1.data), ymin)

        iris.plot.fill_between(
            cube.coord("hour"),
            yval1,
            yval2,
            facecolor="lightblue",
            linewidth=0,
            zorder=1,
            alpha=0.8,
        )

        # Plot horizontal lines
        for hline_kwargs in self.plots[plot_type]["hlines"]:
            axes.axhline(**hline_kwargs)

        # Default plot appearance
        multi_dataset_facets = self._get_multi_dataset_facets(datasets)
        axes.set_title(multi_dataset_facets["long_name"])
        axes.set_xlabel("Hour")
        axes.set_ylabel(
            f"{multi_dataset_facets[self.cfg['group_variables_by']]} "
            f"[{multi_dataset_facets['units']}]"
        )
        axes.set_xticks(range(0, 24), minor=True)
        axes.set_xticks(range(0, 24, 3), [str(m) for m in range(0, 24, 3)])
        gridline_kwargs = self._get_gridline_kwargs(plot_type)
        if gridline_kwargs is not False:
            axes.grid(**gridline_kwargs)

        # Legend
        legend_kwargs = self.plots[plot_type]["legend_kwargs"]
        if legend_kwargs is not False:
            axes.legend(**legend_kwargs)

        # Customize plot appearance
        self._process_pyplot_kwargs(plot_type, multi_dataset_facets)

        # Save plot
        plot_path = self.get_plot_path(plot_type, multi_dataset_facets)
        fig.savefig(plot_path, **self.cfg["savefig_kwargs"])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save netCDF file
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)
        var_attrs = {
            n: datasets[0][n] for n in ("short_name", "long_name", "units")
        }
        dataset = plot_datasets[0]
        cubes[self._get_label(dataset)] = dataset["cube"]
        io.save_1d_data(cubes, netcdf_path, "hour", var_attrs)

        # Provenance tracking
        caption = (
            f"Diurnal cycle of {multi_dataset_facets['long_name']} for "
            f"various datasets."
        )
        provenance_record = {
            "ancestors": ancestors,
            "authors": ["schlund_manuel"],
            "caption": caption,
            "plot_types": ["seas"],
            "long_names": [var_attrs["long_name"]],
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)
            provenance_logger.log(netcdf_path, provenance_record)

    def create_benchmarking_boxplot(self):
        """Create boxplot."""
        plot_type = "benchmarking_boxplot"
        if plot_type not in self.plots:
            return

        dframe = pd.DataFrame(columns=["Variable", "Dataset", "Value"])
        ifile = 0

        cubes = iris.cube.CubeList()
        benchmark_datasets = []
        variables = []

        for var_key, datasets in self.grouped_input_data.items():
            logger.info("Processing variable %s", var_key)

            if not datasets:
                raise ValueError(f"No input data to plot '{plot_type}' given")

            # Get dataset to be benchmarked
            plot_datasets = self._get_benchmark_datasets(datasets)
            benchmark_dataset = plot_datasets[0]

            logger.info(
                "Plotting %s for dataset %s",
                plot_type,
                benchmark_dataset["dataset"],
            )

            # Get datasets for benchmarking
            benchmark_group = self._get_benchmark_group(datasets)
            logger.info(
                "Benchmarking group of %i datasets.", len(benchmark_group)
            )

            ancestors = [benchmark_dataset["filename"]]
            for dataset in benchmark_group:
                ancestors.append(dataset["filename"])

            for dataset in benchmark_group:
                dataset_name = dataset["dataset"]
                cube = iris.load_cube(dataset["filename"])
                dframe.loc[ifile] = [var_key, dataset_name, cube.data]
                ifile = ifile + 1

            dframe["Value"] = dframe["Value"].astype(str).astype(float)

            cubes.append(benchmark_dataset["cube"])
            benchmark_datasets.append(benchmark_dataset)
            variables.append(var_key)

        # order of variables
        if self.plots[plot_type]["var_order"]:
            var_order = self.plots[plot_type]["var_order"]
            if set(variables) == set(var_order):
                ind = [
                    variables.index(var_order[i])
                    for i in range(len(variables))
                ]
                cubes = iris.cube.CubeList([cubes[i] for i in ind])
                benchmark_datasets = [benchmark_datasets[i] for i in ind]
                variables = var_order
            else:
                raise ValueError(
                    "List of ordered variables do not agree with"
                    " processed variables"
                )

        (plot_path, netcdf_paths) = self._plot_benchmarking_boxplot(
            dframe, cubes, variables, benchmark_datasets
        )

        # Save plot
        plt.savefig(plot_path, **self.cfg["savefig_kwargs"])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save netCDF file
        for netcdf_path, cube in netcdf_paths.items():
            io.iris_save(cube, netcdf_path)

            # Provenance tracking
            caption = "Boxplot."
            provenance_record = {
                "ancestors": ancestors,
                "authors": ["bock_lisa", "schlund_manuel"],
                "caption": caption,
                "plot_types": ["box"],
            }
            with ProvenanceLogger(self.cfg) as provenance_logger:
                provenance_logger.log(plot_path, provenance_record)
                provenance_logger.log(netcdf_path, provenance_record)

    def create_map_plot(self, datasets):
        """Create map plot."""
        plot_type = "map"
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        # Get reference dataset if possible
        ref_dataset = self._get_reference_dataset(datasets)
        if ref_dataset is None:
            logger.info("Plotting %s without reference dataset", plot_type)
        else:
            logger.info(
                "Plotting %s with reference dataset '%s'",
                plot_type,
                self._get_label(ref_dataset),
            )

        # Get plot function
        plot_func = self._get_plot_func(plot_type)

        # Create a single plot for each dataset (incl. reference dataset if
        # given)
        for dataset in datasets:
            if dataset == ref_dataset:
                continue
            ancestors = [dataset["filename"]]
            if ref_dataset is None:
                (plot_path, netcdf_paths) = self._plot_map_without_ref(
                    plot_func, dataset
                )
                caption = (
                    f"Map plot of {dataset['long_name']} of dataset "
                    f"{dataset['alias']}."
                )
            else:
                (plot_path, netcdf_paths) = self._plot_map_with_ref(
                    plot_func, dataset, ref_dataset
                )
                caption = (
                    f"Map plot of {dataset['long_name']} of dataset "
                    f"{dataset['alias']} including bias relative to "
                    f"{ref_dataset['alias']}."
                )
                ancestors.append(ref_dataset["filename"])

            # If statistics are shown add a brief description to the caption
            if self.plots[plot_type]["show_stats"]:
                caption += (
                    " The number in the top left corner corresponds to the "
                    "spatial mean (weighted by grid cell areas)."
                )

            # Save plot
            plt.savefig(plot_path, **self.cfg["savefig_kwargs"])
            logger.info("Wrote %s", plot_path)
            plt.close()

            # Save netCDFs
            for netcdf_path, cube in netcdf_paths.items():
                io.iris_save(cube, netcdf_path)

            # Provenance tracking
            provenance_record = {
                "ancestors": ancestors,
                "authors": ["schlund_manuel"],
                "caption": caption,
                "plot_types": ["map"],
                "long_names": [dataset["long_name"]],
            }
            with ProvenanceLogger(self.cfg) as provenance_logger:
                provenance_logger.log(plot_path, provenance_record)
                for netcdf_path in netcdf_paths:
                    provenance_logger.log(netcdf_path, provenance_record)

    def create_benchmarking_map_plot(self, datasets):
        """Create benchmarking map plot."""
        plot_type = "benchmarking_map"
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        # Get reference dataset
        ref_dataset = self._get_benchmarking_reference(datasets)
        # Get dataset to be benchmarked
        plot_datasets = self._get_benchmark_datasets(datasets)
        # Get percentiles from multi-model statistics
        percentile_dataset = self._get_benchmark_percentiles(datasets)
        # Get benchmarking metric
        metric = self._get_benchmark_metric(datasets)

        # Get plot function
        plot_func = self._get_plot_func(plot_type)

        # load data

        percentile_data = []

        for dataset_to_load in percentile_dataset:
            filename = dataset_to_load["filename"]
            logger.info("Loading %s", filename)
            cube = iris.load_cube(filename)
            percentile_data.append(cube)

        for dataset in plot_datasets:
            ancestors = [dataset["filename"]]
            (plot_path, netcdf_paths) = self._plot_benchmarking_map(
                plot_func, dataset, percentile_data, metric
            )
            caption = (
                f"Map plot of {dataset['long_name']} of dataset "
                f"{dataset['alias']}."
            )
            ancestors.append(ref_dataset["filename"])

            # Save plot
            plt.savefig(plot_path, **self.cfg["savefig_kwargs"])
            logger.info("Wrote %s", plot_path)
            plt.close()

            # Save netCDFs
            for netcdf_path, cube in netcdf_paths.items():
                io.iris_save(cube, netcdf_path)

            # Provenance tracking
            provenance_record = {
                "ancestors": ancestors,
                "authors": ["schlund_manuel"],
                "caption": caption,
                "plot_types": ["map"],
                "long_names": [dataset["long_name"]],
            }
            with ProvenanceLogger(self.cfg) as provenance_logger:
                provenance_logger.log(plot_path, provenance_record)
                for netcdf_path in netcdf_paths:
                    provenance_logger.log(netcdf_path, provenance_record)

    def create_zonal_mean_profile_plot(self, datasets):
        """Create zonal mean profile plot."""
        plot_type = "zonal_mean_profile"
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        # Get reference dataset if possible
        ref_dataset = self._get_reference_dataset(datasets)
        if ref_dataset is None:
            logger.info("Plotting %s without reference dataset", plot_type)
        else:
            logger.info(
                "Plotting %s with reference dataset '%s'",
                plot_type,
                self._get_label(ref_dataset),
            )

        # Get plot function
        plot_func = self._get_plot_func(plot_type)

        # Create a single plot for each dataset (incl. reference dataset if
        # given)
        for dataset in datasets:
            if dataset == ref_dataset:
                continue
            ancestors = [dataset["filename"]]
            if ref_dataset is None:
                (plot_path, netcdf_paths) = (
                    self._plot_zonal_mean_profile_without_ref(
                        plot_func, dataset
                    )
                )
                caption = (
                    f"Zonal mean profile of {dataset['long_name']} of dataset "
                    f"{dataset['alias']}."
                )
            else:
                (plot_path, netcdf_paths) = (
                    self._plot_zonal_mean_profile_with_ref(
                        plot_func, dataset, ref_dataset
                    )
                )
                caption = (
                    f"Zonal mean profile of {dataset['long_name']} of dataset "
                    f"{dataset['alias']} including bias relative to "
                    f"{ref_dataset['alias']}."
                )
                ancestors.append(ref_dataset["filename"])

            # If statistics are shown add a brief description to the caption
            if self.plots[plot_type]["show_stats"]:
                caption += (
                    " The number in the top left corner corresponds to the "
                    "spatial mean (weighted by grid cell areas)."
                )

            # Save plot
            plt.savefig(plot_path, **self.cfg["savefig_kwargs"])
            logger.info("Wrote %s", plot_path)
            plt.close()

            # Save netCDFs
            for netcdf_path, cube in netcdf_paths.items():
                io.iris_save(cube, netcdf_path)

            # Provenance tracking
            provenance_record = {
                "ancestors": ancestors,
                "authors": ["schlund_manuel"],
                "caption": caption,
                "plot_types": ["vert"],
                "long_names": [dataset["long_name"]],
            }
            with ProvenanceLogger(self.cfg) as provenance_logger:
                provenance_logger.log(plot_path, provenance_record)
                for netcdf_path in netcdf_paths:
                    provenance_logger.log(netcdf_path, provenance_record)

    def create_benchmarking_zonal_plot(self, datasets):
        """Create benchmarking zonal mean profile plot."""
        plot_type = "benchmarking_zonal"
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        # Get dataset to be benchmarked
        plot_datasets = self._get_benchmark_datasets(datasets)
        # Get percentiles from multi-model statistics
        percentile_dataset = self._get_benchmark_percentiles(datasets)
        # Get benchmarking metric
        metric = self._get_benchmark_metric(datasets)

        # Get plot function
        plot_func = self._get_plot_func(plot_type)

        # Create a single plot for each dataset (incl. reference dataset if
        # given)

        # load data

        percentile_data = []

        for dataset_to_load in percentile_dataset:
            filename = dataset_to_load["filename"]
            logger.info("Loading %s", filename)
            cube = iris.load_cube(filename)
            percentile_data.append(cube)

        for dataset in plot_datasets:
            (plot_path, netcdf_paths) = self._plot_benchmarking_zonal(
                plot_func, dataset, percentile_data, metric
            )
            ancestors = [dataset["filename"]]

            caption = (
                f"Zonal mean profile of {dataset['long_name']} of dataset "
                f"{dataset['alias']}."
            )

            # Save plot
            plt.savefig(plot_path, **self.cfg["savefig_kwargs"])
            logger.info("Wrote %s", plot_path)
            plt.close()

            # Save netCDFs
            for netcdf_path, cube in netcdf_paths.items():
                io.iris_save(cube, netcdf_path)

            # Provenance tracking
            provenance_record = {
                "ancestors": ancestors,
                "authors": ["schlund_manuel"],
                "caption": caption,
                "plot_types": ["vert"],
                "long_names": [dataset["long_name"]],
            }
            with ProvenanceLogger(self.cfg) as provenance_logger:
                provenance_logger.log(plot_path, provenance_record)
                for netcdf_path in netcdf_paths:
                    provenance_logger.log(netcdf_path, provenance_record)

    def create_1d_profile_plot(self, datasets):
        """Create 1D profile plot."""
        plot_type = "1d_profile"
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        logger.info("Plotting %s", plot_type)
        fig = plt.figure(**self.cfg["figure_kwargs"])
        axes = fig.add_subplot()

        multi_dataset_facets = self._get_multi_dataset_facets(datasets)

        # Plot all datasets in one single figure
        ancestors = []
        cubes = {}
        for dataset in datasets:
            ancestors.append(dataset["filename"])
            cube = dataset["cube"]
            cubes[self._get_label(dataset)] = cube
            self._check_cube_dimensions(cube, plot_type)

            # Plot 1D profile
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs["axes"] = axes

            iris.plot.plot(cube, **plot_kwargs)

        # Plot horizontal lines
        for hline_kwargs in self.plots[plot_type]["hlines"]:
            axes.axhline(**hline_kwargs)

        # Default plot appearance
        axes.set_title(multi_dataset_facets["long_name"])
        axes.set_xlabel(
            f"{multi_dataset_facets[self.cfg['group_variables_by']]} "
            f"[{multi_dataset_facets['units']}]"
        )
        z_coord = cube.coord(axis="Z")
        axes.set_ylabel(f"{z_coord.long_name} [{z_coord.units}]")

        # apply logarithmic axes
        if self.plots[plot_type]["log_y"]:
            axes.set_yscale("log")
            axes.get_yaxis().set_major_formatter(FormatStrFormatter("%.1f"))
        if self.plots[plot_type]["show_y_minor_ticklabels"]:
            axes.get_yaxis().set_minor_formatter(FormatStrFormatter("%.1f"))
        else:
            axes.get_yaxis().set_minor_formatter(NullFormatter())
        if self.plots[plot_type]["log_x"]:
            axes.set_xscale("log")
            # major and minor ticks
            x_major = LogLocator(base=10.0, numticks=12)
            axes.get_xaxis().set_major_locator(x_major)
            x_minor = LogLocator(
                base=10.0, subs=np.arange(1.0, 10.0) * 0.1, numticks=12
            )

            axes.get_xaxis().set_minor_locator(x_minor)
            axes.get_xaxis().set_minor_formatter(NullFormatter())

        # gridlines
        gridline_kwargs = self._get_gridline_kwargs(plot_type)
        if gridline_kwargs is not False:
            axes.grid(**gridline_kwargs)
        # nicer aspect ratio
        aspect_ratio = self.plots[plot_type]["aspect_ratio"]
        axes.set_box_aspect(aspect_ratio)

        # Legend
        legend_kwargs = self.plots[plot_type]["legend_kwargs"]
        if legend_kwargs is not False:
            axes.legend(**legend_kwargs)

        # Customize plot appearance
        self._process_pyplot_kwargs(plot_type, multi_dataset_facets)

        # Save plot
        plot_path = self.get_plot_path(plot_type, multi_dataset_facets)
        fig.savefig(plot_path, **self.cfg["savefig_kwargs"])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save netCDF file
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)
        var_attrs = {
            n: datasets[0][n] for n in ("short_name", "long_name", "units")
        }
        io.save_1d_data(cubes, netcdf_path, z_coord.standard_name, var_attrs)

        # Provenance tracking
        caption = (
            "Vertical one-dimensional profile of "
            f"{multi_dataset_facets['long_name']}"
            " for various datasets."
        )
        provenance_record = {
            "ancestors": ancestors,
            "authors": ["schlund_manuel", "winterstein_franziska"],
            "caption": caption,
            "plot_types": ["line"],
            "long_names": [var_attrs["long_name"]],
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)
            provenance_logger.log(netcdf_path, provenance_record)

    def create_variable_vs_lat_plot(self, datasets):
        """Create Variable as a function of latitude."""
        plot_type = "variable_vs_lat"
        if plot_type not in self.plots:
            return
        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")
        logger.info("Plotting %s", plot_type)
        fig = plt.figure(**self.cfg["figure_kwargs"])
        axes = fig.add_subplot()

        # Plot all datasets in one single figure
        ancestors = []
        cubes = {}
        for dataset in datasets:
            ancestors.append(dataset["filename"])
            cube = dataset["cube"]
            cubes[self._get_label(dataset)] = cube
            self._check_cube_dimensions(cube, plot_type)

            # Plot data
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            plot_kwargs["axes"] = axes
            iris.plot.plot(cube, **plot_kwargs)

        # Plot horizontal lines
        for hline_kwargs in self.plots[plot_type]["hlines"]:
            axes.axhline(**hline_kwargs)

        # Default plot appearance
        multi_dataset_facets = self._get_multi_dataset_facets(datasets)
        axes.set_title(multi_dataset_facets["long_name"])
        axes.set_xlabel("latitude [°N]")
        axes.set_ylabel(
            f"{multi_dataset_facets[self.cfg['group_variables_by']]} "
            f"[{multi_dataset_facets['units']}]"
        )
        gridline_kwargs = self._get_gridline_kwargs(plot_type)
        if gridline_kwargs is not False:
            axes.grid(**gridline_kwargs)

        # Legend
        legend_kwargs = self.plots[plot_type]["legend_kwargs"]
        if legend_kwargs is not False:
            axes.legend(**legend_kwargs)

        # Customize plot appearance
        self._process_pyplot_kwargs(plot_type, multi_dataset_facets)

        # Save plot
        plot_path = self.get_plot_path(plot_type, multi_dataset_facets)
        fig.savefig(plot_path, **self.cfg["savefig_kwargs"])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save netCDF file
        netcdf_path = get_diagnostic_filename(Path(plot_path).stem, self.cfg)
        var_attrs = {
            n: datasets[0][n] for n in ("short_name", "long_name", "units")
        }
        io.save_1d_data(cubes, netcdf_path, "latitude", var_attrs)

        # Provenance tracking
        caption = (
            f"{multi_dataset_facets['long_name']} vs. latitude for "
            f"various datasets."
        )
        provenance_record = {
            "ancestors": ancestors,
            "authors": ["sarauer_ellen"],
            "caption": caption,
            "plot_types": ["line"],
            "long_names": [var_attrs["long_name"]],
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)
            provenance_logger.log(netcdf_path, provenance_record)

    def create_hovmoeller_z_vs_time_plot(self, datasets):
        """Create Hovmoeller Z vs. time plot."""
        plot_type = "hovmoeller_z_vs_time"
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        # Get reference dataset if possible
        ref_dataset = self._get_reference_dataset(datasets)
        if ref_dataset is None:
            logger.info("Plotting %s without reference dataset", plot_type)
        else:
            logger.info(
                "Plotting %s with reference dataset '%s'",
                plot_type,
                self._get_label(ref_dataset),
            )

        # Get plot function
        plot_func = self._get_plot_func(plot_type)

        # Create a single plot for each dataset (incl. reference dataset if
        # given)
        for dataset in datasets:
            if dataset == ref_dataset:
                continue
            ancestors = [dataset["filename"]]
            if ref_dataset is None:
                (plot_path, netcdf_paths) = (
                    self._plot_hovmoeller_z_vs_time_without_ref(
                        plot_func, dataset
                    )
                )
                caption = (
                    f"Hovmoeller Z vs. time plot of {dataset['long_name']} "
                    f"of dataset {dataset['alias']}."
                )
            else:
                (plot_path, netcdf_paths) = (
                    self._plot_hovmoeller_z_vs_time_with_ref(
                        plot_func, dataset, ref_dataset
                    )
                )
                caption = (
                    f"Hovmoeller Z vs. time plot of {dataset['long_name']} "
                    f"of dataset {dataset['alias']} including bias relative "
                    f"to {ref_dataset['alias']}."
                )
                ancestors.append(ref_dataset["filename"])

            # If statistics are shown add a brief description to the caption
            if self.plots[plot_type]["show_stats"]:
                caption += (
                    " The number in the top left corner corresponds to the "
                    "spatiotemporal mean."
                )

            # Save plot
            plt.savefig(plot_path, **self.cfg["savefig_kwargs"])
            logger.info("Wrote %s", plot_path)
            plt.close()

            # Save netCDFs
            for netcdf_path, cube in netcdf_paths.items():
                io.iris_save(cube, netcdf_path)

            # Provenance tracking
            provenance_record = {
                "ancestors": ancestors,
                "authors": ["kuehbacher_birgit", "heuer_helge"],
                "caption": caption,
                "plot_types": ["vert"],
                "long_names": [dataset["long_name"]],
            }
            with ProvenanceLogger(self.cfg) as provenance_logger:
                provenance_logger.log(plot_path, provenance_record)
                for netcdf_path in netcdf_paths:
                    provenance_logger.log(netcdf_path, provenance_record)

    def create_hovmoeller_time_vs_lat_or_lon_plot(self, datasets):
        """Create the Hovmoeller plot with time vs latitude or longitude."""
        plot_type = "hovmoeller_time_vs_lat_or_lon"
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        # Get reference dataset if possible
        ref_dataset = self._get_reference_dataset(datasets)
        if ref_dataset is None:
            logger.info("Plotting %s without reference dataset", plot_type)
        else:
            logger.info(
                "Plotting %s with reference dataset '%s'",
                plot_type,
                self._get_label(ref_dataset),
            )

        # Get plot function
        plot_func = self._get_plot_func(plot_type)

        # Create a single plot for each dataset (incl. reference dataset if
        # given)
        for dataset in datasets:
            if dataset == ref_dataset:
                continue
            ancestors = [dataset["filename"]]
            if ref_dataset is None:
                (plot_path, netcdf_paths) = (
                    self._plot_hovmoeller_time_vs_lat_or_lon_without_ref(
                        plot_func, dataset
                    )
                )
                caption = (
                    f"Hovmoeller plot of {dataset['long_name']} of dataset "
                    f"{dataset['alias']}."
                )
            else:
                (plot_path, netcdf_paths) = (
                    self._plot_hovmoeller_time_vs_lat_or_lon_with_ref(
                        plot_func, dataset, ref_dataset
                    )
                )
                caption = (
                    f"Hovmoeller plot of {dataset['long_name']} of dataset "
                    f"{dataset['alias']} including bias relative to "
                    f"{ref_dataset['alias']}."
                )
                ancestors.append(ref_dataset["filename"])

            # Save plot
            plt.savefig(plot_path, **self.cfg["savefig_kwargs"])
            logger.info("Wrote %s", plot_path)
            plt.close()

            # Save netCDFs
            for netcdf_path, cube in netcdf_paths.items():
                io.iris_save(cube, netcdf_path)

            # Provenance tracking
            provenance_record = {
                "ancestors": ancestors,
                "authors": [
                    "schlund_manuel",
                    "kraft_jeremy",
                    "lindenlaub_lukas",
                ],
                "caption": caption,
                "plot_types": ["zonal"],
                "long_names": [dataset["long_name"]],
            }
            with ProvenanceLogger(self.cfg) as provenance_logger:
                provenance_logger.log(plot_path, provenance_record)
                for netcdf_path in netcdf_paths:
                    provenance_logger.log(netcdf_path, provenance_record)

    def create_hovmoeller_anncyc_vs_lat_or_lon_plot(self, datasets):
        """Create the Hovmoeller plot with annual cycle vs latitude or longitude."""
        plot_type = "hovmoeller_anncyc_vs_lat_or_lon"
        if plot_type not in self.plots:
            return

        if not datasets:
            raise ValueError(f"No input data to plot '{plot_type}' given")

        # Get reference dataset if possible
        ref_dataset = self._get_reference_dataset(datasets)
        if ref_dataset is None:
            logger.info("Plotting %s without reference dataset", plot_type)
        else:
            logger.info(
                "Plotting %s with reference dataset '%s'",
                plot_type,
                self._get_label(ref_dataset),
            )

        # Get plot function
        plot_func = self._get_plot_func(plot_type)

        # Create a single plot for each dataset (incl. reference dataset if
        # given)
        for dataset in datasets:
            if dataset == ref_dataset:
                continue
            ancestors = [dataset["filename"]]
            if ref_dataset is None:
                (plot_path, netcdf_paths) = (
                    self._plot_hovmoeller_anncyc_vs_lat_or_lon_without_ref(
                        plot_func, dataset
                    )
                )
                caption = (
                    f"Hovmoeller plot of {dataset['long_name']} of dataset "
                    f"{dataset['alias']}."
                )
            else:
                (plot_path, netcdf_paths) = (
                    self._plot_hovmoeller_anncyc_vs_lat_or_lon_with_ref(
                        plot_func, dataset, ref_dataset
                    )
                )
                caption = (
                    f"Hovmoeller plot of {dataset['long_name']} of dataset "
                    f"{dataset['alias']} including bias relative to "
                    f"{ref_dataset['alias']}."
                )
                ancestors.append(ref_dataset["filename"])

            # Save plot
            plt.savefig(plot_path, **self.cfg["savefig_kwargs"])
            logger.info("Wrote %s", plot_path)
            plt.close()

            # Save netCDFs
            for netcdf_path, cube in netcdf_paths.items():
                io.iris_save(cube, netcdf_path)

            # Provenance tracking
            provenance_record = {
                "ancestors": ancestors,
                "authors": [
                    "schlund_manuel",
                    "hassler_birgit",
                ],
                "caption": caption,
                "plot_types": ["zonal"],
                "long_names": [dataset["long_name"]],
            }
            with ProvenanceLogger(self.cfg) as provenance_logger:
                provenance_logger.log(plot_path, provenance_record)
                for netcdf_path in netcdf_paths:
                    provenance_logger.log(netcdf_path, provenance_record)

    def compute(self):
        """Plot preprocessed data."""
        with mpl.rc_context(self.cfg["matplotlib_rc_params"]):
            self.create_benchmarking_boxplot()
            for var_key, datasets in self.grouped_input_data.items():
                logger.info("Processing variable %s", var_key)
                self.create_timeseries_plot(datasets)
                self.create_annual_cycle_plot(datasets)
                self.create_benchmarking_annual(datasets)
                self.create_benchmarking_diurnal(datasets)
                self.create_benchmarking_map_plot(datasets)
                self.create_benchmarking_timeseries(datasets)
                self.create_benchmarking_zonal_plot(datasets)
                self.create_diurnal_cycle_plot(datasets)
                self.create_map_plot(datasets)
                self.create_zonal_mean_profile_plot(datasets)
                self.create_1d_profile_plot(datasets)
                self.create_variable_vs_lat_plot(datasets)
                self.create_hovmoeller_z_vs_time_plot(datasets)
                self.create_hovmoeller_time_vs_lat_or_lon_plot(datasets)
                self.create_hovmoeller_anncyc_vs_lat_or_lon_plot(datasets)


def main():
    """Run diagnostic."""
    with run_diagnostic() as config:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Using DEFAULT_SPHERICAL_EARTH_RADIUS",
                category=UserWarning,
                module="iris",
            )
            MultiDatasets(config).compute()


if __name__ == "__main__":
    main()
