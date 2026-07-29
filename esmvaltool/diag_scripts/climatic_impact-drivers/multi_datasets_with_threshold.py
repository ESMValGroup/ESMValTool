"""Diagnostic vizualizing climate impact drivers.

Description
-----------
This diagnostic vizualizes climate impact drivers that were defined
in Elling et al. (2026) for various models.

This diagnostic adapts the monitor/multi_datasets.py diagnostic adding
pre-processing routines to allow for plots of the number of days in a
year that certain thresholds are exceeded. With this diagnostic,
multiple datasets can be visualized in a single plot.

Plot types can be specified with the recipe option 'plots' and pre-
processing options can be accessed with the recipe option 'options'.

Supported pre-processing options
--------------------------------
-   ``threshold_conversion``: Replace the given dataset by the count of
                              on how many days the data exceeds a
                              certian threshold at some point of time.

    Additional options
    ------------------
    threshold: float
        The threshold which should be exceeded for days to add to the
        count.
    inverted: bool, optional
        Optional boolean that can be set to true to count the days on
        which the values drop below the threshold at some point of time
        instead. The default value is ``false``.
    accumulated: bool, optional
        Optional boolean that must be set to true to check if the daily
        accumulated amount is exceeding the threshold for variables
        that are accumulating over time (e.g. absolute amount of
        precipitation or snow) if the temporal resolution is shorter
        than daily. For all other variables (i.e. variables for which
        the number of days on which the absolute value of the variable
        exceeds the threshold at some point in time should be counted)
        use the default value. The default value is ``false``.

Supported plot types
--------------------
-   ``timeseries`` (1D plot): Plot time series. Input data needs single
                              dimension `time`. For each variable
                              separately, all datasets are plotted in
                              one single figure. Input data needs to
                              be 1D.
-   ``map`` (2D plot): Plot map plot. Input data needs dimensions
                       `(longitude, latitude)`. For each variable and
                       dataset, an individual figure is plotted. Input
                       data needs to be 2D. A single reference dataset
                       can be defined by setting the facet
                       ``reference_for_monitor_diags: True`` in the
                       dataset definition in the recipe. In this case,
                       three panels are plotted, incl. a bias. Note
                       that if a reference dataset is defined, all
                       input datasets need to be given on the same
                       horizontal and vertical grid (you can use the
                       preprocessors :func:`esmvalcore.preprocessor.regrid`
                       and :func:`esmvalcore.preprocessor.extract_levels`
                       for this).

    Additional options for timeseries
    ---------------------------------
    aspect_ratio: float, optional (default: None)
        Aspect ratio of the plot.
    axes_kwargs: dict, optional
        Optional calls to methods of the corresponding
        :class:`matplotlib.axes.Axes` instance. Dictionary keys are
        functions of :class:`matplotlib.axes.Axes`. Dictionary values
        are used as argument(s) for these functions (if values are
        dictionaries, these are interpreted as keyword arguments;
        otherwise a single argument is assumed). String arguments can
        include facets in curly brackets which will be derived from the
        corresponding dataset, e.g., ``{project}``, ``{short_name}``,
        ``{exp}``. Examples: ``{set_title: 'Awesome Plot of {long_name}'}``,
        ``{set_xlabel: '{short_name}'}``, ``{set_xlim: [0, 5]}``.
    caption: str, optional
        Figure caption used for provenance tracking. Can include facets
        in curly brackets which will be derived from the corresponding
        dataset, e.g., ``{project}``, ``{short_name}``, ``{exp}``.
        By default, uses a very basic caption.
    gridline_kwargs: dict, optional
        Optional keyword arguments for grid lines. By default, uses
        ``{color: 'lightgrey', alpha: 0.5}``. Use ``gridline_kwargs:
        False`` to not show gridlines.
    hlines: list of dict, optional
        Horizontal lines to show in plot. Each list element corresponds
        to one line, and each list element should contain a dictionary
        with keywords arguments passed to
        :meth:`matplotlib.axes.Axes.axhline`.
        Example: ``[{y: 0}, {y: 1, color: 'red'}]``.
    legend_kwargs: dict, optional
        Optional keyword arguments for :func:`matplotlib.pyplot.legend`.
        Use ``legend_kwargs: False`` to not show legends.
    log_x: bool, optional (default: False)
        Use logarithmic X-axis.
    log_y: bool, optional (default: False)
        Use logarithmic Y-axis.
    plot_kwargs: dict, optional
        Optional keyword arguments for :func:`iris.plot.plot`.
        Dictionary keys are elements identified by
        ``facet_used_for_labels`` or ``'default'``, e.g., ``'CMIP6'``
        if ``facet_used_for_labels: 'project'`` or ``'historical'`` if
        ``facet_used_for_labels: 'exp'``. Dictionary values are
        dictionaries used as keyword arguments for
        :func:`iris.plot.plot`. String arguments can include facets in
        curly brackets which will be derived from the corresponding
        dataset, e.g., ``{project}``, ``{short_name}``, ``{exp}``.
        Examples: ``{default: {linestyle: '-', label: '{project}'},
        CMIP6: {color: 'red', linestyle: '--'}, OBS: {color: 'black'}}``.
    pyplot_kwargs: dict, optional
        Optional calls to functions of :mod:`matplotlib.pyplot`.
        Dictionary keys are functions of :mod:`matplotlib.pyplot`.
        Dictionary values are used as argument(s) for these functions
        (if values are dictionaries, these are interpreted as keyword
        arguments; otherwise a single argument is assumed). String
        arguments can include facets in curly brackets which will be
        derived from the corresponding dataset, e.g., ``{project}``,
        ``{short_name}``, ``{exp}``.
        Examples: ``{title: 'Awesome Plot of {long_name}'}``,
        ``{xlabel: '{short_name}'}``, ``{xlim: [0, 5]}``.
    rasterize: bool, optional (default: False)
        If ``True``, use rasterization_ for plots to produce smaller
        files.  This is only relevant for vector graphics
        (e.g., ``output_file_type: 'pdf'``).
    time_format: str, optional (default: None)
        :func:`~datetime.datetime.strftime` format string that is used
        to format the time axis using :class:`matplotlib.dates.DateFormatter`.
        If ``None``, use the default formatting imposed by the iris
        plotting function.
    transpose_axes: bool, optional (default: False)
        Swap X- and Y-axis.
    x_major_formatter: str, optional (default: None)
        Format string for :class:`matplotlib.ticker.FormatStrFormatter`
        used to format major tick labels of X-axis.
    x_minor_formatter: str, optional (default: None)
        Format string for :class:`matplotlib.ticker.FormatStrFormatter`
        used to format minor tick labels of X-axis.
    y_major_formatter: str, optional (default: None)
        Format string for :class:`matplotlib.ticker.FormatStrFormatter`
        used to format major tick labels of Y-axis.
    y_minor_formatter: str, optional (default: None)
        Format string for :class:`matplotlib.ticker.FormatStrFormatter`
        used to format minor tick labels of Y-axis.

    Additional options for map
    --------------------------
    aspect_ratio: float, optional (default: None)
        Aspect ratio of the plot.
    axes_kwargs: dict, optional
        Optional calls to methods of the corresponding
        :class:`matplotlib.axes.Axes` instance. Dictionary keys are
        functions of :class:`matplotlib.axes.Axes`. Dictionary values
        are used as argument(s) for these functions (if values are
        dictionaries, these are interpreted as keyword arguments;
        otherwise a single argument is assumed). String arguments can
        include facets in curly brackets which will be derived from the
        corresponding dataset, e.g., ``{project}``, ``{short_name}``,
        ``{exp}``. Examples: ``{set_title: 'Plot of {long_name}'}``,
        ``{set_xlabel: '{short_name}'}``, ``{set_xlim: [0, 5]}``.
    caption: str, optional
        Figure caption used for provenance tracking. Can include facets
        in curly brackets which will be derived from the corresponding
        dataset, e.g., ``{project}``, ``{short_name}``, ``{exp}``.
        By default, uses a very basic caption.
    cbar_label: str, optional (default: '{short_name} [{units}]')
        Colorbar label. Can include facets in curly brackets which will
        be derived from the corresponding dataset, e.g., ``{project}``,
        ``{short_name}``, ``{exp}``.
    cbar_label_bias: str, optional (default: 'Δ{short_name} [{units}]')
        Colorbar label for plotting biases. Can include facets in curly
        brackets which will be derived from the corresponding dataset,
        e.g., ``{project}``, ``{short_name}``, ``{exp}``. Only relevant
        for plots including reference datasets.
    cbar_kwargs: dict, optional
        Optional keyword arguments for :func:`matplotlib.pyplot.colorbar`.
        By default, uses ``{orientation: 'vertical'}``.
    cbar_kwargs_bias: dict, optional
        Optional keyword arguments for :func:`matplotlib.pyplot.colorbar`
        for plotting biases. These keyword arguments update (and
        potentially overwrite) the ``cbar_kwargs`` for the bias plot.
        Only relevant for plots including reference datasets.
    common_cbar: bool, optional (default: False)
        Use a common colorbar for the top panels (i.e., plots of the
        dataset and the corresponding reference dataset) when using a
        reference dataset. If neither ``vmin`` and ``vmax`` nor
        ``levels`` is given in ``plot_kwargs``, the colorbar bounds are
        inferred from the dataset in the top left panel, which might
        lead to an inappropriate colorbar for the reference dataset
        (top right panel). Thus, specify ``vmin`` and ``vmax`` or
        ``levels`` when using ``common_cbar: True``. Only relevant for
        plots including reference datasets.
    fontsize: int, optional (default: None)
        Fontsize used for ticks, labels and titles. For the latter, use
        the given fontsize plus 2. Does not affect suptitles. If not
        given, use default matplotlib values. For a more fine-grained
        definition of fontsizes, use the option ``matplotlib_rc_params``
        (see above).
    gridline_kwargs: dict, optional (default: False)
        Optional keyword arguments for grid lines. Use
        ``gridline_kwargs: False`` to not show grid lines.
    log_x: bool, optional (default: False)
        Use logarithmic X-axis.
    log_y: bool, optional (default: False)
        Use logarithmic Y-axis.
    plot_func: str, optional (default: 'contourf')
        Plot function used to plot the maps. Must be a function of
        :mod:`iris.plot` that supports plotting of 2D data.
    plot_kwargs: dict, optional
        Optional keyword arguments for the plot function defined by
        ``plot_func``. Dictionary keys are elements identified by
        ``facet_used_for_labels`` or ``'default'``, e.g., ``'CMIP6'``
        if ``facet_used_for_labels: 'project'`` or ``'historical'`` if
        ``facet_used_for_labels: 'exp'``. Dictionary values are
        dictionaries used as keyword arguments for the plot function
        defined by ``plot_func``. String arguments can include facets
        in curly brackets which will be derived from the corresponding
        dataset, e.g., ``{project}``, ``{short_name}``, ``{exp}``.
        Examples: ``{default: {levels: 2}, CMIP6: {vmin: 200, vmax: 250}}``.
        In addition to the normalization_ options supported by the plot
        function, the option ``{norm: 'centered'}`` can be specified.
        In this case, the keywords ``vcenter`` and ``halfrange`` should
        be used instead of ``vmin`` or ``vmax`` (see
        :class:`~matplotlib.colors.CenteredNorm`).
    plot_kwargs_bias: dict, optional
        Optional keyword arguments for the plot function defined by
        ``plot_func`` for plotting biases. These keyword arguments
        update (and potentially overwrite) the ``plot_kwargs`` for the
        bias plot. By default, uses ``{cmap: 'bwr', norm: 'centered'}``.
        Only relevant for plots including reference datasets.
    projection: str, optional (default: None)
        Projection used for the plot. Needs to be a valid projection
        class of :mod:`cartopy.crs`. Keyword arguments can be specified
        using the option ``projection_kwargs``. For map plots,
        ``'Robinson'`` is used as default.
    projection_kwargs: dict, optional
        Optional keyword arguments for the projection given by
        ``projection``. For map plots, the default keyword arguments
        ``{central_longitude: 10}`` are used.
    pyplot_kwargs: dict, optional
        Optional calls to functions of :mod:`matplotlib.pyplot`.
        Dictionary keys are functions of :mod:`matplotlib.pyplot`.
        Dictionary values are used as argument(s) for these functions
        (if values are dictionaries, these are interpreted as keyword
        arguments; otherwise a single argument is assumed). String
        arguments can include facets in curly brackets which will be
        derived from the corresponding dataset, e.g., ``{project}``,
        ``{short_name}``, ``{exp}``. Examples:
        ``{title: 'Plot {long_name}'}``, ``{xlabel: '{short_name}'}``,
        ``{xlim: [0, 5]}``.
    rasterize: bool, optional (default: False)
        If ``True``, use rasterization_ for plots to produce smaller
        files.  This is only relevant for vector graphics (e.g.,
        ``output_file_type: 'pdf'``).
    show_stats: bool, optional (default: True)
        Show basic statistics on the plots.
    time_format: str, optional (default: None)
        :func:`~datetime.datetime.strftime` format string that is used
        to format the time axis using
        :class:`matplotlib.dates.DateFormatter`. If ``None``, use the
        default formatting imposed by the iris plotting function.
    transpose_axes: bool, optional (default: False)
        Swap X- and Y-axis.
    x_major_formatter: str, optional (default: None)
        Format string for :class:`matplotlib.ticker.FormatStrFormatter`
        used to format major tick labels of X-axis.
    x_minor_formatter: str, optional (default: None)
        Format string for :class:`matplotlib.ticker.FormatStrFormatter`
        used to format minor tick labels of X-axis.
    x_pos_stats_avg: float, optional (default: 0.01)
        Text X-position of average (shown on the left) in Axes
        coordinates. Can be adjusted to avoid overlap with the figure.
        Only relevant if ``show_stats: True``.
    x_pos_stats_bias: float, optional (default: 0.7)
        Text X-position of bias statistics (shown on the right) in Axes
        coordinates. Can be adjusted to avoid overlap with the figure.
        Only relevant if ``show_stats: True``.
    y_major_formatter: str, optional (default: None)
        Format string for :class:`matplotlib.ticker.FormatStrFormatter`
        used to format major tick labels of Y-axis.
    y_minor_formatter: str, optional (default: None)
        Format string for :class:`matplotlib.ticker.FormatStrFormatter`
        used to format minor tick labels of Y-axis.

Recipe configuration options
----------------------------

facet_used_for_labels: str, optional (default: 'dataset')
    Facet used to label different datasets in plot titles and legends.
    For example, ``facet_used_for_labels: 'dataset'`` will use dataset
    names in plot titles and legends; ``facet_used_for_labels: 'exp'``
    will use experiments in plot titles and legends. In addition,
    ``facet_used_for_labels`` is used to select the correct
    ``plot_kwargs`` for the different datasets (see configuration
    options for the different plot types below).
figure_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.figure`. By
    default, uses ``{constrained_layout: True}``.
group_variables_by: str, optional (default: 'short_name')
    Facet or coordinate which is used to create variable groups. For
    each variable group, an individual plot is created. Specifying a
    coordinate allows to create one plot for each point along a
    dimension. For example, when used in combination with the
    preprocessor function :func:`esmvalcore.preprocessor.extract_shape`
    the `shape_id` coordinate can be used to create one plot for each
    shape.
matplotlib_rc_params: dict, optional
    Optional :class:`matplotlib.RcParams` used to customize matplotlib
    plots. Options given here will be passed to
    :func:`matplotlib.rc_context` and used for all plots produced with
    this diagnostic. Note: fontsizes specified here might be overwritten
    by the plot-type-specific option ``fontsize`` (see below).
options: dict, optional
    Additional pre-processing options applied by this diagnostic (see
    list above). Dictonary values are dictonaries used as options for
    the corresponding pre-processing option.
plots: dict
    Plot types plotted by this diagnostic (see list above). Dictionary
    keys must be elements of the list above.  Dictionary values are
    dictionaries used as options for the corresponding plot.
plot_filename: str, optional
    Filename pattern for the plots. By default, uses
    ``'{plot_type}_{real_name}_{dataset}_{mip}_{exp}_{ensemble}'``.
    All tags (i.e., the entries in curly brackets, e.g., ``'{dataset}'``,
    are replaced with the corresponding tags).
plot_folder: str, optional
    Path to the folder to store figures. By default, uses
    ``'{plot_dir}/../../{dataset}/{exp}/{modeling_realm}/{real_name}'``.
    All tags (i.e., the entries in curly brackets, e.g., ``'{dataset}'``,
    are replaced with the corresponding tags). ``'{plot_dir}'`` is
    replaced with the default ESMValTool plot directory (i.e.,
    ``output_dir/plots/diagnostic_name/script_name/``, see
    :ref:`esmvalcore:outputdata`).
savefig_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.savefig`.
    By default, uses ``{bbox_inches: 'tight', dpi: 300, orientation:
    'landscape'}``.
seaborn_settings: dict, optional
    Options for :func:`seaborn.set_theme` (affects all plots). By
    default, uses ``{style: 'ticks'}``.
"""

from __future__ import annotations

import logging
import warnings
from copy import deepcopy
from pathlib import Path
from pprint import pformat
from typing import TYPE_CHECKING, Any

import iris
import iris.analysis
import iris.coord_categorisation as cat
import iris.pandas
import iris.plot
import matplotlib.lines as mlines
import numpy as np
import seaborn as sns
from esmvalcore.preprocessor import area_statistics
from iris.coords import AuxCoord
from iris.cube import Cube
from iris.exceptions import ConstraintMismatchError

import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts.monitor.multi_datasets import MultiDatasets
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    group_metadata,
    io,
    run_diagnostic,
)

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

logger = logging.getLogger(Path(__file__).stem)


class MultiDatasetsThreshold(MultiDatasets):
    """Diagnostic to plot multi-dataset plots."""

    @property
    def options_settings(self) -> dict[str, dict[str, Any]]:
        """pre-plotting settings."""
        default_settings = {
            "threshold": np.nan,
            "inverted": False,
            "accumulated": False,
            "operators": [],
        }
        return {
            "threshold_conversion": {
                "default_settings": {
                    **default_settings,
                },
            },
        }

    def __init__(self, cfg: dict) -> None:
        """Initialize class member."""
        super(MultiDatasets, self).__init__(cfg)

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

        # Check for options/preproc options and initialize them
        if "options" in self.cfg:
            self.options = self.cfg["options"]
        else:
            self.options = {}

        for options_type, option_options in self.options.items():
            if options_type not in self.options_settings:
                msg = (
                    f"Got unexpected options type '{options_type}' for option "
                    f"'options', expected one of {list(self.options_settings)}"
                )
                raise ValueError(msg)
            if option_options is None:
                option_options = {}  # noqa: PLW2901
                self.options[options_type] = option_options

            default_settings_opt = self.options_settings[options_type][
                "default_settings"
            ]
            for key, val in default_settings_opt.items():
                self.options[options_type].setdefault(key, val)

        # Check given plot types and set default settings for them
        for plot_type, plot_options in self.plots.items():
            if plot_type not in self.plot_settings:
                msg = (
                    f"Got unexpected plot type '{plot_type}' for option "
                    f"'plots', expected one of {list(self.plot_settings)}"
                )
                raise ValueError(msg)
            if plot_options is None:
                plot_options = {}  # noqa: PLW2901
                self.plots[plot_type] = plot_options

            # Only use default projection options if no projection is specified
            if "projection" in plot_options:
                self.plots[plot_type].setdefault("projection_kwargs", {})

            default_settings = self.plot_settings[plot_type][
                "default_settings"
            ]

            for key, val in default_settings.items():
                self.plots[plot_type].setdefault(key, val)

        # Load input data
        self.input_data = self._load_and_preprocess_data()
        self.grouped_input_data = group_metadata(
            self.input_data,
            self.cfg["group_variables_by"],
            sort=self.cfg["facet_used_for_labels"],
        )

        # Check that facet_used_for_labels is present for every dataset
        for dataset in self.input_data:
            if self.cfg["facet_used_for_labels"] not in dataset:
                msg = (
                    f"facet_used_for_labels "
                    f"'{self.cfg['facet_used_for_labels']}' not present for "
                    f"the following dataset:\n{pformat(dataset)}"
                )
                raise ValueError(msg)

        # Load seaborn settings
        sns.set_theme(**self.cfg["seaborn_settings"])

    def _check_timeframe(self, cube: Cube) -> bool:
        """Check that the timeframe is a full year period."""
        # Function ``convert_data_thresholded`` not implemented for
        # a start/end date within the year:
        # For partial years the count is biased, since the value
        # corresponds to the scenario where the variable is 0 for the
        # days in the rest of the year. Thus, this leads to an error.

        time = cube.coord("time")
        start = cube.coord("time").bounds[0, 0]
        end = cube.coord("time").bounds[-1, -1]
        startc = time.units.num2date(start)
        endc = time.units.num2date(end)
        start_day = startc.day
        start_month = startc.month
        end_day = endc.day
        end_month = endc.month

        check_timeframe = (
            (start_day == 1 and start_month == 1)
            or (start_day == 31 and start_month == 12)
        ) and (
            (end_day == 31 and end_month == 12)
            or (end_day == 1 and end_month == 1)
        )
        if not check_timeframe:
            msg = (
                "Currently, support for not full year periods is not "
                "implemented cause they could produce missleading results."
            )
            raise ValueError(msg)
        return check_timeframe

    def convert_data_thresholded(
        self,
        cube,
    ):
        """Count the number of days per year on which the threshold is exceeded."""
        # Preventing that this option is executed several times
        if cube.coords("day_of_year"):
            msg = "Reusing already aggregated cube"
            warnings.warn(msg, UserWarning, stacklevel=2)

        else:
            var_unit = cube.units

            cat.add_day_of_year(cube, "time")
            cat.add_year(cube, "time")

            for options_type in self.options:
                # Ensuring that the data is daily, by regridding to daily
                # timestep eventually. Note that for absolute values like
                # temperature one should take the max (accumulated: false),
                # and for cummulated values like total precipitation one
                # should accumulate the values (accumulated: true).
                if self.options[options_type]["accumulated"]:
                    cube = cube.aggregated_by(
                        ["year", "day_of_year"],
                        iris.analysis.SUM,
                    )
                elif self.options[options_type]["inverted"]:
                    cube = cube.aggregated_by(
                        ["year", "day_of_year"],
                        iris.analysis.MIN,
                    )

                else:
                    cube = cube.aggregated_by(
                        ["year", "day_of_year"],
                        iris.analysis.MAX,
                    )

                threshold = self.options[options_type]["threshold"]

                # Count the number of days with values above or below threshold
                if self.options[options_type]["inverted"]:
                    cube = cube.aggregated_by(
                        "year",
                        iris.analysis.COUNT,
                        function=lambda values, threshold=threshold: values
                        < threshold,
                    )
                else:
                    cube = cube.aggregated_by(
                        "year",
                        iris.analysis.COUNT,
                        function=lambda values, threshold=threshold: values
                        > threshold,
                    )

            var_name = cube.var_name
            var_long = cube.long_name

            cube.standard_name = None
            cube.rename(f"{var_name}geq{threshold}count")
            cube.long_name = (
                f"Average number of days per year on which the {var_long} "
                f"exceeds {threshold} {var_unit} at some point"
            )

            cube.units = "days/year"
            for plot_type in self.plots:
                self.plots[plot_type]["cbar_label"] = (
                    f"days with {var_name} exceeding {threshold} {var_unit} [{cube.units}]"
                )
                self.plots[plot_type]["cbar_label_bias"] = (
                    f"Δ days with {var_name} exceeding {threshold} {var_unit} [{cube.units}]"
                )

        return cube

    def _get_netcdf_path(
        self,
        plot_path: Path | str,
        suffix: str | None = None,
        option: str | None = None,
    ) -> str:
        """Get netCDF path."""
        basename = Path(plot_path).stem
        if option is not None:
            basename += "_spacial" + option
        if suffix is not None:
            basename += suffix
        return get_diagnostic_filename(basename, self.cfg)

    def _load_and_preprocess_data(self) -> list[dict]:  # noqa: PLR0912
        """Load and preprocess data."""
        input_data = list(self.cfg["input_data"].values())

        if not input_data:
            msg = "No input data given"
            raise ValueError(msg)

        slices = not any(
            self.cfg["group_variables_by"] in ds for ds in input_data
        )
        datasets = []
        for dataset in input_data:
            filename = dataset["filename"]
            logger.info("Loading %s", filename)
            cubes = iris.load(filename)
            if len(cubes) == 1:
                cube: Cube = cubes[0]
            else:
                var_name = dataset["short_name"]
                try:
                    cube = cubes.extract_cube(
                        iris.NameConstraint(var_name=var_name),
                    )
                except ConstraintMismatchError as exc:
                    var_names = [c.var_name for c in cubes]
                    msg = (
                        f"Cannot load data: multiple variables ({var_names}) "
                        f"are available in file {filename}, but not the "
                        f"requested '{var_name}'"
                    )
                    raise ValueError(msg) from exc

            # Fix time coordinate if present
            if cube.coords("time", dim_coords=True):
                ih.unify_time_coord(cube)

            # Check if the time period are only full years
            self._check_timeframe(cube)

            dims = {"latitude": "degrees_north", "longitude": "degrees_east"}
            for dim, deg in dims.items():
                # Add scalar latitude and longitude coordinates if these are
                # not present (necessary for calculation of area weights). The
                # exact values for the points/bounds of these coordinates do
                # not matter since they don't change the weights.
                if not cube.coords(dim):
                    lon_coord = AuxCoord(
                        0.0,
                        bounds=[-90.0, 90.0],
                        var_name=dim[:3],
                        standard_name=dim,
                        long_name=dim,
                        units=deg,
                    )
                    cube.add_aux_coord(lon_coord, ())

                # Remove additional coordinate systems to avoid errors
                # calculating bias datasets. In particular,  removing the
                # additional coordinate system of the ESACCI obssevational
                # Dataset, introducing a neglegible error.
                if cube.coord(dim).coord_system:
                    msg = (
                        "Removing the coordinate system "
                        f"{cube.coord(dim).coord_system} from the dataset "
                        f"{dataset}  in the dimension {dim}"
                    )
                    warnings.warn(msg, UserWarning, stacklevel=1)
                    cube.coord(dim).coord_system = None

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

            if "threshold_conversion" in self.options:
                cube = self.convert_data_thresholded(cube)

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

        def plot_1d_data(
            cube: Cube,
            label_dataset: str,
            plot_type: str,
            axes: Axes,
            operator: str | None = None,
            linestyle: dict | None = None,
            dataset_colors: dict | None = None,
        ) -> None:
            """Plot single dataset (optional: associated to one operator) in the plot."""
            plot_kwargs = self._get_plot_kwargs(plot_type, dataset)
            coords = self._check_cube_coords(cube, plot_type)
            coord = cube.coord(coords[0], dim_coords=True)

            if operator:
                linestyle_op = linestyle[operator]

                label = f"{label_dataset} - {operator}"
                plot_kwargs.setdefault("linestyle", linestyle_op)
            else:
                label = f"{label_dataset}"

            if dataset_colors:
                plot_kwargs.setdefault("color", dataset_colors[label_dataset])

            plot_kwargs.setdefault("label", label)
            plot_kwargs["axes"] = axes

            if self.plots[plot_type]["transpose_axes"]:
                iris.plot.plot(
                    cube,
                    coord,
                    **plot_kwargs,
                )
            else:
                iris.plot.plot(
                    coord,
                    cube,
                    **plot_kwargs,
                )

        operators = []
        linestyle = {}
        linestyle_iter = iter(["--", "-.", ":", (0, (3, 5, 1, 5, 1, 5))])

        datasets_labels = list(
            dict.fromkeys(self._get_label(d) for d in datasets),
        )

        colors = sns.color_palette("husl", len(datasets_labels))
        dataset_colors = dict(zip(datasets_labels, colors, strict=True))

        multi_dataset_facets = self._get_multi_dataset_facets(datasets)

        if "threshold_conversion" in self.options:
            threshold = self.options["threshold_conversion"]["threshold"]
            operators = self.options["threshold_conversion"]["operators"] or [
                "mean",
            ]
            for operator in operators:
                linestyle[operator] = (
                    "-" if operator == "mean" else next(linestyle_iter, "--")
                )

            axes.set_title(
                "Average number of days per year on which the "
                f"{multi_dataset_facets['long_name']} exceeds "
                f"{threshold} {multi_dataset_facets['units']} at some point",
            )
            var_label = (
                f"{multi_dataset_facets[self.cfg['group_variables_by']]} "
                f"exceeding {threshold} {multi_dataset_facets['units']} "
                "[days/year]"
            )

        else:
            axes.set_title(multi_dataset_facets["long_name"])
            var_label = (
                f"{multi_dataset_facets[self.cfg['group_variables_by']]} "
                f"[{multi_dataset_facets['units']}]"
            )

        for dataset in datasets:
            label_dataset = self._get_label(dataset)
            cube = dataset["cube"]

            # Plotting the observations in black
            for val in dataset.values():
                if "OBS" in str(val):
                    dataset_colors[label_dataset] = "black"

            if "threshold_conversion" in self.options:
                oldcube = cube
                operators = self.options["threshold_conversion"][
                    "operators"
                ] or ["mean"]

                for operator in operators:
                    cube = area_statistics(oldcube, operator=operator)
                    plot_1d_data(
                        cube,
                        label_dataset,
                        plot_type,
                        axes,
                        operator,
                        linestyle,
                        dataset_colors,
                    )

            else:
                operator_list = [
                    cm.method
                    for cm in cube.cell_methods
                    if "latitude" in cm.coord_names
                ]
                if len(operator_list) > 1:
                    msg = (
                        "There are multiple operations accociated with the "
                        "time coordinate, expected is only one. Continuing "
                        "with the first one, but results might be not "
                        "accurate."
                    )
                    warnings.warn(msg, UserWarning, stacklevel=1)
                if len(operator_list) > 0:
                    operator = operator_list[0]

                    if operator not in operators:
                        linestyle[operator] = (
                            "-"
                            if operator == "mean"
                            else next(linestyle_iter, "--")
                        )
                        operators.append(operator)

                    plot_1d_data(
                        cube,
                        label_dataset,
                        plot_type,
                        axes,
                        operator,
                        linestyle,
                        dataset_colors,
                    )

                else:
                    plot_1d_data(
                        cube,
                        label_dataset,
                        plot_type,
                        axes,
                    )

        # Plot horizontal lines
        for hline_kwargs in self.plots[plot_type]["hlines"]:
            axes.axhline(**hline_kwargs)

        # Axis labels
        coords = self._check_cube_coords(cube, plot_type)
        coord = cube.coord(coords[0], dim_coords=True)
        coord_label = f"{coord.name()} [{coord.units}]"

        if self.plots[plot_type]["transpose_axes"]:
            axes.set_xlabel(var_label)
            axes.set_ylabel(coord_label)
        else:
            axes.set_xlabel(coord_label)
            axes.set_ylabel(var_label)

        # Customize plot with user-defined settings
        self._customize_plot(plot_type, axes, multi_dataset_facets)

        # Legend
        col_handles = [
            mlines.Line2D(
                [],
                [],
                color=dataset_colors[dl],
                label=dl,
            )
            for dl in datasets_labels
        ]
        style_handles = [
            mlines.Line2D(
                [],
                [],
                color="gray",
                linestyle=linestyle[o],
                label=o,
            )
            for o in operators
        ]
        handles = col_handles + style_handles
        legend_kwargs = self.plots[plot_type]["legend_kwargs"]
        plt_k = self._get_plot_kwargs(plot_type, dataset).keys()
        if (
            len(style_handles) > 1
            and "linestyle" not in plt_k
            and "color" not in plt_k
        ):
            # Legend for default colors and linestyles
            axes.legend(handles=handles, **legend_kwargs)
        elif legend_kwargs is not False:
            # Legend for custom colors or linestyles
            axes.legend(**legend_kwargs)

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
        fig: Figure = axes.get_figure()

        # Some options are not supported for map plots
        if "map" in plot_type:
            self.plots[plot_type]["aspect_ratio"] = None

        # Plot data
        cube = dataset["cube"]
        if "threshold_conversion" in self.options:
            cube = cube.collapsed("time", iris.analysis.MEAN)
        plot_kwargs = self._get_plot_kwargs(plot_type, dataset, bias=bias)
        plot_kwargs.update(additional_plot_kwargs)
        plot_kwargs["axes"] = axes
        plot_output = self._plot_2d(plot_type, cube, **plot_kwargs)

        # Show coastlines for map plots
        if "map" in plot_type:
            axes.coastlines()

        # Title and axis labels
        fig.suptitle(dataset["long_name"])
        axes.set_title(self._get_label(dataset))
        (x_coord, y_coord) = self._get_coords_for_2d_plotting(plot_type, cube)
        axes.set_xlabel(f"{x_coord.name()} [{x_coord.units}]")
        axes.set_ylabel(f"{y_coord.name()} [{y_coord.units}]")

        # Customize plot with user-defined settings
        self._customize_plot(plot_type, axes, dataset)

        return plot_output

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
        if "threshold_conversion" in self.options:
            operators = self.options["threshold_conversion"]["operators"]
            # for operator in operators:
            cubes_threshold: dict[str, dict[str, Cube]] = {
                operator: {
                    label: area_statistics(cube, operator=operator)
                    for label, cube in cubes.items()
                }
                for operator in operators
            }

        cube_0 = datasets[0]["cube"]
        if "threshold_conversion" in self.options:
            cube_0 = area_statistics(cube_0, operator="mean")
        coord_name = cube_0.coord(dim_coords=True).name()
        var_attrs = {
            n: datasets[0][n] for n in ("short_name", "long_name", "units")
        }

        if "threshold_conversion" in self.options:
            for operator in operators:
                netcdf_path = self._get_netcdf_path(plot_path, option=operator)
                io.save_1d_data(
                    cubes_threshold[operator],
                    netcdf_path,
                    coord_name,
                    var_attrs,
                    attributes={"operator": operator},
                )
        else:
            netcdf_path = self._get_netcdf_path(plot_path)
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


def main(cfg: dict) -> None:
    """Run diagnostic."""
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Using DEFAULT_SPHERICAL_EARTH_RADIUS",
            category=UserWarning,
            module="iris",
        )
        MultiDatasetsThreshold(cfg).compute()


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
