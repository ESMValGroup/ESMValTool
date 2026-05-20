"""
Plot types can be specified with the recipe option 'plots'; Pre-processing options can be accessed with the recipe option 'options' 

"""


from __future__ import annotations

import inspect
import logging
import warnings
from copy import copy, deepcopy
from functools import partial
from pathlib import Path
from pprint import pformat
from typing import TYPE_CHECKING, Any

import cartopy.crs as ccrs
import dask.array as da
import iris
import iris.analysis
import iris.pandas
import iris.plot
import iris.coord_categorisation as cat
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
from matplotlib.colors import CenteredNorm
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

#TODO: make this import from esmvalcore unnecessary:
from esmvalcore.preprocessor._shared import (
    apply_mask,
    get_dims_along_axes,
    get_iris_aggregator,
    get_normalized_cube,
    preserve_float_dtype,
    try_adding_calculated_cell_area,
    update_weights_kwargs,
)

from esmvalcore.iris_helpers import (
    has_regular_grid,
    ignore_iris_vague_metadata_warnings,
)


if TYPE_CHECKING:
    from collections.abc import Callable, Iterable

    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

logger = logging.getLogger(Path(__file__).stem)

class MultiDatasets(MonitorBase):
    """Diagnostic to plot multi-dataset plots."""

    @property
    def options_settings(self) -> dict[str, dict[str, Any]]:
        """pre-plotting settings."""
        default_settings = {
                "threshold": np.nan,
                "inverted": False,
                "accumulated": False,  # This should only be swiched to true for variables which are accumulating over time.
                "operators": [],
        }
        return{
            "threshold_conversion": {
                "function": partial(self.convert_data_thresholded, "threshold_conversion", ),
                "provenance": {
                    "authors": ["caspers_laura"],
                    "caption": "Converting daily data to count of days per year on which thresholds are exceeded.",
                },
                "default_settings": {
                    **default_settings,
                    
                },
            },
            

        }

    
    @property
    def plot_settings(self) -> dict[str, dict[str, Any]]:
        """Plot settings."""
        default_settings_1d = {
            "aspect_ratio": None,
            "axes_kwargs": {},
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
            "axes_kwargs": {},
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
        return {
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

        }

    def __init__(self, cfg: dict) -> None:
        """Initialize class member."""
        super().__init__(cfg)

        ###TODO maybe edit here for labels...
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

       
       


        #load input data
        self.input_data = self._load_and_preprocess_data()
        self.grouped_input_data = group_metadata(
            self.input_data,
            self.cfg["group_variables_by"],
            sort=self.cfg["facet_used_for_labels"],
        )

        #check for options/preproc options and initialize them
        if "options" in self.cfg:
            self.options = self.cfg["options"]
        else:
            self.options = {}


        for options_type, option_options in self.options.items():
            if options_type not in self.options_settings:
                msg = (
                    f"Got unexpected options type '{options_type}' for option "
                    f"'plots', expected one of {list(self.options_settings)}"
                )
                raise ValueError(msg) 
            if option_options is None:
                option_options = {}  # noqa: PLW2901
                self.options[options_type] = option_options
        
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

            default_settings_opt = self.options_settings[options_type][
                "default_settings"
            ]
            for key, val in default_settings_opt.items():
                self.options[options_type].setdefault(key, val)

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

    def _add_colorbar(  # noqa: PLR0913
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
            msg = f"plot_type '{plot_type}' not supported"
            raise NotImplementedError(msg)

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
        if np.abs(mean.data) >= 0.1:  # noqa: PLR2004
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
        if np.abs(rmse.data) >= 0.1:  # noqa: PLR2004
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
        msg = (
            f"Expected cube with dimensional coordinates "
            f"{expected_dims_str}, got {cube.summary(shorten=True)}"
        )
        raise ValueError(msg)

    def _customize_plot(  # noqa: PLR0912
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
                subs=np.arange(1.0, 10.0) * 0.1,
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
        self._process_axes_kwargs(
            axes,
            self.plots[plot_type]["axes_kwargs"],
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
            msg = (
                f"Not all necessary facets in {description} available for "
                f"dataset\n{pformat(dataset)}"
            )
            raise ValueError(msg) from exc
        return string

    def _get_benchmark_datasets(self, datasets: list[dict]) -> list[dict]:
        """Get dataset to be benchmarked."""
        variable = datasets[0][self.cfg["group_variables_by"]]
        benchmark_datasets = [
            d for d in datasets if d.get("benchmark_dataset", False)
        ]
        if len(benchmark_datasets) >= 1:
            return benchmark_datasets

        msg = (
            f"Expected at least 1 benchmark dataset (with 'benchmark_dataset: "
            f"True' for variable '{variable}'), got "
            f"{len(benchmark_datasets):d}"
        )
        raise ValueError(msg)

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
            msg = (
                f"Could not create benchmarking mask, unknown benchmarking "
                f"metric: '{metric}'"
            )
            raise ValueError(msg)

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
            msg = f"Unknown benchmarking metric: '{metric}'."
            raise ValueError(msg)

        if len(percentile_datasets) >= n_percentiles[metric]:
            return percentile_datasets

        variable = datasets[0][self.cfg["group_variables_by"]]
        msg = (
            f"Expected at least {n_percentiles[metric]} percentile datasets "
            f"(created with multi-model statistics preprocessor for variable "
            f"'{variable}'), got {len(percentile_datasets):d}"
        )
        raise ValueError(msg)

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
        return self._fill_facet_placeholders(cbar_label, dataset, descr)

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
            msg = (
                f"Got invalid plot function '{plot_func}' for plotting "
                f"{plot_type}, expected function of iris.plot"
            )
            raise AttributeError(msg)
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
                val = self._fill_facet_placeholders(  # noqa: PLW2901
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
            msg = (
                f"Got invalid projection '{projection}' for plotting "
                f"{plot_type}, expected class of cartopy.crs"
            )
            raise AttributeError(msg)

        return getattr(ccrs, projection)(**projection_kwargs)

#TODO: Update the next function here:
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
            msg = (
                f"Expected at most 1 reference dataset (with "
                f"'reference_for_monitor_diags: True' for variable "
                f"'{variable}', got {len(ref_datasets):d}"
            )
            raise ValueError(msg)
        if ref_datasets:
            return ref_datasets[0]
        return None

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
            if self.options["threshold_conversion"]:
                oldcube = cube
                cube = self.thr_area_statistics(cube, operator = "mean")
                for operator in self.options["threshold_conversion"]["operators"]:
                    cube = self.thr_area_statistics(oldcube, operator = operator) #TODO: append instead of overwrite...
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
            cube = cube.copy(npx.ma.filled(cube.core_data(), np.nan))

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
        fig: Figure = axes.get_figure()

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


























            
    def convert_data_thresholded(
        self, 
        options_type: str,
        datasets: list[dict],
        collapsed = False,
    ) -> list[dict]:
        for dataset in datasets:
            cube = dataset["cube"]
            var_unit = cube.units
            #TODO: this condition is a bit hacky fix to prevent for this option to be executed several times, would be better to fix it nicely, directly calling this option only once before the plotting calls.
            # if len(cube.coords('day_of_year')) > 0:
            #     msg = (
            #         f"WARNING: Reusing already aggregated cube..."
            #     )
            #     warnings.warn(msg, ESMValToolDeprecationWarning, stacklevel=2)

            # else:
            cat.add_day_of_year(cube, "time")
            cat.add_year(cube, "time")

            # Ensuring that the data is daily, by regridding to daily timestep eventually. 
            # Note that for absolute values like temperature one should take the max (accumulated: false), and for cummulated values like total precipitation one should accumulate the values (accumulated: true).  
            if self.options[options_type]["accumulated"]:
                cube = cube.aggregated_by(["year", "day_of_year"], iris.analysis.SUM)
            else:
                cube = cube.aggregated_by(["year", "day_of_year"], iris.analysis.MAX)

            threshold = self.options[options_type]["threshold"]
        
            # Count the number of days with values above or below threshold
            if self.options[options_type]["inverted"]:
                cube = cube.aggregated_by("year", iris.analysis.COUNT, function = lambda values: values < threshold)
            else:
                cube = cube.aggregated_by("year", iris.analysis.COUNT, function = lambda values: values > threshold)
        
            var_name = cube.var_name
            var_long = cube.long_name
            #var_comment = cube.comment

            # ###todo, add this for maps:
            # #count_cube = oldcube.collapsed("time", iris.analysis.COUNT, function = lambda values: values > threshold)
            # #count_cube.data = (count_cube.data / time_length) * 360
            # # count_cube.long_name = f"count_of_days_per_year_with_{var_name}_above_{threshold}_ADDUNITHERE"
            #cube.var_name = f"{var_name}geq{threshold}count"
            cube.standard_name = None
            cube.rename(f"{var_name}geq{threshold}count")
            cube.long_name = f"Average number of days per year on which the near-surface (usually, 2 meter) air temperature exceeds {threshold} {var_unit} at some point" #TODO: make auto insertions, also for °C...
            # cube.comment = f"TODO..."
            cube.units = "days/year"
        
            #if len(cubes) == 1:
            #    cube: Cube = count_cube
            # cube: Cube = cubes[0]
            dataset["cube"] = cube
            dataset["standard_name"] = None
            dataset["var_name"] =  f"{var_name}geq{threshold}count"
            dataset["long_name"] = f"Number of days per year on which the {var_long} exceeds the threshold of {threshold} {var_unit}" #TODO: make auto insertions, also for °C...
            #dataset["comment"] = f"Average number of days per year on which the {var_comment} exceeds the threshold of {threshold} {var_unit} at some point"
            dataset["units"] = "days/year"
        
        ######################################################################   
           
            #TODO: put this in the map etc routine, such that the original cube remains intact for other kinds of plots...:
            #for maps etc. use the mean of the yearly datapoints:
            if collapsed:
                cube = cube.collapsed("time", iris.analysis.MEAN)
                dataset["cube"] = cube
 

          
           

            return datasets




    def _process_axes_kwargs(
        self,
        axes: Axes,
        axes_kwargs: dict[str, Any],
        dataset: dict,
    ) -> None:
        """Process functions for :class:`matplotlib.axes.Axes`."""
        for func, arg in axes_kwargs.items():
            # For set_extent, make sure to specify coordinate system that is
            # used (will always be PlateCarree; see
            # https://stackoverflow.com/questions/43470238/cartopy-set-extent-extending-requested-boundary#43505490)
            if func == "set_extent":
                if not isinstance(arg, dict):
                    arg = {"extents": arg}  # noqa: PLW2901
                arg["crs"] = ccrs.PlateCarree()
            if isinstance(arg, str):
                arg = self._fill_facet_placeholders(  # noqa: PLW2901
                    arg,
                    dataset,
                    f"axes_kwargs '{func}: {arg}'",
                )
            if arg is None:
                getattr(axes, func)()
            elif isinstance(arg, dict):
                getattr(axes, func)(**arg)
            else:
                getattr(axes, func)(arg)

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
                func = func.replace("x", "y", 1)  # noqa: PLW2901
            elif transpose_axes and func.startswith("y"):
                func = func.replace("y", "x", 1)  # noqa: PLW2901
            if isinstance(arg, str):
                arg = self._fill_facet_placeholders(  # noqa: PLW2901
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
        if self.options["threshold_conversion"]:
            for label, cube in cubes.items():
                cubes[label] = self.thr_area_statistics(cube, operator = "mean")
        netcdf_path = self._get_netcdf_path(plot_path)
        #coord_name = datasets[0]["cube"].coord(dim_coords=True).name()
        cube_0 = datasets[0]["cube"]
        if self.options["threshold_conversion"]:
            cube_0 = self.thr_area_statistics(cube_0, operator = "mean")
        coord_name = cube_0.coord(dim_coords=True).name()
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


# @register_supplementaries(
#     variables=["areacella", "areacello"],
#     required="prefer_at_least_one",
# )
# Adapted preprocessor function "area_statistics".
    #@preserve_float_dtype #-> function in preprocessor/_shared.py; TODO: can we just do without?
    def thr_area_statistics(
        self,
        cube: Cube,
        operator: str,
        normalize: Literal["subtract", "divide"] | None = None,
        **operator_kwargs: Any,
    ) -> list[dict]:
        # for dataset in datasets:
        #     cube = dataset["cube"]

        has_cell_measure = bool(cube.cell_measures("cell_area"))

        # Get aggregator and correct kwargs (incl. weights)
        (agg, agg_kwargs) = get_iris_aggregator(operator, **operator_kwargs)
        agg_kwargs = update_weights_kwargs(
            operator,
            agg,
            agg_kwargs,
            "cell_area",
            cube,
            try_adding_calculated_cell_area,
        )

        with ignore_iris_vague_metadata_warnings():
            result = cube.collapsed(["latitude", "longitude"], agg, **agg_kwargs)
        if normalize is not None:
            result = get_normalized_cube(cube, result, normalize)

        # Make sure input cube has not been modified
        if not has_cell_measure and cube.cell_measures("cell_area"):
            cube.remove_cell_measure("cell_area")

        # dataset["cube"] = result
        return result






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

        """Adjust data."""
        for options_type in self.options:
            options_settings = self.options_settings[options_type]
            options_function = options_settings["function"]
            logger.info("Executing option %s", options_type)
            
            benchmark_datasets = options_function(datasets)
        
        # For threshold data, do the spatial mean here:
        if self.options["threshold_conversion"]:
            benchmark_datasets = self.thr_area_statistics(datasets = benchmark_datasets, operator = "mean")
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
        """Adjust data."""
        for options_type in self.options:
            options_settings = self.options_settings[options_type]
            options_function = options_settings["function"]
            logger.info("Executing option %s", options_type)
            
            datasets = options_function(datasets)

        # # For threshold data, do the spatial mean here:
        # if self.options["threshold_conversion"]:
        #     datasets_1d = self.thr_area_statistics(datasets = datasets, operator = "mean")
        #     self._plot_1d_data(plot_type, datasets_1d, axes)
        #     self._save_1d_data(plot_type, datasets_1d, fig)
        # else:
        self._plot_1d_data(plot_type, datasets, axes)
        self._save_1d_data(plot_type, datasets, fig)

    def create_2d_benchmarking_plot(
        self,
        plot_type: str,
        datasets: list[dict],
    ) -> None:
        for options_type in self.options:
            options_settings = self.options_settings[options_type]
            options_function = options_settings["function"]
            logger.info("Executing option %s", options_type)
            
            datasets = options_function(datasets, collapsed = True)
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
        for options_type in self.options:
            options_settings = self.options_settings[options_type]
            options_function = options_settings["function"]
            logger.info("Executing option %s", options_type)
            
            datasets = options_function(datasets, collapsed = True)
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




    def compute(self) -> None:
        """Plot preprocessed data."""
        for plot_type in self.plots:
            plot_settings = self.plot_settings[plot_type]
            plot_function = plot_settings["function"]
            mpl_rc_params = self._get_custom_mpl_rc_params(plot_type)
            logger.info("Plotting %s", plot_type)


                   
            # Handle deprecations -> see manuels orriginal diagnostic for more on this, deleted here cause it is not needed
            
            # Inspect plot function to determine arguments
            plot_parameters = inspect.signature(plot_function).parameters



            # Plot types where only one plot in total is created
            if not plot_parameters:
                with mpl.rc_context(mpl_rc_params):


    # """Adjust data."""
  #                  for options_type in self.options:
   #                     options_settings = self.options_settings[options_type]
    #                    options_function = options_settings["function"]
     #                   logger.info("Executing option %s", options_type)
      #      
       #                 adj_datas = options_function(datasets)
        #                plot_function(adj_datas)
         #           else
                        plot_function()





            # Plot types where multiple plots might be created
            else:
           #     msg = (
            #        f"WARNING: This is not implemented yet..."
            #    )
            #    warnings.warn(msg, ESMValToolDeprecationWarning, stacklevel=2)
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
        









