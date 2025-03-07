""" Plot timeseries of historical period and different ssps for each variable.
Observations will be shown as reference when given. Shaded area shows MM stddv.
Works with combined or individual inputs for historical/ssp experiments.

Configuration options in recipe
-------------------------------
plot_mmm: bool, optional (default: True)
smooth: bool, optional (default: True)
    Yearly averages for each variable, before mm operations and plot.
combined_split_years: int, optional (default: 65)
    If experiments are already combined, this is the number of full years that
    is used to split the data into two seperated experiments. Historical data 
    is only plotted once.
plot_properties: dict, optional
    Additional properties to set on the plot. Passed to ax.set().
figsize: tuple, optional (default: (9, 2))
reuse_mm: bool, optional (default: False)
subplots: bool, optional (default: False)
    Plot all time series as subplots in one figure with shared x-axis.
legend: dict, optional (default: {})
    Names to rename the default legend labels. Keys are the original labels.
"""

from cProfile import label
import logging
from copy import deepcopy
from pathlib import Path
from esmvalcore import preprocessor as pp
from esmvaltool.diag_scripts.droughtindex import utils as ut
from esmvaltool.diag_scripts.droughtindex import styles
import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.axes import Axes
import iris
from iris import plot as iplt
from iris.iterate import izip
from iris.analysis import MEAN, STD_DEV, cartography
from iris.plot import _fixup_dates
from iris.coord_categorisation import add_year
import datetime as dt
from esmvaltool.diag_scripts.shared import (
    # ProvenanceLogger,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
    select_metadata,
)

# import nc_time_axis  # noqa  allow cftime axis to be plotted by mpl
# nc_time_axis works but seems to show wrong days on axis when using cftime
logger = logging.getLogger(Path(__file__).stem)


def convert_units(cube):
    """Convert units of some variables for display"""
    if cube.var_name in ["tas", "tasmax", "tasmin"]:
        cube.convert_units("Celsius")
    if cube.var_name == "pr":
        logger.info("Converting pr units to mm/day")
        cube.units = "mm s-1"
        cube.convert_units("mm day-1")
    if cube.var_name == "evspsblpot":
        # cube.units = 'mm mon-1'
        ut.monthly2daily(cube)
        cube.long_name = "Potential Evapotranspiration"  # NOTE: not working?
        cube.rename("Potential Evapotranspiration")
    return None


def global_mean(cfg, cube):
    """Calculate global mean."""
    ut.guess_lat_lon_bounds(cube)
    if "regions" in cfg:
        print("Extracting regions")
        cube = pp.extract_shape(cube, shapefile='ar6', ids={"Acronym": cfg["regions"]})
    area_weights = cartography.area_weights(cube)
    mean = cube.collapsed(
        ["latitude", "longitude"], MEAN, weights=area_weights
    )
    return mean


def yearly_average(cube):
    add_year(cube, "time")
    return cube.aggregated_by("year", MEAN)


def plot_experiment(cfg, mean, std_dev, experiment, ax):
    time = mean.coord("time")
    # times = time.units.num2date(time.points)
    exp_color = getattr(styles, experiment)
    iplt.fill_between(
        time,
        mean - std_dev,
        mean + std_dev,
        color=exp_color,
        alpha=0.2,
        axes=ax,
    )
    iplt.plot(time, mean, color=exp_color, label=experiment, axes=ax)
    y_label = f"{mean.var_name} [{mean.units}]"
    if mean.long_name:
        y_label = f"{mean.long_name}\n[{mean.units}]"
    y_label = cfg.get("ylabels", {}).get(mean.var_name, y_label)
    y_label = cfg.get("ylabels", {}).get(mean.long_name, y_label)
    ax.set_ylabel(y_label)


def plot_each_model(cubes, metas, cfg, experiment, smooth=False):
    fig, ax = plt.subplots(figsize=cfg["figsize"], dpi=150)
    ax.grid(axis="y", color="0.95")
    time = cubes[0].coord("time")
    for cube, meta in zip(cubes, metas):
        iplt.plot(time, cube, label=meta["dataset"])
    if not cfg.get("strip_plots", False):
        ax.legend()
    basename = f"timeseries_scenarios_{meta['short_name']}_{experiment}"
    fig.savefig(get_plot_filename(basename, cfg), bbox_inches="tight")


def plot_models(cfg, metas, ax, smooth=False):
    historical_plotted = False
    for experiment, models in group_metadata(metas, "exp").items():
        if experiment == "historical" and historical_plotted:
            continue

        basename = f"{experiment}_{metas[0]['short_name']}"
        fname = os.path.join(cfg["work_dir"], f"{basename}")
        recalc = not cfg.get("reuse_mm", False)
        if cfg.get("reuse_mm", False):
            try:
                std_dev = iris.load_cube(fname + "_stddev.nc")
                mean = iris.load_cube(fname + "_mean.nc")
            except FileNotFoundError:
                recalc = True
        if recalc or cfg.get("plot_models", False):
            cubes = [
                global_mean(cfg, iris.load_cube(meta["filename"]))
                for meta in models
            ]
            if smooth:
                cubes = [yearly_average(cube) for cube in cubes]
        if cfg.get("plot_models", False):
            # TODO: convert units for single models?
            plot_each_model(cubes, models, cfg, experiment, smooth=smooth)
        # mean = global_mean(mm["mean"])
        if recalc:
            mm = pp.multi_model_statistics(
                cubes, "overlap", ["mean", "std_dev"]
            )
            std_dev = mm["std_dev"]  # global_mean(mm["std_dev"])
            mean = mm["mean"]
        if recalc and cfg.get("save_mm", True):
            iris.save(mm["mean"], fname + "_mean.nc")
            iris.save(mm["std_dev"], fname + "_stddev.nc")
        convert_units(mean)
        convert_units(std_dev)

        if experiment.startswith("historical-"):
            experiment = experiment.split("-")[1]
            steps = 65 if smooth else 65 * 12
            plot_experiment(
                cfg,
                mean[(steps - 1) :],
                std_dev[(steps - 1) :],
                experiment,
                ax,
            )
            if not historical_plotted:
                plot_experiment(
                    cfg, mean[:steps], std_dev[:steps], "historical", ax
                )
                historical_plotted = True
            continue
        if experiment == "historical":
            historical_plotted = True
        plot_experiment(cfg, mean, std_dev, experiment, ax)


def plot_obs(cfg, metas, ax, smooth=False):
    for meta in metas:
        cube = iris.load_cube(meta["filename"])
        if smooth:
            cube = yearly_average(cube)
        mean = global_mean(cfg, cube)
        convert_units(mean)
        time = mean.coord("time")
        iplt.plot(time, mean, linestyle="--", label=meta["dataset"], axes=ax)


def process_variable(cfg, metas, short_name, fig=None, ax: Axes = None):
    """Process variable."""
    project = cfg.get("project", "CMIP6")
    model_metas = select_metadata(metas, project=project)
    obs_metas = [meta for meta in metas if meta["project"] != project]
    if not cfg.get("subplots", False):
        fig, ax = plt.subplots(figsize=cfg.get("figsize", (9, 2)), dpi=300)
    plot_models(cfg, model_metas, ax, smooth=cfg.get("smooth", False))
    plot_obs(cfg, obs_metas, ax, smooth=cfg.get("smooth", False))
    basename = f"timeseries_scenarios_{short_name}"
    if not cfg.get("subplots", False):
        ax.set_xlabel("Time")
        ax.legend()
    else:
        ax.set_xticklabels([])
        ax.set_xticks([])
    ax.grid(
        axis="both", color="0.6", which="major", linestyle="--", linewidth=0.5
    )
    # ax.set_frame_on(False)
    ax.set_xlim([dt.datetime(1950, 1, 1), dt.datetime(2100, 1, 1)])
    ax.xaxis.set_major_locator(mdates.YearLocator(10))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    for label in ax.get_xticklabels(which="major"):
        label.set(rotation=40, horizontalalignment="right")
    ax.xaxis.set_minor_locator(mdates.YearLocator())
    if "plot_properties" in cfg.keys():
        ax.set(**cfg["plot_properties"])
    if not cfg.get("subplots", False):
        fig.savefig(get_plot_filename(basename, cfg), bbox_inches="tight")


def main(cfg):
    """Run diagnostic."""
    cfg = deepcopy(cfg)
    groups = group_metadata(cfg["input_data"].values(), "short_name")
    fig, axs = None, [None] * len(groups)
    if cfg["subplots"]:
        figsize = cfg.get("figsize", (9, 2))
        height = len(groups) * (figsize[1] + 0.5)
        fig, axs = plt.subplots(
            len(groups), 1, figsize=(figsize[0], height), dpi=300
        )
    for i, (short_name, metas) in enumerate(groups.items()):
        process_variable(cfg, metas, short_name, fig=fig, ax=axs[i])
    if cfg["subplots"]:
        basename = "timeseries_scenarios"
        for ax in axs:
            # ax.set_xticklabels([])
            ax.tick_params(
                axis="x",
                which="both",
                bottom=False,
                top=False,
                labelbottom=False,
            )
            ax.tick_params(
                axis="y",
                which="both",
                left=True,
                right=True,
                labelleft=False,
                labelright=True,
            )
            for spine in ["top", "bottom"]:
                ax.spines[spine].set_visible(False)
        axs[-1].tick_params(
            axis="x", which="both", bottom=True, top=False, labelbottom=True
        )
        axs[-1].spines["bottom"].set_visible(True)
        axs[0].spines["top"].set_visible(True)
        lines, labels = axs[-1].get_legend_handles_labels()
        if cfg.get(
            "legend", cfg.get("subplots", False)
        ):  # rename and reorder handles and labels
            leg_dict = dict(zip(labels, lines))
            print(labels)
            labels = list(cfg["legend"].values())
            handles = [leg_dict[lab] for lab in cfg["legend"].keys()]
            axs[-1].legend(handles, labels)
        fig.subplots_adjust(hspace=0.02)
        fig.tight_layout()
        fig.savefig(get_plot_filename(basename, cfg), bbox_inches="tight")


def set_defaults(cfg):
    cfg.setdefault("plot_mmm", True)
    cfg.setdefault("smooth", True)
    cfg.setdefault("combined_split_years", 65)
    cfg.setdefault("plot_properties", {})
    cfg.setdefault("figsize", (9, 2))
    cfg.setdefault("reuse_mm", False)
    cfg.setdefault("subplots", False)
    cfg.setdefault("legend", {})
    cfg.setdefault("ylabels", {})


if __name__ == "__main__":
    with run_diagnostic() as cfg:
        set_defaults(cfg)
        main(cfg)
