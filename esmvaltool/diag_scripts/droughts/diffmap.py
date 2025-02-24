#!/usr/bin/env  python
# -*- coding: utf-8 -*-
"""Creates a difference map for any given drought index.

A global map is plotted for each dataset with an index (must be unique).
The map shows the difference of the first and last N years
(N = comparison_period).
For multiple datasets a multi-model mean is calculated by default. This can be
disabled using `plot_mmm: False`. To plot only mmm and skip maps for individual
datasets use `plot_models: False`.
The diagnostic is applied to each variable by default, but for single variables
another meta key can be chosen for grouping like `group_by: project` to treat
observations and models seperatly.
The produced maps can be clipped to non polar landmasses (220, 170, -55, 90)
with `clip_land: True`.

TODO: Can MMM be preprocessed in this case? all operations should be linear.
also percent?

TODO: rename metric and group to their real keys in plotkwargs (allow
multi match?) and make group_by accept a list. make sure diffmap_metrics is
always added so that it can be consiedered as extra facet for this diagnostic.
see yml for more notes..

TODO: plot_kwargs overwrites from cfg seems not to overwrite those from
diffmap.yml rename them to 'extra_plot_kwargs' or 'meta_plot_kwargs' and add a
meta key to match any selection (order matters). This allows more flexibility.

Configuration options in recipe
-------------------------------
plot_mmm: bool, optional (default: True)
    Calculate and plot the average over all datasets.
plot_models: bool, optional (default: True)
    Plot maps for each dataset.
basename: str, optional
    Format string for the plot filename. Can use meta keys and diffmap_metric.
    For multi-model mean the dataset will be set to "MMM". Data will be saved
    as same name with .nc extension.
    By default: "{short_name}_{exp}_{diffmap_metric}_{dataset}"
plot_kwargs: dict, optional
    Kwargs passed to diag_scripts.shared.plot.global_contourf function.
    The "cbar_label" parameter is formatted with meta keys. So placeholders
    like "{short_name}" or "{units}" can be used.
    By default {"cmap": "RdYlBu", "extend": "both"}
plot_kwargs_overwrite: list, optional (default: [])
    List of plot_kwargs dicts for specific metrics (diff, first, latest, total)
    and group_by values (ie. pr, tas for group_by: short_name).
    `group` and `metric` can either be strings or lists of strings to be
    applied to all matching plots. Leave any of them empty to apply to all.
    All other given keys are applied to the plot_kwargs dict for this plot.
    Settings will be applied in order of the list, so later entries can
    overwrite previous ones.
comparison_period: int, optional (default: 10)
    Number of years to compare (first and last N years). Must be less or equal
    half of the total time period.
group_by: str, optional (default: short_name)
    Meta key to loop over for multiple datasets.
clip_land: bool, optional (default: False)
    Clips map plots to non polar land area (220, 170, -55, 90).
strip_plots: bool, optional (default: False)
    Removes titles, margins and colorbars from plots (to use them in panels).
mdtol: float, optional (default: 0.5)
    Tolerance for missing data in multi-model mean calculation. 0 means no
    missing data is allowed. For 1 mean is calculated if any data is available.
metrics: list, optional
    List of metrics to calculate and plot. For the difference ("percent" and
    "diff") the mean over two comparison periods ("first" and "last") is
    calculated. The "total" periods mean can be calculated and plotted as well. 
    By default ["first", "last", "diff", "total", "percent"]
"""

import iris
import os
import numpy as np
import yaml
import matplotlib as mpl
from esmvalcore import preprocessor as pp
from iris.analysis import MEAN
import logging
from collections import defaultdict
import esmvaltool.diag_scripts.droughts.utils as ut
import matplotlib.pyplot as plt
import esmvaltool.diag_scripts.shared as e
from cartopy.util import add_cyclic_point

# from esmvaltool.diag_scripts.droughts import colors  # noqa: F401

log = logging.getLogger(__file__)


TITLES = {
    "first": "Mean Historical",
    "last": "Mean Future",
    "trend": "Future - Historical",
    "diff": "Future - Historical",
    "total": "Mean Full Period",
    "percent": "Relative Change",
}

METRICS = ["first", "last", "diff", "total", "percent"]


def plot_colorbar(
    cfg: dict,
    plotfile: str,
    plot_kwargs: dict,
    orientation="vertical",
    mappable=None,
) -> None:
    # fig, ax = plt.subplots(figsize=(1, 4), layout="constrained")
    fig = plt.figure(figsize=(1.5, 3))
    # fixed size axes in fixed size figure
    cbar_ax = fig.add_axes([0.01, 0.04, 0.2, 0.92])
    if mappable is None:
        cmap = plot_kwargs.get("cmap", "RdYlBu")
        norm = mpl.colors.Normalize(
            vmin=plot_kwargs.get("vmin"), vmax=plot_kwargs.get("vmax")
        )
        mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(
        mappable,
        cax=cbar_ax,
        orientation=orientation,
        label=plot_kwargs["cbar_label"],
        pad=0.0,
    )
    if "cbar_ticks" in plot_kwargs:
        cb.set_ticks(plot_kwargs["cbar_ticks"], minor=False)
    fontsize = plot_kwargs.get("cbar_fontsize", 14)
    cb.ax.tick_params(labelsize=fontsize)
    cb.set_label(
        plot_kwargs["cbar_label"], fontsize=fontsize, labelpad=fontsize
    )

    if plotfile.endswith(".png"):
        plotfile = plotfile[:-4]
    fig.savefig(plotfile + "_cb.png")  # , bbox_inches="tight")


def plot(cfg, meta, cube, basename, kwargs=None):
    """Plot map using diag_scripts.shared module."""
    plotfile = e.get_plot_filename(basename, cfg)
    plot_kwargs = cfg.get("plot_kwargs", {}).copy()
    if kwargs is not None:
        plot_kwargs.update(kwargs)
    if "vmax" in plot_kwargs and "vmin" in plot_kwargs:
        plot_kwargs["levels"] = np.linspace(
            plot_kwargs["vmin"], plot_kwargs["vmax"], 9
        )
    label = plot_kwargs.get("cbar_label", "{short_name} ({units})")
    plot_kwargs["cbar_label"] = label.format(**meta)
    for coord in cube.coords(dim_coords=True):
        if not coord.has_bounds():
            log.info("NO BOUNDS GUESSING: %s", coord.name())
            cube.coord(coord.name()).guess_bounds()
    cyclic_data, cyclic_lon = add_cyclic_point(
        cube.data, cube.coord("longitude").points
    )
    if (
        meta["dataset"] == "ERA5"
        and meta["short_name"] == "evspsblpot"
        and len(cube.data[0]) == 360
    ):
        # NOTE: fill missing gap at 360 for era5 pet calculation
        cube.data[:, 359] = cube.data[:, 0]
    mapplot = e.plot.global_contourf(cube, **plot_kwargs)
    if cfg.get("clip_land", False):
        plt.gca().set_extent((220, 170, -55, 90))
        # plt.gcf().set_size_inches(6, 3)
    plt.title(meta.get("title", basename))
    if cfg.get("strip_plots", False):
        plt.gca().set_title(None)
        plt.gca().set_ylabel(None)
        plt.gca().set_xlabel(None)
        cb_mappable = mapplot.colorbar.mappable
        mapplot.colorbar.remove()
        plot_colorbar(cfg, plotfile, plot_kwargs, mappable=cb_mappable)
    fig = mapplot.get_figure()
    fig.savefig(plotfile, bbox_inches="tight")
    plt.close()
    log.info("saved figure: %s", plotfile)


def apply_plot_kwargs_overwrite(kwargs, overwrites, metric, group):
    """Apply plot_kwargs_overwrite to kwargs dict for selected plots."""
    for overwrite in overwrites:
        # print(overwrite)
        new_kwargs = overwrite.copy()
        groups = new_kwargs.pop("group", [])
        if not isinstance(groups, list):
            groups = [groups]
        if len(groups) > 0 and group not in groups:
            continue
        metrics = new_kwargs.pop("metric", [])
        if not isinstance(metrics, list):
            metrics = [metrics]
        if len(metric) > 0 and metric not in metrics:
            continue
        kwargs.update(new_kwargs)
    return kwargs



def calculate_diff(cfg, meta, mm, output_meta, group, norm):
    """absolute difference between first and last years of a cube.
    Calculates the absolut difference between the first and last period of
    a cube. Writing data to mm and plotting each dataset depends on cfg.
    """
    fname = meta["filename"]
    cube = iris.load_cube(fname)
    if meta["short_name"] in cfg.get("convert_units", {}):
        pp.convert_units(cube, cfg["convert_units"][meta["short_name"]])
    try:  # TODO: maybe don't keep this from cmorizer
        cube.remove_coord("Number of stations")  # dropped by unit conversions
    except Exception:
        pass
    if "start_year" in cfg.keys() or "end_year" in cfg.keys():
        log.info("selecting time period")
        # print(cfg.keys())
        cube = pp.extract_time(
            cube, cfg["start_year"], 1, 1, cfg["end_year"], 12, 31
        )
        print(cube.data.shape)
    # if meta["short_name"] in ["evpsblpot"]:
    #     cube.convert_units("1.e-5 kg m-2 s-1")
    dtime = cfg.get("comparison_period", 10) * 12
    cubes = {}
    cubes["total"] = cube.collapsed("time", MEAN)
    do_metrics = cfg.get("metrics", METRICS)
    calc_metrics = ["first", "last", "diff", "percent"]
    if any(m in do_metrics for m in calc_metrics):
        cubes["first"] = cube[0:dtime].collapsed("time", MEAN)
        cubes["last"] = cube[-dtime:].collapsed("time", MEAN)
    if any(m in do_metrics for m in ["diff", "percent"]):
        cubes["diff"] = cubes["last"] - cubes["first"]
        cubes["diff"].data /= norm
        if cubes["diff"].data[0, 0] != np.nan:
            print(cubes["diff"])
        cubes["diff"].units = str(cubes["diff"].units) + " / 10 years"
        cubes["percent"] = cubes["diff"] / cubes["first"] * 100
        cubes["percent"].units = "% / 10 years"
    if cfg.get("plot_mmm", True):
        for key in do_metrics:
            mm[key].append(cubes[key])
    for key, cube in cubes.items():
        if key not in do_metrics:
            continue  # i.e. first/last if only diff is needed
        meta["diffmap_metric"] = key
        meta["exp"] = meta.get("exp", "exp")
        basename = cfg["basename"].format(**meta)
        meta["title"] = f" {basename} ({TITLES[key]})"
        if cfg.get("plot_models", True):
            plot_kwargs = cfg.get("plot_kwargs", {}).copy()
            overwrites = cfg.get("plot_kwargs_overwrite", [])
            apply_plot_kwargs_overwrite(plot_kwargs, overwrites, key, group)
            plot(cfg, meta, cube, basename, kwargs=plot_kwargs)
            plt.close()
        if cfg.get("save_models", True):
            work_file = os.path.join(cfg["work_dir"], f"{basename}.nc")
            iris.save(cube, work_file)
            meta["filename"] = work_file
            output_meta[work_file] = meta.copy()


def calculate_mmm(cfg, meta, mm, output_meta, group, key="diff"):
    """Calculate multi-model mean for a given metric."""
    drop = cfg.get("dropcoords", ["time", "height"])
    meta = meta.copy()  # don't modify meta in place:
    meta["dataset"] = "MMM"
    meta["diffmap_metric"] = key
    basename = cfg["basename"].format(**meta)
    mmm, _ = ut.mmm(
        mm[key],
        dropcoords=drop,
        dropmethods=key != "diff",
        mdtol=cfg.get("mdtol", 0.3),
        # mdtol=0,
    )
    meta["title"] = f"Multi-model Mean ({cfg['titles'][key]})"
    if cfg.get("plot_mmm", True):
        plot_kwargs = cfg.get("plot_kwargs", {}).copy()
        overwrites = cfg.get("plot_kwargs_overwrite", [])
        apply_plot_kwargs_overwrite(plot_kwargs, overwrites, key, group)
        plot(cfg, meta, mmm, basename, kwargs=plot_kwargs)
    if cfg.get("save_mmm", True):
        work_file = os.path.join(cfg["work_dir"], f"{basename}.nc")
        meta["filename"] = work_file
        meta["diffmap_metric"] = key
        output_meta[work_file] = meta.copy()
        iris.save(mmm, work_file)


def set_defaults(cfg):
    """update cfg with default values from diffmap.yml"""
    config_file = os.path.realpath(__file__)[:-3] + ".yml"
    with open(config_file, "r", encoding="utf-8") as f:
        defaults = yaml.safe_load(f)
    for key, val in defaults.items():
        cfg.setdefault(key, val)
    if cfg["plot_kwargs_overwrite"] is not defaults["plot_kwargs_overwrite"]:
        cfg["plot_kwargs_overwrite"].extend(defaults["plot_kwargs_overwrite"])


def main(cfg):
    """Main function."""
    # cfg["group_by"] = cfg.get("group_by", "short_name")
    set_defaults(cfg)
    groups = e.group_metadata(cfg["input_data"].values(), cfg["group_by"])
    output = {}
    for group, metas in groups.items():
        mm = defaultdict(list)
        skipped = 0
        for meta in metas:
            # TODO: fix diag_spei output to contain all relevant meta data
            ut.guess_experiment(meta)
            if "end_year" not in meta:
                try:
                    meta.update(ut.get_time_range(meta["filename"]))
                except Exception:
                    log.error(
                        "failed to get time range for %s", meta["filename"]
                    )
                    skipped += 1
                    log.error("skipped datasets: %s", skipped)
                    continue
            # adjust norm for selected time period
            meta["end_year"] = cfg.get("end_year", meta["end_year"])
            meta["start_year"] = cfg.get("start_year", meta["start_year"])
            norm = (
                int(meta["end_year"])
                - int(meta["start_year"])
                + 1  # count full end year
                - cfg.get("comparison_period", 10)  # decades center to center
            ) / 10
            calculate_diff(cfg, meta, mm, output, group, norm)
        do_mmm = cfg.get("plot_mmm", True) or cfg.get("save_mmm", True)
        if do_mmm and len(metas) > 1:
            for metric in cfg.get("metrics", METRICS):
                calculate_mmm(cfg, metas[0], mm, output, group, metric)
    ut.save_metadata(cfg, output)
    # TODO close all and everything to free up memory
    # if "panels" in cfg:
    #     for grid in cfg["panels"]:
    #         create_panels(cfg, output, grid)


if __name__ == "__main__":
    with e.run_diagnostic() as config:
        main(config)
