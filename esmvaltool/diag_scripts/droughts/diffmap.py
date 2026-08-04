"""Plot relative and absolute differences between two time intervals.

A global map is plotted for each dataset with an index (must be unique).
The map shows the difference of the first and last N years
(N = comparison_period).
For multiple datasets a multi-model mean is calculated by default. This can be
disabled using `plot_mmm: False`. To plot only mmm and skip maps for individual
datasets use `plot_models: False`.
The diagnostic is applied to each variable by default, but for single variables
another meta key can be chosen for grouping like `group_by: project` to treat
observations and models separately.
The produced maps can be clipped to non polar landmasses (220, 170, -55, 90)
with `clip_land: True`.


Configuration options in recipe
-------------------------------
basename: str, optional
    Format string for the plot filename. Can use meta keys and diffmap_metric.
    For multi-model mean the dataset will be set to "MMM". Data will be saved
    as same name with .nc extension.
    By default: "{short_name}_{exp}_{diffmap_metric}_{dataset}"
clip_land: bool, optional (default: False)
    Clips map plots to non polar land area (220, 170, -55, 90).
comparison_period: int, optional (default: 10)
    Number of years to compare (first and last N years). Must be less or equal
    to half of the total time period.
filters: dict or list, optional
    Filter for metadata keys to select datasets. Only datasets with matching
    values will be processed. This can be useful, if ancestors or preprocessed
    data is available, that should not be processed by the diagnostic.
    If a list of dicts is given, all datasets matching any of the filters will
    be considered.
    By default None.
group_by: str, optional (default: short_name)
    Meta key to loop over for multiple datasets.
metrics: list, optional
    List of metrics to calculate and plot. For the difference ("percent" and
    "diff") the mean over two comparison periods ("first" and "last") is
    calculated. The "total" periods mean can be calculated and plotted as well.
    By default ["first", "last", "diff", "total", "percent"]
mdtol: float, optional (default: 0.5)
    Tolerance for missing data in multi-model mean calculation. 0 means no
    missing data is allowed. At 1 the mean is calculated if any data is
    available.
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
plot_mmm: bool, optional (default: True)
    Calculate and plot the average over all datasets.
plot_models: bool, optional (default: True)
    Plot maps for each dataset.
strip_plots: bool, optional (default: False)
    Removes titles, margins and colorbars from plots (to use them in panels).
titles: dict, optional
    Customize plot titles for different metrics. Possible dict keys are
    "first", "last", "trend", "diff", "total", "percent". The values are
    formatted using meta data. Placeholders like "{short_name}" can be used.
    By default {"first": "Mean Historical", "last": "Mean Future",
    "trend": "Future - Historical", "diff": "Future - Historical",
    "total": "Mean Full Period", "percent": "Relative Change"}.
"""

from __future__ import annotations

import logging
import os
from collections import defaultdict
from pathlib import Path

import iris
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import yaml
from cartopy.util import add_cyclic_point
from esmvalcore import preprocessor as pp
from iris.analysis import MEAN
from iris.cube import Cube

from esmvaltool.diag_scripts.droughts import utils
from esmvaltool.diag_scripts.shared import (
    get_plot_filename,
    group_metadata,
    run_diagnostic,
    select_metadata,
)
from esmvaltool.diag_scripts.shared.plot import global_contourf

log = logging.getLogger(__file__)


METRICS = ["first", "last", "diff", "total", "percent"]
PROVENANCE = {
    "authors": ["lindenlaub_lukas"],
    "domains": ["global"],
    "plot_types": ["map"],
}


def _get_provenance(cfg: dict, meta: dict) -> dict:
    """Create provenance dict for single model plots."""
    prov = PROVENANCE.copy()
    prov["statistics"] = ["mean"]
    dataset = meta.get("dataset", "unknown")
    if dataset == "MMM":
        dataset = "Multi-Model Mean"
    if meta["diffmap_metric"] == "diff":
        prov["statistics"] = ["diff", "mean"]
        prov["caption"] = (
            f"Absolute difference in {meta['long_name']} between first and "
            f"last {cfg['comparison_period']} years of the period "
            f"{meta['start_year']}-{meta['end_year']}, based on "
            f"{meta['dataset']}."
        )
    elif meta["diffmap_metric"] == "percent":
        prov["statistics"] = ["diff", "mean"]
        prov["caption"] = (
            f"Relative difference in {meta['long_name']} between first and "
            f"last {cfg['comparison_period']} years of the period "
            f"{meta['start_year']}-{meta['end_year']}, based on "
            f"{meta['dataset']}."
        )
    elif meta["diffmap_metric"] == "first":
        prov["caption"] = (
            f"Average {meta['long_name']} over the period "
            f"{meta['start_year']}-"
            f"{meta['start_year'] + cfg['comparison_period'] - 1}, based on "
            f"{meta['dataset']}."
        )
    elif meta["diffmap_metric"] == "last":
        prov["caption"] = (
            f"Average {meta['long_name']} over the period "
            f"{meta['end_year'] - cfg['comparison_period'] + 1}-"
            f"{meta['end_year']}, based on "
            f"{meta['dataset']}."
        )
    elif meta["diffmap_metric"] == "total":
        prov["caption"] = (
            f"Average {meta['long_name']} over the period "
            f"{meta['start_year']}-{meta['end_year']}, based on "
            f"{meta['dataset']}."
        )
    return prov


def plot_colorbar(
    cfg: dict,
    plotfile: str,
    plot_kwargs: dict,
    orientation: str = "vertical",
    mappable: mpl.cm.ScalarMappable | None = None,
) -> None:
    """Plot colorbar in its own figure for strip_plots."""
    _ = cfg  # we might need this in the future
    fig = plt.figure(figsize=(1.5, 3))
    # fixed size axes in fixed size figure
    cbar_ax = fig.add_axes([0.01, 0.04, 0.2, 0.92])
    if mappable is None:
        cmap = plot_kwargs.get("cmap", "RdYlBu")
        norm = mpl.colors.Normalize(
            vmin=plot_kwargs.get("vmin"),
            vmax=plot_kwargs.get("vmax"),
        )
        mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(
        mappable,
        cax=cbar_ax,
        orientation=orientation,
        label=plot_kwargs["cbar_label"],
        pad=0.0,
    )
    if "cbar_ticks" in plot_kwargs:
        cbar.set_ticks(plot_kwargs["cbar_ticks"], minor=False)
    fontsize = plot_kwargs.get("cbar_fontsize", 14)
    cbar.ax.tick_params(labelsize=fontsize)
    cbar.set_label(
        plot_kwargs["cbar_label"],
        fontsize=fontsize,
        labelpad=fontsize,
    )
    plotfile = plotfile.removesuffix(".png")
    fig.savefig(plotfile + "_cb.png")


def fill_era5_gap(meta: dict, cube: Cube) -> None:
    """Fill missing gap at 360 for era5 pet calculation."""
    if (
        meta["dataset"] == "ERA5"
        and meta["short_name"] == "evspsblpot"
        and len(cube.data[0]) == 360  # noqa: PLR2004
    ):
        cube.data[:, 359] = cube.data[:, 0]


def plot(
    cfg: dict,
    meta: dict,
    cube: Cube,
    basename: str,
    kwargs: dict | None = None,
) -> str:
    """Plot map using diag_scripts.shared module."""
    plotfile = get_plot_filename(basename, cfg)
    plot_kwargs = cfg.get("plot_kwargs", {}).copy()
    if kwargs is not None:
        plot_kwargs.update(kwargs)
    if "vmax" in plot_kwargs and "vmin" in plot_kwargs:
        plot_kwargs["levels"] = np.linspace(
            plot_kwargs["vmin"],
            plot_kwargs["vmax"],
            9,
        )
    label = plot_kwargs.get("cbar_label", "{short_name} ({units})")
    plot_kwargs["cbar_label"] = label.format(**meta)
    for coord in cube.coords(dim_coords=True):
        if not coord.has_bounds():
            log.warning("NO BOUNDS. GUESSING: %s", coord.name())
            cube.coord(coord.name()).guess_bounds()
    fill_era5_gap(meta, cube)
    add_cyclic_point(cube.data, cube.coord("longitude").points)
    mapplot = global_contourf(cube, **plot_kwargs)
    if cfg.get("clip_land", False):
        plt.gca().set_extent((220, 170, -55, 90))
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
    return plotfile


def apply_plot_kwargs_overwrite(
    kwargs: dict,
    overwrites: dict,
    metric: str,
    group: str,
) -> dict:
    """Apply plot_kwargs_overwrite to kwargs dict for selected plots."""
    for overwrite in overwrites:
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


def calculate_diff(cfg, meta, mm_data, output_meta, group) -> None:
    """Absolute difference between first and last years of a cube.

    Calculates the absolute and relative difference between the first and last
    period of a cube. Add data to the multi model collection ``mm_data`` and
    optionally plot each dataset.
    """
    cube = iris.load_cube(meta["filename"])
    if meta["short_name"] in cfg.get("convert_units", {}):
        pp.convert_units(cube, cfg["convert_units"][meta["short_name"]])
    if "start_year" in cfg or "end_year" in cfg:
        log.info("selecting time period")
        cube = pp.extract_time(
            cube,
            cfg["start_year"],
            1,
            1,
            cfg["end_year"],
            12,
            31,
        )
    dtime = cfg.get("comparison_period", 10) * 12
    cubes = {}
    cubes["total"] = cube.collapsed("time", MEAN)
    do_metrics = cfg.get("metrics", METRICS)
    norm = (
        int(meta["end_year"])
        - int(meta["start_year"])
        + 1  # count full end year
        - cfg.get("comparison_period", 10)  # decades center to center
    ) / 10
    cubes["first"] = cube[0:dtime].collapsed("time", MEAN)
    cubes["last"] = cube[-dtime:].collapsed("time", MEAN)
    cubes["diff"] = cubes["last"] - cubes["first"]
    cubes["diff"].data /= norm
    cubes["diff"].units = str(cubes["diff"].units) + " / 10 years"
    cubes["percent"] = cubes["diff"] / cubes["first"] * 100
    cubes["percent"].units = "% / 10 years"
    if cfg.get("plot_mmm", True):
        for key in do_metrics:
            mm_data[key].append(cubes[key])
    for key, cube in cubes.items():
        if key not in do_metrics:
            continue  # i.e. first/last if only diff is needed
        meta["diffmap_metric"] = key
        meta["exp"] = meta.get("exp", "exp")
        basename = cfg["basename"].format(**meta)
        meta["title"] = cfg["titles"][key].format(**meta)
        prov = _get_provenance(cfg, meta)
        prov["ancestors"] = [meta["filename"]]
        if cfg.get("plot_models", True):
            plot_kwargs = cfg.get("plot_kwargs", {}).copy()
            apply_plot_kwargs_overwrite(
                plot_kwargs,
                cfg.get("plot_kwargs_overwrite", []),
                key,
                group,
            )
            plotfile = plot(cfg, meta, cube, basename, kwargs=plot_kwargs)
            plt.close()
            utils.log_provenance(cfg, plotfile, prov)
        if cfg.get("save_models", True):
            meta["filename"] = str(Path(cfg["work_dir"]) / f"{basename}.nc")
            iris.save(cube, meta["filename"])
            output_meta[meta["filename"]] = meta.copy()
            utils.log_provenance(cfg, meta["filename"], prov)


def calculate_mmm(cfg, meta, mm_data, output_meta, group) -> None:
    """Calculate multi-model mean for a given metric."""
    for metric in cfg.get("metrics", METRICS):
        drop = cfg.get("dropcoords", ["time", "height"])
        meta = meta.copy()  # don't modify meta in place:
        meta["dataset"] = "MMM"
        meta["diffmap_metric"] = metric
        basename = cfg["basename"].format(**meta)
        mmm, _ = utils.mmm(
            mm_data[metric],
            dropcoords=drop,
            dropmethods=metric != "diff",
            mdtol=cfg.get("mdtol", 0.3),
        )
        meta["title"] = f"Multi-model Mean ({cfg['titles'][metric]})"
        if cfg.get("plot_mmm", True):
            plot_kwargs = cfg.get("plot_kwargs", {}).copy()
            overwrites = cfg.get("plot_kwargs_overwrite", [])
            apply_plot_kwargs_overwrite(plot_kwargs, overwrites, metric, group)
            plotfile = plot(cfg, meta, mmm, basename, kwargs=plot_kwargs)
            prov = _get_provenance(cfg, meta)
            prov["ancestors"] = meta["ancestors"]
            utils.log_provenance(cfg, plotfile, prov)
        if cfg.get("save_mmm", True):
            work_file = str(Path(cfg["work_dir"]) / f"{basename}.nc")
            meta["filename"] = work_file
            meta["diffmap_metric"] = metric
            output_meta[work_file] = meta.copy()
            iris.save(mmm, work_file)


def set_defaults(cfg: dict) -> None:
    """Update cfg with default values from diffmap.yml in place."""
    config_fpath = os.path.realpath(__file__)[:-3] + ".yml"
    with open(config_fpath, encoding="utf-8") as config_file:
        defaults = yaml.safe_load(config_file)
    for key, val in defaults.items():
        cfg.setdefault(key, val)
    if cfg["plot_kwargs_overwrite"] is not defaults["plot_kwargs_overwrite"]:
        cfg["plot_kwargs_overwrite"].extend(defaults["plot_kwargs_overwrite"])
    titles = defaults.get("titles", {})
    titles.update(cfg["titles"])
    cfg["titles"] = titles


def filter_metas(metas: list, filters: dict | list) -> list:
    """Filter metas by filter dicts."""
    if isinstance(filters, dict):
        filters = [filters]
    filtered = {}
    for selection in filters:
        for meta in select_metadata(metas, **selection):
            filtered[meta["filename"]] = meta
    return list(filtered.values())


def main(cfg) -> None:
    """Execute Diagnostic."""
    set_defaults(cfg)
    metas = cfg["input_data"].values()
    if cfg.get("filters") is not None:
        metas = filter_metas(metas, cfg["filters"])
    groups = group_metadata(metas, cfg["group_by"])
    output = {}
    for group, g_metas in groups.items():
        mm_data = defaultdict(list)
        for meta in g_metas:
            if "end_year" not in meta:
                meta.update(utils.get_time_range(meta["filename"]))
            # adjust norm for selected time period
            meta["end_year"] = cfg.get("end_year", meta["end_year"])
            meta["start_year"] = cfg.get("start_year", meta["start_year"])
            calculate_diff(cfg, meta, mm_data, output, group)
        do_mmm = cfg.get("plot_mmm", True) or cfg.get("save_mmm", True)
        if do_mmm and len(g_metas) > 1:
            # copy meta from first dataset and add all ancestors
            keep_keys = [
                "short_name",
                "long_name",
                "units",
                "start_year",
                "end_year",
                "exp",
            ]
            meta = {k: g_metas[0][k] for k in keep_keys}
            meta["ancestors"] = [met["filename"] for met in g_metas]
            meta["dataset"] = "MMM"
            calculate_mmm(cfg, meta, mm_data, output, group)
    utils.save_metadata(cfg, output)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
