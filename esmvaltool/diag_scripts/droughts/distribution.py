#!/usr/bin/env  python
# -*- coding: utf-8 -*-
"""Creates histograms and regional boxplots for given timeperiods.

Global histograms are plotted to compare distributions of any variable in given
intervals and/or by experiment. Combined experiments (historical-ssp*) will
be splitted into individual ones. 

NOTE: This diagnostic loads all values from all datasets into memory. If you
provide a lot of data make sure enough memory is available.

TODO: Free up memory each split for global histograms only save bins and fit.

Configuration options in recipe
-------------------------------
group_by: str, optional (default: short_name)
    All input datasets are grouped by this metadata key. The diagnostic runs
    for each of this groups. Consider adding group_by to the basename.
split_experiments: bool, optional (default: False)
    If true, combined experiments like "historical-ssp126" get treated
    individually as "historical" and "ssp126" when grouped.
    NOTE: Only works for historical + scenarioMIP (cutted at 2014) yet.
split_by: str, optional (default: exp)
    Create individual distributions for each split in the same figure.
    Can be any metadata key or 'interval'. For 'interval' the corresponding
    parameter must be given.
    TODO: For 'exp' consider changing split_experiments.
intervals: list of dicts, optional (default: [])
    List of dicts containing a `label` (optional) and a `range` 
    (timerange in ISO 8601 format. For example YYYY/YYYY) or 
    `start` and `end` (ISO 8601).
    In the diagnostics `interval_label` and `interval_range` will be added
    to metadata and can be used as split_by option and in basename.
    NOTE: if splitted by experiment, the first interval is used for historical,
    the second for scenarioMIP. For other splits this option is ignored.
regions: list of str, optional (default: [])
    List of AR6 regions to extract data for histograms (all given regions)
    and boxplots (each region).
sort_regions_by: str, optional (default: "")
    Sort regions by median value of given split_by value. If empty no sorting.
basename: str, optional (default: '{plot_type}_{group}')
    Filename for plot files. Can contain placeholders for group, plot_type
    and interval_label or any other metadata key.
comparison_period: int, optional (default: 10)
    Number of years to compare (first and last N years). Must be less or equal
    half of the total time period.
plot_mmm: bool, optional (default: True)
    Calculate and plot the average over all datasets.
plot_models: bool, optional (default: True)
    Plot maps for each dataset.
plot_properties: dict, optional (default: {})
    Kwargs passed to all axes.set() functions. Can be styling of ticks, labels,
    limits, grid, etc. See matplotlib documentation for all options:
    https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set.html
    Use `histogram.plot_properties` or `regional_stats.plot_properties`
    to specify kwargs for corresponding plots.
plot_kwargs: dict, optional (default: {})
    Kwargs to all plot functions. Use `histogram.plot_kwargs` or 
    `regional_stats.plot_kwargs` to specify kwargs for corresponding plots.
regional_stats.fig_kwargs: dict, optional (default: 
    {'figsize': (15, 3), 'dpi': 300})
    Kwargs passed to the figure function for the scenarios plot. The best size
    for the figure depends on the number of splits and regions and the wanted
    aspect ratio. 15,3 is the default for a large number of regions.
colors: dict, optional (default: {})
    Define colors as dict. Keys are the split_by values. For experiment or
    dataset splits colors are used from ipcc colors, if available.
    TODO: add ipcc dataset colors.
labels: dict, optional (default: {})
    Define labels as dict. Keys are the split_by values. Values are strings 
    shown in the legend.
strip_plots: bool, optional (default: False)
    Removes titles, margins and colorbars from plots (to use them in panels).
grid_lines: list of float, optional (default: [])
    Draw helper lines to the plots to indicate threshholds or ranges. At
    given locations hlines are drawn in boxplots and vlines in histograms.
grid_line_style: dict, optional (default: 
    {'color': 'gray', 'linestyle': '--', 'linewidth': 0.5})
    Style of the grid lines. See matplotlib documentation for all options.
"""

import iris
import yaml
import copy
import numpy as np
import xarray as xr
from cycler import cycler
from esmvalcore import preprocessor as pp
# from matplotlib.lines import Line2D
from matplotlib import cbook
from iris.analysis.cartography import area_weights
from iris.time import PartialDateTime
import logging
from collections import defaultdict
import matplotlib.pyplot as plt
from esmvaltool.diag_scripts.shared import (
    # get_plot_filename,
    group_metadata,
    run_diagnostic,
    get_diagnostic_filename,
)
from scipy.stats import norm
from esmvaltool.diag_scripts.droughtindex import (
    colors as ipcc_colors,
    utils as ut,
)

log = logging.getLogger(__file__)


SER_KEYS = ["mean", "med", "q1", "q3", "iqr", 
            "whislo", "whishi", "cihi", "cilo"]
BOX_PLOT_KWARGS = {}

def load_constrained(cfg, meta, regions=[]):
    """Load cube with constraints in meta."""
    cube = iris.load_cube(meta["filename"], constraint=meta.get("cons", None))
    ut.guess_lat_lon_bounds(cube)
    if regions != []:
        log.debug("Extracting regions: %s", regions)
        cube = pp.extract_shape(cube, shapefile="ar6", ids={"Acronym": regions})
    if "interval_range" in meta:
        log.debug("clipping timerange to %s", meta["interval_range"])
        cube = pp.clip_timerange(cube, meta["interval_range"])
    return cube


def group_by_interval(cfg, metas: list):
    """Group metadata by intervals."""
    groups = defaultdict(list)
    for interval in cfg["intervals"]:
        imetas = copy.deepcopy(metas)
        for m in imetas:
            m["interval_range"] = interval["range"]
            m["interval_label"] = interval["label"]
        groups[interval["label"]] = imetas
    return groups


def group_by_exp(cfg, metas, historical_first=False):
    """similar to shared.group_metadata but splits combined experiments.
    meta data will be added to individual experiments.
    To keep this function lazy an iris constraint is added (`meta["cons"]`),
    to be applied when loading the file.
    If comparison_period is given, last part of its length in each experiment
    is used. For historical the first part is used if `historical_first=True`.
    """
    groups = group_metadata(metas, "exp")
    historical_added = []
    year = PartialDateTime(year=2015)
    groups = defaultdict(list, groups)
    remove_keys = set()
    # split combined experiments:
    for exp in list(groups.keys()):
        metas = groups[exp]
        if not exp.startswith("historical-"):
            continue
        if exp not in historical_added:
            for meta in metas:
                hmeta = copy.deepcopy(meta)
                hmeta["cons"] = iris.Constraint(time=lambda cell: cell < year)
                hmeta["exp"] = "historical"
                if len(cfg["intervals"]) > 0:
                    log.warning("set interval range for historical")
                    hmeta["interval_range"] = cfg["intervals"][0]["range"]
                groups["historical"].append(hmeta)
            historical_added.append(exp)
        for meta in metas:
            meta["cons"] = iris.Constraint(time=lambda cell: cell >= year)
            meta["exp"] = exp.split("-")[1]
            if len(cfg["intervals"]) > 0:
                log.warning("set interval range for ssp")
                meta["interval_range"] = cfg["intervals"][1]["range"]
            groups[meta["exp"]].append(meta)
            remove_keys.add(exp)
    for rmkey in remove_keys:
        del groups[rmkey]
    return groups



def calculate_histogram(cfg, splits, output, group):
    """load data for each split and calculate counts, bins and fit parameters.
    Safe parameters to netcdf file, to optionally skip this part on rerun.
    """
    labels = []
    # weights = []
    fits = []
    counts = []
    bins = []
    log.info("start split loading")
    for split, metas in splits.items():
        labels.append(cfg.get("labels", {}).get(split, split))
        split_data = []
        split_weights = []
        for meta in metas:  # merge everything else
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                cube = load_constrained(cfg, meta, regions=cfg["regions"])
                split_weights.append(area_weights(cube))
            applied_mask = cube.data.filled(np.nan)  # type: ignore
            split_data.append(applied_mask.ravel())
        flat = np.array(split_data).ravel()
        flat_weights = np.array(split_weights).ravel()
        fits.append(norm.fit(flat[~np.isnan(flat)]))  # mu, std
        mybins = np.arange(-8, 6.5, step=0.5)
        s_counts, s_bins = np.histogram(flat, mybins, weights=flat_weights)
        split_bins_center = (s_bins[1:] + s_bins[:-1]) / 2
        counts.append(s_counts)
        bins.append(split_bins_center)
        log.info("split done")
    # save output
    data = xr.Dataset({
        "fits": xr.DataArray(fits, dims=["split", "param"], name="fit"),
        "counts": xr.DataArray(counts, dims=["split", "bin"], name="counts"),
        "bins": xr.DataArray(bins, dims=["split", "bin"], name="bins"),
        "labels": xr.DataArray(labels, dims=["split"], name="labels")
    })
    fname = get_diagnostic_filename(f"histogram_{group}", cfg)
    output["fname"] = {"filename": fname, "plottype": "histogram", "group": group}
    data.to_netcdf(fname)
    return data


def load_histogram(cfg, group):
    """Load histogram data from netcdf file."""
    fname = get_diagnostic_filename(f"histogram_{group}", cfg)
    return xr.open_dataset(fname)


def plot_histogram(cfg, splits, output, group, fit=True):
    """Plot one combined histogram of all given regions for each split."""
    plt.figure(**ut.sub_cfg(cfg, "histogram", "fig_kwargs"))
    colors = list(get_split_colors(cfg, splits).values())
    if ut.sub_cfg(cfg, "histogram", "reuse"):
        data = load_histogram(cfg, group)
    else:
        data = calculate_histogram(cfg, splits, output, group)
    bins = data["bins"].values.T
    counts = data["counts"].values.T
    labels = data["labels"].values.T
    plot_kwargs = {
        "bins": np.arange(-8, 6.5, step=0.5),
        "label": labels,
        "density": True,
        "alpha": 0.75,
    }
    log.info("plotting histogram (%s)", group)
    _, bins, patches = plt.hist(bins, weights=counts, color=colors, **plot_kwargs, zorder=3)
    legend = plt.legend()
    for patch in legend.get_patches():
        patch.set_alpha(1)
    plot_props = ut.sub_cfg(cfg, "histogram", "plot_properties")
    plt.gca().set(**plot_props)
    for line in cfg["grid_lines"]:
        plt.axvline(line, **cfg["grid_line_style"])

    meta = next(iter(splits.values()))[0].copy()  # first meta from dict
    meta.update({
        "plot_type": "histogram",
        "group": group,
    })
    filename = ut.get_plot_filename(cfg, cfg["basename"], meta, {"/": "_"})
    plt.savefig(filename)
    log.info("saved %s", filename)

    # plot normal fit
    plot_kwargs_fit = {
        "linewidth": 2,
        "linestyle": "-",
        "alpha": 1,
    }
    for patch in patches:  # type: ignore
        for rect in patch:
            rect.set_alpha(0.2)
    for iii, fit in enumerate(data["fits"].values):
        x = np.linspace(bins[0], bins[-1], 200)
        p = norm.pdf(x, fit[0], fit[1])
        plot_kwargs_fit.update(cfg["histogram"].get("plot_kwargs", {}))
        plt.plot(x, p, color=colors[iii], **plot_kwargs_fit)
    # fit_line = Line2D([], [], color="black", linestyle="-", alpha=1)
    # handles, labels = plt.gca().get_legend_handles_labels()
    # handles.append(fit_line)
    # plt.legend(handles=handles, labels=labels + ["normal fit"])
    # ax.legend(handles, labels)
    meta["plot_type"] = "histogram_fit"
    filename = ut.get_plot_filename(cfg, cfg["basename"], meta, {"/": "_"})
    log.info("saved %s", filename)
    plt.savefig(filename)
    plt.close()


def sort_regions(data, regions, by="ssp585", inplace=True):
    """Sort regions by median value."""
    if not inplace:
        data = copy.deepcopy(data)
    medians = [np.nanmedian(np.array(reg).ravel()) for reg in data[by]]
    order = np.argsort(medians)
    regions = [regions[i] for i in order]
    for split in data.keys():
        data[split] = [data[split][i] for i in order]
    return data, regions


def calculate_regional_stats(cfg, splits, output, group):
    """Calculate regional statistics for given metadata."""
    data = {}
    stats = {}
    regions = cfg["regions"]
    if regions == []:
        regions = list(ut.get_hex_positions().keys())
    for split, metas in splits.items():  # default each experiment
        region_data = defaultdict(list)
        for meta in metas:
            cube = load_constrained(cfg, meta)
            ut.guess_lat_lon_bounds(cube)
            for region in regions:
                ids = {"Acronym": [region]}
                reg = pp.extract_shape(cube, shapefile="ar6", ids=ids)
                flat = reg.data[~reg.data.mask]  # cheaper  # type: ignore
                region_data[region].extend(flat)
        data[split] = [region_data[r] for r in regions]
    if cfg["sort_regions_by"]:
        data, regions = sort_regions(data, regions, by=cfg["sort_regions_by"])
    for split, dat in data.items():
        stats[split] = cbook.boxplot_stats(dat, whis=(2.3, 97.7), labels=regions)
    # save stats
    fname = get_diagnostic_filename(f"regional_stats_{group}", cfg)
    fname = fname.replace(".nc", ".yml")
    output["fname"] = {
        "filename": fname, "plottype": "regional_stats", "group": group
    }
    for split_stats in stats.values():
        for stat in split_stats:
            del stat["fliers"]
            for key in SER_KEYS:
                stat[key] = stat[key].tolist()
    with open(fname, "w") as fstream:
        yaml.dump(stats, fstream)
    log.info("saved: %s", fname)
    return stats


def load_regional_stats(cfg, group):
    """Load regional statistics from netcdf file."""
    fname = get_diagnostic_filename(f"regional_stats_{group}", cfg)
    fname = fname.replace(".nc", ".yml")
    with open(fname, "r") as f:
        stats = yaml.load(f, Loader=yaml.SafeLoader)
    for split in stats:
        for stat in stats[split]:
            for key in SER_KEYS:
                stat[key] = np.array(stat[key])
    return stats


def plot_regional_stats(cfg, splits, output, group):
    if ut.sub_cfg(cfg, "regional_stats", "reuse"):
        log.info("loading boxplot data from file")
        stats = load_regional_stats(cfg, group)
    else:
        stats = calculate_regional_stats(cfg, splits, output, group)
    regions = [sdata["label"] for sdata in list(stats.values())[0]]
    fig = plt.figure(**ut.sub_cfg(cfg, "regional_stats", "fig_kwargs"))
    positions = np.array(np.arange(len(regions)))
    colors = get_split_colors(cfg, splits)
    for line in cfg["grid_lines"]:
        plt.hlines(line, -1, len(regions), **cfg["grid_line_style"])
    for i, (split, stat) in enumerate(stats.items()):
        n = len(stats)
        group_width = 0.9
        width = group_width / (n + 1)
        offset = ((i + 1) / (n + 1) - 0.5) * group_width
        plot_kwargs_box = {
            "widths": width,
            "positions": positions + offset,
            "showfliers": False,
            "patch_artist": True,
            "medianprops": {"color": "white", "linewidth": 1},
            "boxprops": {
                "facecolor": colors[split],
                "edgecolor": "white",
                "linewidth": 1,
            },
            "whiskerprops": {"color": colors[split], "linewidth": 0.8},
            "capprops": {"color": colors[split], "linewidth": 1},
        }
        # plt.boxplot(sdata, **plot_kwargs_box)
        plt.gca().bxp(stats[split], **plot_kwargs_box)
        plt.plot([], c=colors[split], label=split)  # dummy for legend
    fig.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.xticks(positions, regions, rotation=90)
    plt.tight_layout()
    plt.xlim(-1, len(regions))
    # save plot
    fname = ut.get_plot_filename(cfg, f"regional_stats_{group}")
    plt.savefig(fname)
    log.info("saved %s", fname)


def get_split_colors(cfg, splits):
    """adds colors for each split to the config if not yet present.
    ipcc colors are used if available, matplotlib defaults otherwise.
    TODO: add ipcc model colors.
    """
    colors = cfg.get("colors", {}).copy() # new instance for each group
    mpl_default = plt.rcParams["axes.prop_cycle"].by_key()['color']
    cycle = iter(cycler(color=mpl_default))
    missing_splits = [s for s in splits if s not in colors]
    for split in missing_splits:
        colors[split] = next(cycle)
        if cfg["split_by"] == "exp" and hasattr(ipcc_colors, split):
            colors[split] = getattr(ipcc_colors, split)
    return colors


def set_defaults(cfg):
    cfg.setdefault("group_by", "short_name")
    cfg.setdefault("split_by", "exp")
    cfg.setdefault("intervals", [])
    cfg["intervals"] = [ut.fix_interval(i) for i in cfg["intervals"]]
    cfg.setdefault("basename", "{plot_type}_{group}")
    cfg.setdefault("regions", [])
    cfg.setdefault("sort_regions_by", "")
    cfg.setdefault("plot_models", False)
    cfg.setdefault("plot_mmm", True)
    cfg.setdefault("plot_global", True)
    cfg.setdefault("plot_properties", {})
    cfg.setdefault("plot_kwargs", {})
    cfg.setdefault("histogram", {})
    cfg.setdefault("regional_stats", {})
    cfg["regional_stats"].setdefault("fig_kwargs", 
                                     {"figsize": (15, 3), "dpi": 300})
    cfg.setdefault("reuse", False)
    cfg.setdefault("strip_plots", False)
    cfg.setdefault("grid_lines", [])
    cfg.setdefault("grid_line_style", {
        "color": "gray", "linestyle": "--", "linewidth": 0.5
    })


def main(cfg):
    """main function. executing all plots for each group."""
    set_defaults(cfg)
    groups = group_metadata(cfg["input_data"].values(), cfg["group_by"])
    output = {}
    for group, gmetas in groups.items():
        log.info("running diagnostic for: %s", group)
        if cfg["split_by"] == "exp":
            splits = group_by_exp(cfg, gmetas)
        elif cfg["split_by"] == "interval":
            splits = group_by_interval(cfg, gmetas)
        else:
            splits = group_metadata(gmetas, cfg["split_by"])
        cfg["split_colors"] = get_split_colors(cfg, splits)
        if cfg["plot_global"]:
            log.info("plotting global histogram")
            plot_histogram(cfg, splits, output, group)
            plt.close()
        if cfg["plot_regions"]:
            log.info("plotting regional statistics")
            plot_regional_stats(cfg, splits, output, group)
            plt.close()
    ut.save_metadata(cfg, output)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
