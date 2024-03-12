#!/usr/bin/env  python
# -*- coding: utf-8 -*-
"""Creates histograms for first and last given interval of any variable.

Global histograms are plotted to compare distributions of any variable in given
intervals and/or by experiment. Combined experiments (historical-ssp*) will
be splitted into individual ones. 


Configuration options in recipe
-------------------------------
plot_mmm: bool, optional (default: True)
    Calculate and plot the average over all datasets.
plot_models: bool, optional (default: True)
    Plot maps for each dataset.
plot_properties: dict, optional (default: {})
    Kwargs passed to the axes.set() function. Can be styling of ticks, labels,
    limits, grid, etc. See matplotlib documentation for all options:
    https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set.html
plot_kwargs: dict, optional (default: {})
    Kwargs to the plot function.
plot_kwargs_overwrite: list, optional (default: [])
    List of plot_kwargs dicts for specific metrics (diff, first, latest, total)
    and group_by values (ie. pr, tas for group_by: short_name).
    `group` and `metric` can either be strings or lists of strings to be
    applied to all matching plots. Leave any of them empty to apply to all.
    All other given keys are applied to the plot_kwargs dict for this plot.
    Settings will be applied in order of the list, so later entries can
    overwrite previous ones.
colors: dict, optional (default: {})
    Define colors as dict. Keys are the split_by values. For experiment or
    dataset splits colors are used from style.py, if available.
labels: dict, optional (default: {})
    Define labels as dict. Keys are the split_by values. Values are strings 
    shown in the legend.
comparison_period: int, optional (default: 10)
    Number of years to compare (first and last N years). Must be less or equal
    half of the total time period.
group_by: str, optional (default: short_name)
    Meta key to loop over for multiple datasets.
split_experiments: bool, optional (default: False)
    If true histograms are plotted split by experiment. If comparison period is
    given as well only last period from each experiment is used.
strip_plots: bool, optional (default: False)
    Removes titles, margins and colorbars from plots (to use them in panels).
"""

import iris
import copy
import numpy as np
import os
import matplotlib as mpl
from esmvalcore import preprocessor as pp
from matplotlib.lines import Line2D
from iris.analysis.cartography import area_weights
import logging
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import iris.plot as ipl
from esmvaltool.diag_scripts.shared import (
    get_plot_filename,
    group_metadata,
    run_diagnostic
)
from scipy.stats import norm
# from esmvaltool.diag_scripts.droughtindex import (
#     colors as ipcc_colors,
#     utils as ut
# )
from collections import defaultdict
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

logger = logging.getLogger(__file__)


def load_constrained(meta, comparison_period=None):
    """Load cube with constraints in meta."""
    cube = iris.load_cube(meta["filename"], constraint=meta.get("cons", None))
    if comparison_period is not None:
        cube = cube[-comparison_period*12:]
    return cube


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
    year = iris.time.PartialDateTime(year=2015)
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
                groups["historical"].append(hmeta)
            historical_added.append(exp)
        for meta in metas:
            meta["cons"] = iris.Constraint(time=lambda cell: cell >= year)
            meta["exp"] = exp.split("-")[1]
            groups[meta["exp"]].append(meta)
            remove_keys.add(exp)
    for rmkey in remove_keys:
        del groups[rmkey]
    return groups


def plot_histogram(cfg, splits, output, group, fit=True):
    data = []
    weights = []
    colors = []
    labels = []
    fits = []
    for split, metas in splits.items():
        labels.append(cfg.get("labels", {}).get(split, split))
        split_data = []
        split_weights = []
        for meta in metas:  # merge everything else
            cube = iris.load_cube(meta["filename"])
            # ut.guess_lat_lon_bounds(cube)
            split_weights.append(area_weights(cube))
            applied_mask = cube.data.filled(np.nan)  # set masked to nan
            split_data.append(applied_mask.ravel())
        if cfg.get("split_by", "exp") == "exp":
            colors =  [ "blue", "green"] #!
            # color = getattr(ipcc_colors, split)
            # colors.append(cfg.get("colors", {}).get(split, color))
        else:
            colors =  [ "blue", "green"] 
        # TODO: elif dataset colors
        flat = np.array(split_data).ravel()
        fits.append(norm.fit(flat[~np.isnan(flat)]))  # mu, std
        data.append(flat)
        weights.append(np.array(split_weights).ravel())
    # np_data = np.array(data)
    plot_kwargs = {
        "bins": np.arange(0, 1e-1, step=5e-3),
        "label": labels,
        # "color": colors,
        "density": True,
        "weights": weights,
        "alpha": 0.75
    }
    n, bins, patches = plt.hist(data, **plot_kwargs, zorder=3)
    plt.legend()
    # plt.gca().set(**cfg.get("plot_properties", {}))
    plt.gca().set(ylim=(0, 100))
    # For log-scale
    # plt.gca().set(xscale='log') 
    # plt.gca().set(yscale='log')
    roi = cfg.get("highlight_region", [-2, +2])
    if roi is not None:
        ylim = plt.gca().get_ylim()[1]
        plt.gca().patch.set_facecolor('#f2f2f2')
        rois = PatchCollection(
            [Rectangle((roi[0], 0), roi[1] - roi[0], ylim)],
            facecolor="white", edgecolor="white", zorder=0)
        plt.gca().add_collection(rois)
    basename = f"histogram_{group}"
    plt.savefig(get_plot_filename(basename, cfg))

    # plot normal fit
    plot_kwargs_fit = {
        "linewidth": 2,
        "linestyle": "-",
        "alpha": 1,
    }
    for patch in patches:
        try:
            for rect in patch:
                rect.set_alpha(0.25)
        except:
            print('patch is already a rectangle itself')
            patch.set_alpha(0.25)
    for iii, fit in enumerate(fits):
        x = np.linspace(bins[0], bins[-1], 100)
        p = norm.pdf(x, fit[0], fit[1])
        plt.plot(x, p, color=colors[iii], **plot_kwargs_fit)
    fit_line = Line2D([], [], color="black", linestyle="-", alpha=1)
    handles, labels = plt.gca().get_legend_handles_labels()
    handles.append(fit_line)
    plt.legend(handles=handles, labels=labels + ["normal fit"])
    # ax.legend(handles, labels)
    basename = f"histogram_{group}_fit"
    plt.savefig(get_plot_filename(basename, cfg))

def main(cfg):
    """Main function."""
    cfg["group_by"] = cfg.get("group_by", "short_name")
    groups = group_metadata(cfg["input_data"].values(), cfg["group_by"])
    output = {}
    for group, gmetas in groups.items():
        if cfg.get("split_by", "exp") == "exp":
            splits = group_by_exp(cfg, gmetas)
        else:
            splits = group_metadata(gmetas, cfg["split_by"])
        if cfg.get("plot_global", True):
            plot_histogram(cfg, splits, output, group)
        # if cfg.get("plot_regions", True):
        #     regional_stats(cfg, splits, output, group)
    # ut.save_metadata(cfg, output)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)