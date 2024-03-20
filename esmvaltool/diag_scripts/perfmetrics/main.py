# This diagnostic provides plot functionalities 
# for performance metrics.
# The multi model overview heatmap might be usefull for different
# tasks and therefore this diagnostic tries to be as flexible as possible.
# X and Y axis, grouping parameter and slits for each rectangle can be
# configured in the recipe. All *_by parameters can be set to any metadata
# key. To split by 'reference' this key needs to be set as extra_facet in recipe.
# NOTE: should different references be done by aliases or specific extra_facet in recipe?
# 
# x_by: [alias]
# y_by: [variable_group]
# group_by: [project] gaps always applied in x direction
# split_by: [None]
# plot_kwargs: {}
# cbar_kwargs: {}


import matplotlib as mpl
from esmvalcore import preprocessor as pp
import itertools
import logging
import matplotlib.pyplot as plt
import esmvaltool.diag_scripts.shared as e
from esmvaltool.diag_scripts.shared import (
    # ProvenanceLogger,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
    select_metadata,
)

import iris
import numpy as np
import matplotlib.pyplot as plt



def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw=None, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if ax is None:
        ax = plt.gca()

    if cbar_kw is None:
        cbar_kw = {} 

    im = ax.imshow(data, **kwargs)
    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
        va="center", rotation_mode="anchor")

    # Turn spines off and create white grid.
    # ax.spines[:].set_visible(False)
    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="black", linestyle='-', linewidth=1)
    ax.tick_params(which="both", bottom=False, left=False)
    return ax, cbar


def remove_reference(metas):
    """"return metadata that are not reference for metric"""
    return [meta for meta in metas 
        if not meta.get("reference_for_metric", False)]


def plot(cfg, metas, groups=None):
    # print(metas)
    y_labels = list(group_metadata(metas, cfg["y_by"]).keys())
    x_labels = list(group_metadata(metas, cfg["x_by"]).keys())
    # print(x_labels)
    data = np.zeros((len(y_labels), len(x_labels)))
    if groups is not None:
        # add empty column as spacer
        data = np.zeros((len(y_labels), len(x_labels)+len(groups)-1))
    print(data)
    # add empty labels at certain positions for gaps
    # drop last position (not needed)
    print(groups)
    gap_positions = list(groups.values())[:-1]
    gap_positions.reverse()
    for gap in gap_positions:
        x_labels.insert(gap, "")
    data = np.zeros((len(y_labels), len(x_labels)))
    print(x_labels)
    for x, x_label in enumerate(x_labels): 
        if x_label == "":
            data[:, x] = 0
            continue
        for y, y_label in enumerate(y_labels):
            selection = {cfg["x_by"]: x_label, cfg["y_by"]: y_label}
            meta = select_metadata(metas, **selection)[0]
            cube = iris.load_cube(meta["filename"])
            data[y, x] = cube.data
            # if x == list(groups.values())[gaps]: print("GAP ADDED at ", x)
            #     gaps += 1
            #     data[x+gaps][y] = 0

    cbar_kw = {"extend": "both"}
    cbar_kw.update(cfg.get("cbar_kw", {}))
    ax, cbar = heatmap(data, y_labels, x_labels, **cfg["plot_kwargs"], cbar_kw=cbar_kw)

    # hide gap by drawing white rectangle
    for x, x_label in enumerate(x_labels):
        if x_label == "":
            ax.add_patch(mpl.patches.Rectangle((x-0.45, 0-0.55), 0.9, 9, fill=True, color="white", zorder=10))
    # ax.add_patch(mpl.patches.Rectangle((0, 0), 1, 1, fill=True, color="white", zorder=10))
    # set any axis properties
    ax.set(**cfg["axes_properties"])
    basename = "performance_metrics.png"
    fname = get_plot_filename(basename, cfg)
    plt.savefig(fname)
    print("Figure saved:")
    print(fname)


def apply_grouping(cfg, metas):
    """returns sorted metadata, and group labels with positions"""
    group_by = cfg.get("group_by", "project")
    grouped = group_metadata(metas, group_by)
    counts = []
    labels = []
    for group, metas in grouped.items():
        labels.append(group)
        x_entries = len(group_metadata(metas, cfg["x_by"]))
        counts.append(x_entries)
    positions = np.cumsum(counts)
    groups = dict(zip(labels, positions))
    metas = list(itertools.chain(*grouped.values()))
    return metas, groups


def set_defaults(cfg):
    """set default values for most important config parameters"""
    cfg.setdefault("x_by", "dataset")
    cfg.setdefault("y_by", "variable_group")
    cfg.setdefault("group_by", "project")
    cfg.setdefault("split_by", None)
    cfg.setdefault("axes_properties", {})
    cfg.setdefault("plot_kwargs", {})
    cfg["plot_kwargs"].setdefault("cmap", "RdYlBu_r")


def main(cfg):
    set_defaults(cfg)
    metas = cfg["input_data"].values()
    metas = remove_reference(metas)
    metas, groups = apply_grouping(cfg, metas)
    plot(cfg, metas, groups=groups)
    # grouped = group_metadata(metas, "project")
    # for group, group_metas in grouped.items():
    #     plot_group()


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)