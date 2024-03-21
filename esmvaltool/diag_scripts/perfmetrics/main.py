"""Overview plot for performance metrics.

Description
-----------
This diagnostic provides plot functionalities for performance metrics.
The multi model overview heatmap might be usefull for different
tasks and therefore this diagnostic tries to be as flexible as possible.
X and Y axis, grouping parameter and slits for each rectangle can be
configured in the recipe. All *_by parameters can be set to any metadata
key. To split by 'reference' this key needs to be set as extra_facet in recipe.

Author
------
Lukas Ruhe (Universität Bremen, Germany)
Diego Cammarano

Configuration parameters through recipe:
----------------------------------------
x_by: str, optional
    Metadata key for x coordinate.
    By default 'alias'.
y_by: str, optional
    Metadata key for y coordinate.
    By default 'variable_group'.
group_by: str, optional
    Metadata key for grouping.
    Grouping is always applied in x direction. Can be set to None to skip
    grouping into subplots.
    By default 'project'.
split_by: str, optional
    Not implemented yet.
    By default None.
plot_kwargs: dict, optional
    Dictionary that gets passed as kwargs to `matplotlib.pyplot.imshow()`.
    Colormaps will be converted to 11 discrete steps automatically. Default
    colormap is RdYlBu_r but can be changed with cmap.
    Other common keywords: vmin, vmax
    By default {}.
cbar_kwargs: dict, optional
    Dictionary that gets passed to `matplotlib.pyplot.colorbar()`.
    E.g. label, ticks...
    By default {}.
plot_properties: dict, optional
    Dictionary that gets passed to `matplotlib.axes.Axes.set()`.
    Subplots can be widely customized. For a full list of
    properties see:
    https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set.html
    E.g. xlabel, ylabel, yticklabels, xmargin...
    By default {}.
figsize: list(float), optional
   [width, height] of the figure in inches. The final figure will be saved with
   bbox_inches="tight", which can change the resulting aspect ratio.
   By default [5, 3].
"""

import itertools
import iris
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from esmvaltool.diag_scripts.shared import (
    get_plot_filename,
    group_metadata,
    run_diagnostic,
    select_metadata,
)
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.patches as patches


def unify_limits(cfg, grid):
    """set same limits for all subplots"""
    vmin, vmax = np.inf, -np.inf
    images = [ax.get_images()[0] for ax in grid]
    for im in images:
        vmin = min(vmin, im.get_clim()[0])
        vmax = max(vmax, im.get_clim()[1])
    for im in images:
        im.set_clim(vmin, vmax)


def plot_matrix(data, row_labels, col_labels, ax, plot_kwargs):
    """Create an image for given data."""
    im = ax.imshow(data, **plot_kwargs)
    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)
    # Rotate the tick labels and set their alignment.
    plt.setp(
        ax.get_xticklabels(),
        rotation=90,
        ha="right",
        va="center",
        rotation_mode="anchor",
    )
    # Turn spines off and create white grid.
    # ax.spines[:].set_visible(False)
    ax.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="black", linestyle="-", linewidth=1)
    ax.tick_params(which="both", bottom=False, left=False)
    return im


def remove_reference(metas):
    """ "return metadata that are not reference for metric"""
    return [
        meta for meta in metas if not meta.get("reference_for_metric", False)
    ]


def load_data(cfg, metas):
    """load nc files into dictionary of numpy arrays
    dict keys are the values for each split (i.e. reference).
    """
    splitted = group_metadata(metas, cfg["split_by"])
    y_labels = list(group_metadata(metas, cfg["y_by"]).keys())
    x_labels = list(group_metadata(metas, cfg["x_by"]).keys())
    data = {key: np.zeros((len(y_labels), len(x_labels))) for key in splitted}
    for x, x_label in enumerate(x_labels):
        for y, y_label in enumerate(y_labels):
            for split, split_metas in splitted.items():
                selection = {cfg["x_by"]: x_label, cfg["y_by"]: y_label}
                try:
                    meta = select_metadata(split_metas, **selection)[0]
                except IndexError:
                    print(f"No data found for {selection}")
                    data[split][y, x] = np.nan
                    continue
                cube = iris.load_cube(meta["filename"])
                data[split][y, x] = cube.data
    return data, x_labels, y_labels


def overlay_reference(cfg, ax, data, triangle):
    """create triangular overlays for given data and axes."""
    # use same colors as in main plot
    cmap = ax.get_images()[0].get_cmap()
    norm = ax.get_images()[0].norm
    # print("plotting overlay for", split)
    # print(data)
    # print(ax.get_images()[0])
    for i, j in itertools.product(*map(range, data.shape)):
        if np.isnan(data[i, j]):
            continue
        print(data[i, j])
        color = cmap(norm(data[i, j]))
        edges = [(e[0] + j, e[1] + i) for e in triangle]
        patch = patches.Polygon(
            edges,
            closed=True,
            facecolor=color,
            edgecolor="black",
            linewidth=0.5,
            fill=True,
        )
        ax.add_patch(patch)


def plot_additional_splits(cfg, grid, datas):
    """add references as triangular overlays"""
    # if len(datas) < 2:
    #     return
    # if len(datas) > 4:
    #     print("Too many splits for overlay, only 4 will be plotted.")
    half = [
        [(0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)],  # lower right
    ]
    quarters = [
        [(-0.5, -0.5), (0, 0), (-0.5, 0.5)],  # left
        [(0.5, -0.5), (0, 0), (0.5, 0.5)],  # right
        [(-0.5, 0.5), (0, 0), (0.5, 0.5)],  # bottom
        [(-0.5, -0.5), (0, 0), (0.5, -0.5)],  # top (?)
    ]
    if len(datas) == 2:
        quarters = half

    for ggg, data in enumerate(datas):
        data.pop(cfg["reference"])  # TODO: fix me
        for sss, (split, array) in enumerate(data.items()):
            # if split == cfg["reference"]:
            #     continue
            print("plotting overlay for", split, array.shape)
            overlay_reference(cfg, grid[ggg], array, quarters[sss])


def plot_group(cfg, ax, metas, title=None):
    """create matrix for one subplot in ax using plt.imshow()
    returns image object and dict of numpy array for each reference (split).
    """
    data, x_labels, y_labels = load_data(cfg, metas)
    # plot matrix with simple rectangles
    im = plot_matrix(
        data[cfg["reference"]], y_labels, x_labels, ax, cfg["plot_kwargs"]
    )
    if title is not None:
        ax.set_title(title)
    # overlay splits (additional references)
    # cmap might change after plot_group. set vmin and vmax before plotting?
    # if len(data) > 1:
    #     plot_additional_splits(cfg, ax, data)
    ax.set(**cfg["axes_properties"])
    return im, data


def plot(cfg, grouped_metas):
    """creates figure with subplots for each group"""
    fig = plt.figure(1, cfg.get("figsize", (5.5, 3.5)))
    grid = ImageGrid(
        fig,
        111,  # similar to subplot(111)
        cbar_mode="single",
        cbar_location="right",
        cbar_pad=0.1,
        cbar_size=0.2,
        nrows_ncols=(1, len(grouped_metas)),
        axes_pad=0.1,
    )
    # remap colorbar to 10 discrete steps
    cmap = mpl.cm.get_cmap(cfg.get("cmap", "RdYlBu_r"), 10)
    cfg["plot_kwargs"]["cmap"] = cmap
    datas = []
    for i, (group, metas) in enumerate(grouped_metas.items()):
        title = group if len(grouped_metas) > 1 else None
        im, data = plot_group(cfg, grid[i], metas, title=title)
        datas.append(data)  # collect data for overlay
    # use same colorrange and colorbar for all subplots:
    unify_limits(cfg, grid)
    grid.cbar_axes[0].colorbar(im, **cfg["cbar_kw"])
    # cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # plot overlay with same colors
    plot_additional_splits(cfg, grid, datas)
    basename = "performance_metrics.png"
    fname = get_plot_filename(basename, cfg)
    # plt.tight_layout()
    plt.savefig(fname, bbox_inches="tight")
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
    cfg.setdefault("x_by", "alias")
    cfg.setdefault("y_by", "variable_group")
    cfg.setdefault("group_by", "project")
    cfg.setdefault("split_by", None)
    cfg.setdefault("cbar_kw", {})
    cfg.setdefault("axes_properties", {})
    cfg.setdefault("plot_kwargs", {})
    cfg["plot_kwargs"].setdefault("cmap", "RdYlBu_r")
    cfg["plot_kwargs"].setdefault("vmin", 0)
    cfg["plot_kwargs"].setdefault("vmax", 1)


def main(cfg):
    """run the diagnostic"""
    set_defaults(cfg)
    metas = cfg["input_data"].values()
    metas = remove_reference(metas)
    grouped_metas = group_metadata(metas, cfg["group_by"])
    cfg["figsize"] = (7.5, 5.5)
    cfg["y_by"] = "short_name"
    cfg[
        "split_by"
    ] = "ref"  # TODO: hardcode split as extra facet? and use None as default?
    cfg[
        "reference"
    ] = "ref1"  # TODO: this becomes obsolet if we use extra facet
    plot(cfg, grouped_metas)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
