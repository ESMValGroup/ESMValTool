"""Overview plot for performance metrics.

Description
-----------
This diagnostic provides plot functionalities for performance metrics.
The multi model overview heatmap might be useful for different
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
    The rectangles can be split into 2-4 triangles. This is used to show
    metrics for different references. For this case there is no need to change
    this parameter. Multiple variables can be set in the recipe with `split`
    assigned as extra_facet to label the different references. Data without
    a split assigned will be plotted as main rectangles, this can be changed
    by setting default_split parameter.
    By default 'split'.
default_split: str, optional
    Data labeled with this string, will be used as main rectangles. All other
    splits will be plotted as overlays. This can be used to choose the base
    reference, while all references are labeled for the legend.
legend: dict, optional
    Customize, if, how and where the legend is plotted. The 'best' position
    and size of the legend depends on multiple parameters of the figure
    (i.e. lengths of labels, aspect ratio of the plots...). And might require
    manual adjustment of `x`, `y` and `size` to fit the figure layout.
    Keys (each optional) that will be handled are:
    position: str or None, optional
        Position of the legend. Can be 'right' or 'left'. Or set to None to
        disable plotting the legend. By default 'right'.
    x_offset: float, optional
        Manually adjust horizontal position to save space or fix overlap.
        Number given in Inches. By default 0.
    y_offset: float, optional
        Manually adjust vertical position to save space or fix overlap.
        Number given in Inches. By default 0.
    size: float, optional
        Size of the legend in Inches. By default 0.3.
    TODO: set defaults and add offsets to legend function.
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
import logging

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib import patches
from mpl_toolkits.axes_grid1 import ImageGrid

from esmvaltool.diag_scripts.shared import (
    get_plot_filename,
    group_metadata,
    run_diagnostic,
    select_metadata,
)

log = logging.getLogger(__name__)


def unify_limits(cfg, grid):
    """Set same limits for all subplots."""
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
    """Remove reference for metric from list of metadata setting split=none
    allows to omit it in recipe, but handle it as special case of split data.

    list() creates a copy with same references to allow removing in
    place
    """
    for meta in list(metas):
        if meta.get("reference_for_metric", False):
            metas.remove(meta)


def add_split_none(metas):
    """List of metadata with split=None if no split is given."""
    for meta in metas:
        if "split" not in meta:
            meta["split"] = None


def open_file(metadata, **selection):
    """Try to find a single file for selection and return data.

    If multiple files are found, raise an error. If no file is found,
    return np.nan.
    """
    metas = select_metadata(metadata, **selection)
    if len(metas) > 1:
        raise ValueError(f"Multiple files found for {selection}")
    if len(metas) < 1:
        log.debug("No Metadata found for %s", selection)
        return np.nan
    log.warning("Metadata found for %s", selection)
    ds = xr.open_dataset(metas[0]["filename"])
    varname = list(ds.data_vars.keys())[0]
    return ds[varname].values.item()
    # iris.load_cube(metas[0]["filename"]).data


def load_data(cfg, metas):
    """Load all nc files from metadata into xarray dataset.

    The dataset contains all relevant information for the plot. Coord
    names are metadata keys, ordered as x, y, group, split. The default
    reference is None, or if all references are named the first from the
    list.
    """
    coords = {  # order matters: x, y, group, split
        cfg["x_by"]: list(group_metadata(metas, cfg["x_by"]).keys()),
        cfg["y_by"]: list(group_metadata(metas, cfg["y_by"]).keys()),
        cfg["group_by"]: list(group_metadata(metas, cfg["group_by"]).keys()),
        cfg["split_by"]: list(group_metadata(metas, cfg["split_by"]).keys()),
    }
    shape = [len(coord) for coord in coords.values()]
    var_data = xr.DataArray(np.full(shape, np.nan), dims=list(coords.keys()))
    data = xr.Dataset({"var": var_data}, coords=coords)
    # loop over each cell (coord combination) and load data if existing
    for coord_tuple in itertools.product(*coords.values()):
        selection = dict(zip(coords.keys(), coord_tuple))
        data['var'].loc[selection] = open_file(metas, **selection)
        # data[coord_tuple] = (list(coords.keys(), value))
    if None in data.coords[cfg["split_by"]].values:
        cfg.update({"default_split": None})
    else:
        cfg.update({"default_split": data.coords[cfg["split_by"]].values[0]})
    log.debug("using %s as default split", cfg["default_split"])
    log.debug("Loaded Data:")
    log.debug(data)
    return data


def split_legend(cfg, grid, data):
    """Create legend for references, based on split coordinate in the dataset.

    Mpl handles axes positions in relative figure coordinates. To anchor the
    legend to the origin of the first graph (bottom left) with fixed size,
    without messing up the layout for changing figure sizes a few extra steps
    are required.
    TODO: maybe `mpl_toolkits.axes_grid1.axes_divider.AxesDivider` simplifies
    this a bit by using `append_axes`.
    """

    fig = grid[0].get_figure()
    fig.canvas.draw()  # set axes position in figure (dont call tight_layout())
    size = cfg["legend"].get("size", 0.5)  # rect width in physical size (inch)
    fig_w, fig_h = fig.get_size_inches()  # physical size of figure
    ax_w, ax_h = (size / fig_w, size / fig_h)  # legend size in figure coords
    gap_x, gap_y = (0.3 / fig_w, 0.3 / fig_h)  # margins to plot in fig coords
    # anchor legend on origin of first plot or colorbar
    anchor = grid[0].get_position().bounds  # relative figure coordinates
    if cfg["legend"].get("position", "right") == "right":
        cbar_x = grid.cbar_axes[0].get_position().bounds[0]
        gap_x *= 0.8  # compensate colorbar padding
        anchor = (cbar_x + gap_x, anchor[1] - gap_y - ax_h)
    else:
        anchor = (anchor[0] - gap_x - ax_w, anchor[1] - gap_y - ax_h)
    # create legend as empty imshow like axes in figure coordinates
    legend = fig.add_axes([anchor[0], anchor[1], ax_w, ax_h])
    legend.imshow(np.zeros((1, 1)))  # same axes properties as main plot
    legend.set_xticks([])
    legend.set_yticks([])
    axy = legend.twinx()  # add twins to allow axes labels on all sides
    axy.set_yticks([])
    axx = legend.twiny()
    axx.set_xticks([])

    labels = data.coords[cfg["split_by"]].values
    label_at = [  # order matches get_triangle_nodes (halves and quarters)
        legend.set_ylabel,  # left
        axy.set_ylabel,  # right
        legend.set_xlabel,  # bottom
        axx.set_xlabel,  # top
    ]
    for i, _ in enumerate(labels):
        nodes = get_triangle_nodes(i, len(labels))
        colors = ["#bbb", "#ccc", "#ddd", "#eee"]
        patch = patches.Polygon(
            nodes,
            closed=True,
            facecolor=colors[i],
            edgecolor="black",
            linewidth=0.5,
            fill=True,
        )
        legend.add_patch(patch)
        label_at[i](labels[i])
    print(labels)


def overlay_reference(ax, data, triangle):
    """Create triangular overlays for given data and axes."""
    # use same colors as in main plot
    cmap = ax.get_images()[0].get_cmap()
    norm = ax.get_images()[0].norm
    for i, j in itertools.product(*map(range, data.shape)):
        if np.isnan(data[i, j]):
            continue
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


def plot_group(cfg, ax, data, title=None):
    """Create matrix for one subplot in ax using plt.imshow()

    by default split None is used, if all splits are named the first is
    used. Other splits will be added by overlaying triangles.
    """
    split = data.sel({cfg["split_by"]: cfg["default_split"]})
    print(f"Plotting group {title}")
    print(split)
    plot_matrix(
        split.values.T,  # 2d numpy array
        split.coords[cfg["y_by"]].values,  # y_labels
        split.coords[cfg["x_by"]].values,  # x_labels
        ax,
        cfg["plot_kwargs"],
    )
    if title is not None:
        ax.set_title(title)
    ax.set(**cfg["axes_properties"])


def get_triangle_nodes(position, total_count=2):
    """Returns list of three tuples with relative x, y coordinates for nodes of
    triangle (-0.5 to +0.5) at given quarters (total_count>2) or halves
    (total_count==2).

    NOTE: Order matters. Ensure axis labels for the legend match when changing.
    """
    if total_count < 3:
        halves = [
            [(0.5, -0.5), (-0.5, -0.5), (-0.5, 0.5)],  # top left
            [(0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)],  # bottom right
        ]
        return halves[position]
    else:
        quarters = [
            [(-0.5, -0.5), (0, 0), (-0.5, 0.5)],  # left
            [(0.5, -0.5), (0, 0), (0.5, 0.5)],  # right
            [(-0.5, 0.5), (0, 0), (0.5, 0.5)],  # bottom
            [(-0.5, -0.5), (0, 0), (0.5, -0.5)],  # top
        ]
        return quarters[position]


def plot_overlays(cfg, grid, data):
    """Call overlay_reference for each split in data and each group in grid."""
    split_count = data.shape[3]
    group_count = data.shape[2]
    for i in range(group_count):
        if split_count < 2:
            log.debug("No additional splits for overlay.")
            break
        if split_count > 4:
            log.warning("Too many splits for overlay, only 3 will be plotted.")
        group_data = data.isel({cfg["group_by"]: i})
        group_data = group_data.dropna(cfg["x_by"], how="all")
        for sss in range(split_count):
            split = group_data.isel({cfg["split_by"]: sss})
            split_label = split.coords[cfg["split_by"]].values.item()
            if split_label == cfg["default_split"]:
                log.debug("Skipping default split for overlay.")
                continue
            nodes = get_triangle_nodes(sss, split_count)
            overlay_reference(grid[i], split.values.T, nodes)


def plot(cfg, data):
    """Creates figure with subplots for each group, sets same color range and
    overlays additional references based on the content of data
    (xr.DataArray)"""
    fig = plt.figure(1, cfg.get("figsize", (5.5, 3.5)))
    group_count = len(data.coords[cfg["group_by"]])
    grid = ImageGrid(
        fig,
        111,  # similar to subplot(111)
        cbar_mode="single",
        cbar_location="right",
        cbar_pad=0.1,
        cbar_size=0.2,
        nrows_ncols=(1, group_count),
        axes_pad=0.1,
    )
    # remap colorbar to 10 discrete steps
    cmap = mpl.cm.get_cmap(cfg.get("cmap", "RdYlBu_r"), 10)
    cfg["plot_kwargs"]["cmap"] = cmap
    for i in range(group_count):
        group = data.isel({cfg["group_by"]: i})
        group = group.dropna(cfg["x_by"], how="all")
        title = None
        if group_count > 1:
            title = group.coords[cfg["group_by"]].values.item()
        plot_group(cfg, grid[i], group, title=title)
    # use same colorrange and colorbar for all subplots:
    unify_limits(cfg, grid)
    # set cb of first image as single cb for the figure
    grid.cbar_axes[0].colorbar(grid[0].get_images()[0], **cfg["cbar_kwargs"])
    if data.shape[3] > 1:
        plot_overlays(cfg, grid, data)
    if cfg.get("plot_legend", True):
        split_legend(cfg, grid, data)
    basename = "performance_metrics"
    fname = get_plot_filename(basename, cfg)
    plt.savefig(fname, bbox_inches="tight")
    log.info("Figure saved:")
    log.info(fname)


def set_defaults(cfg):
    """Set default values for most important config parameters."""
    cfg.setdefault("x_by", "alias")
    cfg.setdefault("y_by", "variable_group")
    cfg.setdefault("group_by", "project")
    cfg.setdefault("split_by", "split")  # extra facet
    cfg.setdefault("cbar_kwargs", {})
    cfg.setdefault("axes_properties", {})
    cfg.setdefault("plot_kwargs", {})
    cfg.setdefault("figsize", (7.5, 3.5))
    cfg.setdefault("legend", {})
    cfg["plot_kwargs"].setdefault("cmap", "RdYlBu_r")
    cfg["plot_kwargs"].setdefault("vmin", 0)
    cfg["plot_kwargs"].setdefault("vmax", 1)


def main(cfg):
    """Run the diagnostic."""
    set_defaults(cfg)
    metas = list(cfg["input_data"].values())
    remove_reference(metas)
    add_split_none(metas)
    dataset = load_data(cfg, metas)
    plot(cfg, dataset["var"])


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
