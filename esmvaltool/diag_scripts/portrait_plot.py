"""Portrait Plot Diagnostic.

Plot performance metrics of multiple datasets vs up to four references
A :doc:`documented example recipe </recipes/recipe_portrait>` to use this
diagnostic is provided as ``recipes/recipe_portrait_CMIP.yml``.

Description
-----------
This diagnostic provides plot functionalities for performance metrics,
and is written to be as flexible as possible to be adaptable to further use
cases. X and Y axis, grouping parameter and splits for each rectangle can be
configured in the recipe. All ``*_by`` parameters can be set to any metadata
key. To split by 'reference' this key needs to be set as extra_facet in recipe.

Authors
-------
- Lukas Lindenlaub (Universität Bremen, Germany)
- Diego Cammarano

Configuration parameters through recipe:
----------------------------------------
axes_properties: dict, optional
    Dictionary that gets passed to :meth:`matplotlib.axes.Axes.set`.
    Subplots can be widely customized.
    E.g. xlabel, ylabel, yticklabels, xmargin...
    By default {}.
cbar_kwargs: dict, optional
    Dictionary that gets passed to :meth:`matplotlib.pyplot.colorbar`.
    E.g. label, ticks...
    By default {}.
default_split: str, optional
    Data labeled with this string, will be used as main rectangles. All other
    splits will be plotted as overlays. This can be used to choose the base
    reference, while all references are labeled for the legend. If None, the
    first split will be used as default.
    By default None.
dpi: int, optional
    Dots per inch for the figure. By default 300.
domain: str, optional
    Domain for provenance. By default 'global'.
figsize: tuple of float, optional
    [width, height] of the figure in inches. The final figure will be saved
    with bbox_inches="tight", which can change the resulting aspect ratio.
    By default [7.5, 3.5].
group_by: str, optional
    Split portrait groups into multiple groups (one matrix per group).
    By default 'project'.
legend: dict, optional
    Customize, if, how and where the legend is plotted. The 'best' position
    and size of the legend depends on multiple parameters of the figure
    (i.e. lengths of labels, aspect ratio of the plots...). Might require
    manual adjustment of ``x``, ``y`` and ``size`` to fit the figure layout.
    Keys (each optional) that will be handled are:

    position: str or None, optional
        Position of the legend. Can be 'right' or 'left'.
        Or set to None to disable plotting the legend. By default 'right'.
    size: float, optional
        Size of the legend in Inches. By default 0.3.
    x_offset: float, optional
        Manually adjust horizontal position to save space or fix overlap.
        Number given in Inches. By default 0.
    y_offset: float, optional
        Manually adjust vertical position to save space or fix overlap.
        Number given in Inches. By default 0.

matplotlib_rc_params: dict, optional
    Optional :class:`matplotlib.RcParams` used to customize the portrait plot.
    Options given here will be passed to :func:`matplotlib.rc_context`.
nan_color: str or None, optional
    Matplotlib named color or hexcode for NaN values. If set to None,
    no triangles are plotted for NaN values.
    By default 'white'.
normalize: str or None, optional
    ('mean', 'median', 'centered_mean', 'centered_median', None).
    Divide by median or mean if not None. Subtract median/mean afterwards if
    centered.
    By default 'centered_median'.
plot_kwargs: dict, optional
    Dictionary that gets passed as kwargs to
    :meth:`matplotlib.axes.Axes.imshow`. Colormaps will be converted to 11
    discrete steps automatically.
    Default colormap is ``cmap='RdYlBu_r'`` with limits ``vmin=-0.5`` and
    ``vmax=0.5``.
plot_legend: bool, optional
    If True, a legend is plotted, when multiple splits are given.
    By default True.
split_by: str, optional
    The rectangles can be split into 2-4 triangles. This is used to show
    metrics for different references. For this case there is no need to change
    this parameter. Multiple variables can be set in the recipe with ``split``
    assigned as extra_facet to label the different references. Data without
    a split assigned will be plotted as main rectangles, this can be changed
    by setting default_split parameter.
    By default 'split'.
x_by: str, optional
    Metadata key for x coordinate.
    By default 'alias'.
y_by: str, optional
    Metadata key for y coordinate.
    By default 'variable_group'.

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
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
    select_metadata,
)

log = logging.getLogger(__name__)


def get_provenance(cfg):
    """Return provenance for this diagnostic."""
    return {
        "ancestors": list(cfg["input_data"].keys()),
        "authors": ["lindenlaub_lukas", "cammarano_diego"],
        "caption": "RMSE performance metric",
        "domains": [cfg["domain"]],
        "plot_types": ["portrait"],
        "references": [
            "gleckler08jgr",
        ],
        "statistics": ["rmsd"],
    }


def unify_limits(grid):
    """Ensure same limits for all subplots."""
    vmin, vmax = np.inf, -np.inf
    images = [ax.get_images()[0] for ax in grid]
    for img in images:
        vmin = min(vmin, img.get_clim()[0])
        vmax = max(vmax, img.get_clim()[1])
    for img in images:
        img.set_clim(vmin, vmax)


def plot_matrix(data, row_labels, col_labels, axe, plot_kwargs):
    """Create an image for given data."""
    img = axe.imshow(data, **plot_kwargs)
    # Show all ticks and label them with the respective list entries.
    axe.set_xticks(np.arange(data.shape[1]), labels=col_labels)
    axe.set_yticks(np.arange(data.shape[0]), labels=row_labels)
    # Rotate the tick labels and set their alignment.
    plt.setp(
        axe.get_xticklabels(),
        rotation=90,
        ha="right",
        va="center",
        rotation_mode="anchor",
    )
    axe.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor=True)
    axe.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor=True)
    axe.grid(which="minor", color="black", linestyle="-", linewidth=0.8)
    axe.tick_params(which="both", bottom=False, left=False)
    return img


def remove_reference(metas):
    """Remove reference for metric from list of metadata."""
    for meta in list(metas):  # list() creates a copy to allow remove in place
        if meta.get("reference_for_metric", False):
            metas.remove(meta)


def add_missing_facets(cfg, metas):
    """Ensure that all facets are present in metadata."""
    for meta in metas:
        facet_config = ["x_by", "y_by", "group_by", "split_by"]
        facets = [cfg[key] for key in facet_config]
        for facet in facets:
            meta.setdefault(facet, "unknown")


def open_file(metadata, **selection):
    """Try to find a single file for selection and return data.

    If multiple files are found, raise an error. If no file is found,
    return np.nan.
    """
    metas = select_metadata(metadata, **selection)
    if len(metas) > 1:
        raise ValueError(f"Multiple files found for {selection}")
    if len(metas) < 1:
        log.debug("No files found for %s", selection)
        return np.nan
    log.debug("File found for %s", selection)
    das = xr.open_dataset(metas[0]["filename"])
    varname = list(das.data_vars.keys())[0]
    try:
        return das[varname].values.item()
    except ValueError as exc:
        msg = f"Expected scalar in input file {metas[0]['filename']}."
        raise ValueError(msg) from exc


def load_data(cfg, metas):
    """Load all netcdf files from metadata into xarray dataset.

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
        selection = dict(zip(coords.keys(), coord_tuple, strict=True))
        data["var"].loc[selection] = open_file(metas, **selection)
    if cfg["default_split"] is None:
        cfg["default_split"] = data.coords[cfg["split_by"]].values[0]
    log.debug("Using %s as default split", cfg["default_split"])
    log.debug("Loaded Data: %s", data)
    return data


def split_legend(cfg, grid, data):
    """Create legend for references, based on split coordinate in the dataset.

    Mpl handles axes positions in relative figure coordinates. To anchor the
    legend to the origin of the first graph (bottom left) with fixed size,
    without messing up the layout for changing figure sizes, a few extra steps
    are required.
    NOTE: maybe ``mpl_toolkits.axes_grid1.axes_divider.AxesDivider`` simplifies
    this a bit by using ``append_axes``.
    """
    grid[0].get_figure().canvas.draw()  # set axes position in figure
    size = cfg["legend"]["size"]  # rect width in physical size (inch)
    fig_size = grid[0].get_figure().get_size_inches()  # physical figure size
    ax_size = (size / fig_size[0], size / fig_size[1])  # legend (fig coords)
    gaps = [0.3 / fig_size[0], 0.3 / fig_size[1]]  # margins (fig coords)
    # anchor legend on origin of first plot or colorbar
    anchor = grid[0].get_position().bounds  # relative figure coordinates
    if cfg["legend"]["position"] == "right":
        cbar_x = grid.cbar_axes[0].get_position().bounds[0]
        gaps[0] *= 0.8  # compensate colorbar padding
        anchor = (
            cbar_x + gaps[0] + cfg["legend"]["x_offset"],
            anchor[1] - gaps[1] - ax_size[1] + cfg["legend"]["y_offset"],
        )
    else:
        anchor = (
            anchor[0] - gaps[0] - ax_size[0] + cfg["legend"]["x_offset"],
            anchor[1] - gaps[1] - ax_size[1] + cfg["legend"]["y_offset"],
        )
    # create legend as empty imshow like axes in figure coordinates
    axes = {"main": grid[0].get_figure().add_axes([*anchor, *ax_size])}
    axes["main"].imshow(np.zeros((1, 1)))  # same axes properties as main plot
    axes["main"].set_xticks([])
    axes["main"].set_yticks([])
    axes["twiny"], axes["twinx"] = [axes["main"].twiny(), axes["main"].twinx()]
    axes["twinx"].set_yticks([])
    axes["twiny"].set_xticks([])
    label_at = [  # order matches get_triangle_nodes (halves and quarters)
        axes["main"].set_ylabel,  # left
        axes["twinx"].set_ylabel,  # right
        axes["main"].set_xlabel,  # bottom
        axes["twiny"].set_xlabel,  # top
    ]
    for i, label in enumerate(data.coords[cfg["split_by"]].values):
        nodes = get_triangle_nodes(i, len(data.coords[cfg["split_by"]].values))
        axes["main"].add_patch(
            patches.Polygon(
                nodes,
                closed=True,
                facecolor=["#bbb", "#ccc", "#ddd", "#eee"][i],
                edgecolor="black",
                linewidth=0.5,
                fill=True,
            ),
        )
        label_at[i](label)


def overlay_reference(cfg, axe, data, triangle):
    """Create triangular overlays for given data and axes."""
    # use same colors as in main plot
    cmap = axe.get_images()[0].get_cmap()
    norm = axe.get_images()[0].norm
    if cfg["nan_color"] is not None:
        cmap.set_bad(cfg["nan_color"])
    for i, j in itertools.product(*map(range, data.shape)):
        if np.isnan(data[i, j]) and cfg["nan_color"] is None:
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
        axe.add_patch(patch)


def plot_group(cfg, axe, data, title=None):
    """Create matrix for one subplot in ax using plt.imshow.

    By default split None is used, if all splits are named the first is
    used. Other splits will be added by overlaying triangles.
    """
    split = data.sel({cfg["split_by"]: cfg["default_split"]})
    plot_matrix(
        split.values.T,  # 2d numpy array
        split.coords[cfg["y_by"]].values,  # y_labels
        split.coords[cfg["x_by"]].values,  # x_labels
        axe,
        cfg["plot_kwargs"],
    )
    if title is not None:
        axe.set_title(title)
    axe.set(**cfg["axes_properties"])


def get_triangle_nodes(position, total_count=2):
    """Return list of nodes with relative x, y coordinates.

    The nodes of the triangle are given as list of three tuples. Each tuple
    contains relative coordinates (-0.5 to +0.5). For total of <= 2 a top left
    (position=0) and bottom right (position=1) rectangle is returned.
    For higher counts (3 or 4) one quartile is returned for each position.
    NOTE: Order matters. Ensure axis labels for the legend match when changing.
    """
    if total_count < 3:
        halves = [
            [(0.5, -0.5), (-0.5, -0.5), (-0.5, 0.5)],  # top left
            [(0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)],  # bottom right
        ]
        return halves[position]
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
            overlay_reference(cfg, grid[i], split.values.T, nodes)


def plot(cfg, data):
    """Create figure with subplots for each group and save to NetCDF.

    Sets same color range and overlays additional references based on
    the content of data (xr.DataArray).
    """
    fig = plt.figure(1, cfg["figsize"])
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
    cmap = mpl.cm.get_cmap(cfg["plot_kwargs"]["cmap"], 10)
    cfg["plot_kwargs"]["cmap"] = cmap
    for i in range(group_count):
        group = data.isel({cfg["group_by"]: i})
        group = group.dropna(cfg["x_by"], how="all")
        title = None
        if group_count > 1:
            title = group.coords[cfg["group_by"]].values.item()
        plot_group(cfg, grid[i], group, title=title)
    # use same colorrange and colorbar for all subplots:
    unify_limits(grid)
    # set cb of first image as single cb for the figure
    grid.cbar_axes[0].colorbar(grid[0].get_images()[0], **cfg["cbar_kwargs"])
    if data.shape[3] > 1:
        plot_overlays(cfg, grid, data)
    if cfg["plot_legend"] and data.shape[3] > 1:
        split_legend(cfg, grid, data)
    basename = "portrait_plot"
    fname = get_plot_filename(basename, cfg)
    plt.savefig(fname, bbox_inches="tight", dpi=cfg["dpi"])
    with ProvenanceLogger(cfg) as prov_logger:
        prov_logger.log(fname, get_provenance(cfg))
    log.info("Figure saved:")
    log.info(fname)


def normalize(array, method, dims):
    """Divide and shift values along dims depending on method."""
    shift = 0
    norm = 1
    if "mean" in method:
        norm = array.mean(dim=dims)
    elif "median" in method:
        norm = array.median(dim=dims)
    if "centered" in method:
        shift = norm
    normalized = (array - shift) / norm
    return normalized


def set_defaults(cfg):
    """Set default values for most important config parameters."""
    cfg.setdefault("axes_properties", {})
    cfg.setdefault("cbar_kwargs", {})
    cfg.setdefault("default_split", None)
    cfg.setdefault("dpi", 300)
    cfg.setdefault("domain", "global")
    cfg.setdefault("figsize", (7.5, 3.5))
    cfg.setdefault("group_by", "project")
    cfg.setdefault("legend", {})
    cfg["legend"].setdefault("position", "right")
    cfg["legend"].setdefault("size", 0.3)
    cfg["legend"].setdefault("x_offset", 0)
    cfg["legend"].setdefault("y_offset", 0)
    cfg.setdefault("matplotlib_rc_params", {})
    cfg.setdefault("nan_color", "white")
    cfg.setdefault("normalize", "centered_median")
    cfg.setdefault("plot_kwargs", {})
    cfg["plot_kwargs"].setdefault("cmap", "RdYlBu_r")
    cfg["plot_kwargs"].setdefault("vmin", -0.5)
    cfg["plot_kwargs"].setdefault("vmax", 0.5)
    cfg.setdefault("plot_legend", True)
    cfg.setdefault("split_by", "split")  # extra facet
    cfg.setdefault("x_by", "alias")
    cfg.setdefault("y_by", "variable_group")


def sort_data(cfg, dataset):
    """Sort the dataset along by custom or alphabetical order."""
    dataset = dataset.sortby(
        [
            dataset[cfg["x_by"]].str.lower(),
            dataset[cfg["y_by"]].str.lower(),
            dataset[cfg["group_by"]].str.lower(),
            dataset[cfg["split_by"]].str.lower(),
        ],
    )
    if cfg["x_by"] in ["alias", "dataset"]:
        # NOTE: not clean, but it works for many cases
        mm_stats = [
            v
            for v in dataset[cfg["x_by"]].values
            if "Mean" in v or "Median" in v or "Percentile" in v
        ]
        others = [
            v
            for v in dataset[cfg["x_by"]].values
            if "Mean" not in v and "Median" not in v and "Percentile" not in v
        ]
        new_order = mm_stats + others
        dataset = dataset.reindex({cfg["x_by"]: new_order})
    return dataset


def save_to_netcdf(cfg, data):
    """Save the final dataset to a NetCDF file."""
    basename = "portrait"
    fname = get_diagnostic_filename(basename, cfg, extension="nc")
    data.to_netcdf(fname)
    log.info("NetCDF file saved:")
    log.info(fname)
    with ProvenanceLogger(cfg) as prov_logger:
        prov_logger.log(fname, get_provenance(cfg))


def main(cfg):
    """Run the diagnostic."""
    set_defaults(cfg)
    metas = list(cfg["input_data"].values())
    remove_reference(metas)
    add_missing_facets(cfg, metas)
    dataset = load_data(cfg, metas)
    dataset = sort_data(cfg, dataset)
    if cfg["normalize"] is not None:
        dataset["var"] = normalize(
            dataset["var"],
            cfg["normalize"],
            [cfg["x_by"], cfg["group_by"]],
        )
    with mpl.rc_context(cfg["matplotlib_rc_params"]):
        plot(cfg, dataset["var"])
    save_to_netcdf(cfg, dataset["var"])


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
