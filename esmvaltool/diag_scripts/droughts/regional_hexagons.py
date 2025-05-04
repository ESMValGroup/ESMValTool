#!/usr/bin/env  python
"""Hexagonal overview plot for IPCC AR6 WG1 reference regions.

Configuration options in recipe
-------------------------------
group_by: str, optional (default: "dataset")
    The `group_by` config parameter can be used to create individual figures
    for input data which differs in values of the given key.
split_by: str, optional (default: "exp")
    Metadata key to split the data into different tiles of the hexagon.
    This is ignored for `split_by_statistic: True`.
    Only keys with 6 or less different values are supported
    (1,2,3,6 can be distributed symmetrically).
split_by_statistic: bool, optional (default: False)
    Split the hexagons into different tiles for each statistic,
    rather than different metadata.
statistics: list, optional (default: ['mean'])
    Any of the operators valid for `esmvalcore.preprocessor.area_statistics`
    can be used: 'mean', 'median', 'min', 'max', 'std_dev', 'sum', 'variance'
    or 'rms'. The regions are collapsed using the given operator.
    If not `split_by_statistic: True`, a figure is created for each operator.
plot_mmm: bool, optional (default: True)
    wether to plot multi-model mean
    TODO: support this and plot_models
cmap: string, optional (default: 'YlOrRd')
    colormap to use
vmin: float, optional (default: None)
    minimum value for colormap
vmax: float, optional (default: None)
    maximum value for colormap
cbar: bool, optional (default: True)
    wether to plot colorbar
labels: dict, optional
    dictionary with labels for each split. Dict keys must match the values
    corresponding to the split_by key. Set it to False to disable legend.
    Default: generated from statistics or split.
    TODO: implement dicts (rn only list works)
strip_plot: bool, optional (default: False)
    wether to plot the colorbar to a seperate file
cb_label: string, optional (default: 'Decadal change of SPEI')
    colorbar label
select_metadata: dict, optional
    limit the metadata to use for the plot. Keys must match metadata keys.
    Default: {'short_name': 'spei', 'diffmap_metric': 'diff'}
shapefile: string, optional (default: None)
    For IPCC WGI reference regions use
    `ar6_regions/IPCC-WGI-reference-regions-v4.shp`.
    TODO: move to preprocessor
exclude_regions: list, optional (default: [])
    regions that are excluded from the plot.
    TODO: move to preprocessor as well?
filename: string, optional (default: {group}_{split}_{operator}.png)
    Filename template for plot files. `group`, `split` and `operator`
    and all metadata keys can be used as placeholders.
    TODO: Replace suffix by filename (Not implemented yet).
show_values: bool, optional (default: False)
    Add corresponding value in a new line to each regions label. For multiple
    splits only the value of the first split is shown.
"""

import logging

import iris
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm as mplcm
from matplotlib import colors as mplcolors
from matplotlib.patches import Polygon

import esmvaltool.diag_scripts.droughtindex.utils as ut
import esmvaltool.diag_scripts.shared as e
from esmvaltool.diag_scripts.shared import (
    # ProvenanceLogger,
    group_metadata,
    select_metadata,
)

logger = logging.getLogger(__file__)


def plot_colorbar(
    cfg: dict,
    plotfile: str,
    plot_kwargs: dict,
    orientation="vertical",
    mappable=None,
) -> None:
    """Creates a colorbar in its own figure. Usefull for multi panel plots."""
    fig, ax = plt.subplots(figsize=(1, 4), layout="constrained")
    if mappable is None:
        cmap = plot_kwargs.get("cmap", "RdYlBu")
        norm = mplcolors.Normalize(
            vmin=plot_kwargs.get("vmin"),
            vmax=plot_kwargs.get("vmax"),
        )
        mappable = mplcm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(
        mappable,
        cax=ax,
        orientation=orientation,
        label=plot_kwargs["cbar_label"],
    )
    plotfile = plotfile.removesuffix(".png")
    fig.savefig(plotfile + "_cb.png", bbox_inches="tight")


def hexmap(
    cfg,
    regions,
    values,
    labels=None,
    suffix="",
    r=0.8,
    bg=False,
    filename=None,
    draw_nans=True,
):
    """Plot hexagons for IPPC WG1 reference regions for data pairs

    Parameters
    ----------
    cfg
        config with indexname and optional plot parameters:
        cmap: string, vmin: float, vmax: float, cbar: bool,
    regions
        list of IDs (strings) for the reference regions
    values
        list (each split) of list (each region) of values (floats)
    labels, optional
        names of the splits to create legend. Same length as values
    texts, optional
        array (same length as regions) of strings which are added to each cell,
        by default None
    suffix, optional
        string to include in filename. TODO: replace by filename format str.
    r, optional
        scaling parameter of a hexagon, by default 1
    bg, optional
        draw white background instead of transparent, by default False
    draw_nans, optional
        draw gray polygons for nan values, by default True.
    filename, optional
        filename to save the plot. If not provided the filename is generated
        automatically using shared.get_plot_filename(). By default None
        TODO: format with metadata.

    Raises
    ------
    ValueError
        For missmatching input array lengths
    """
    if not len(regions) == len(values[0]):
        raise ValueError("regions and values must have the same length")
    if labels and not len(labels) == len(values):
        raise ValueError("values and labels must have the same length")
    values = np.array(values)  # np array makes it easier to deal with inf/nan
    figsize = (12, 6) if cfg["cbar"] and not cfg["strip_plot"] else (10, 6)
    fig, axx = plt.subplots(figsize=figsize, dpi=300, frameon=False)
    axx.tick_params(
        bottom=False,
        left=False,
        labelbottom=False,
        labelleft=False,
    )
    axx.set_xlim(-0.5, 19.5)
    axx.set_ylim(0, 12)
    cmap = plt.get_cmap(cfg.get("cmap", "YlOrRd"))
    cmap.set_bad(cfg.get("cmap_nan", "lightgray"), 1.0)
    cmap.set_over(cfg.get("cmap_inf", "dimgray"), 1.0)
    if not bg:
        plt.axis("off")
    # calculate hexagon positions based on figure and scale
    rx = np.sqrt(3) / 2 * r  # 0.8660254037844386
    corners = np.array(
        [
            [0, r],
            [rx, r / 2],
            [rx, -r / 2],
            [0, -r],
            [-rx, -r / 2],
            [-rx, r / 2],
        ],
    )
    cells = ut.get_hex_positions()  # dict of coordinates for hexagons
    cells = {
        a: [c[0] * rx, (8 - c[1]) * (3 / 2 * r)] for a, c in cells.items()
    }
    vmin = cfg.get("vmin", values[np.isfinite(values)].min())
    vmax = cfg.get("vmax", values[np.isfinite(values)].max())
    norm = mplcolors.Normalize(vmin=vmin, vmax=vmax)
    print(vmin)
    print(vmax)
    if "levels" in cfg:
        lvls = cfg["levels"]
        # cmap_colors = [cmap(i/(vmax-vmin)) for i in lvls]
        cmap_colors = [cmap(norm(lvl)) for lvl in lvls]
        cmap = mplcolors.ListedColormap(cmap_colors)
        norm = mplcolors.BoundaryNorm(lvls, cmap.N)
        cmap.set_bad(cfg.get("cmap_nan", "lightgray"), 1.0)
        cmap.set_over(cfg.get("cmap_inf", "lightgray"), 1.0)
    abbrs = ut.get_region_abbrs()
    # draw hexagon for each region
    for iii, name in enumerate(regions):
        if name not in abbrs:
            logger.warning("No hex cells for %s. Skipping.", name)
            continue
        abbr = abbrs[name]
        exclude = cfg.get("exclude_regions", [])
        if abbr in exclude or name in exclude:
            logger.info("Region %s excluded. Skipping hexagon.", abbr)
            continue
        if abbr not in cells:
            logger.warning("Region %s not found. Skipping hexagon.", abbr)
            continue
        if len(cfg["regions"]) > 1 and abbr not in cfg["regions"]:
            logger.info("Region %s not in regions. Skipping hexagon.", abbr)
            # continue
            values[0][iii] = np.nan
        c = cells[abbr]  # center
        hex_corners = np.array([c, c, c, c, c, c]) + corners
        hc = list(hex_corners)
        hc.append(hc[0])  # last corner = first corner for closed polies
        # create polygon vertices depending on the number of values per hex
        if len(values) > 3:  # 6 pieces
            verts = [(hc[i], hc[i + 1], c) for i in range(6)]
        elif len(values) > 2:
            verts = [
                (hc[2 * i], hc[2 * i + 1], hc[2 * i + 2], c) for i in range(3)
            ]
        elif len(values) > 1:
            verts = np.array([hc[0:4], hc[3:]])
        else:
            verts = [hc]
        color = "black"
        for vvv in range(len(values[:])):
            val = values[vvv][iii]
            color = cmap(norm(val))
            if not draw_nans and np.isnan(val):
                continue
            poly = Polygon(verts[vvv], ec="white", fc=color, lw=1)
            axx.add_artist(poly)
        base = Polygon(hex_corners, ec="white", fc=None, lw=2, fill=False)
        axx.add_artist(base)
        tval = values[0][iii]  # first split for text label
        text = abbr
        if cfg["show_values"] and not np.isnan(tval):
            text = f"{abbr}\n{tval:.2f}"
        text_style = {"ha": "center", "va": "center", "fontsize": 10}
        axx.text(c[0], c[1], text, **text_style, color=ut.font_color(color))

    if filename is None:
        filename = ut.get_plot_filename(cfg, "hexmap_regions" + suffix)

    mappable = mplcm.ScalarMappable(norm=norm, cmap=cmap)
    if cfg.get("cbar", True) and not cfg.get("strip_plot", False):
        plt.colorbar(mappable, ax=axx)
    if cfg.get("strip_plot", False):
        cb_kwargs = {"cbar_label": cfg.get("cb_label", "")}
        plot_colorbar(cfg, filename, cb_kwargs, mappable=mappable)
    if labels:
        pos = (1, 2)
        c = pos
        anchors = [  # corner, ha, va
            [(rx, rx), "bottom", "left"],
            [(1.5 * rx, 0), "center", "left"],
            [(rx, -rx), "top", "left"],
            [(-rx, -rx), "top", "right"],
            [(-1.5 * rx, 0), "center", "right"],
            [(-rx, rx), "bottom", "right"],
        ]
        # TODO: repeated below.. make function
        hex_corners = np.array([c, c, c, c, c, c]) + corners
        hc = list(hex_corners)
        hc.append(hc[0])  # last corner = first corner for closed polies
        # create polygon vertices depending on the number of values per hex
        if len(values) > 3:  # 6 pieces
            verts = [(hc[i], hc[i + 1], c) for i in range(6)]
        elif len(values) > 2:
            verts = [
                (hc[2 * i], hc[2 * i + 1], hc[2 * i + 2], c) for i in range(3)
            ]
        elif len(values) > 1:
            verts = np.array([hc[0:4], hc[3:]])
        else:
            verts = [hc]
        for vvv in range(len(values[:])):
            color = ["#fff", "#aaa", "#777", "#555", "#333", "#000"][vvv]
            poly = Polygon(verts[vvv], ec="white", fc=color, lw=1)
            axx.add_artist(poly)
            anc = anchors[vvv]
            lpos = (anc[0][0] + pos[0], anc[0][1] + pos[1])
            lab = labels[vvv]
            axx.text(lpos[0], lpos[1], lab, fontsize=10, va=anc[1], ha=anc[2])
    fig.savefig(filename, bbox_inches="tight")


def ensure_single_meta(meta, txt):
    """Raise error if there is not exactly one entry in meta list."""
    if len(meta) == 0:
        raise ValueError(f"No files for {txt}.")
    if len(meta) > 1:
        raise ValueError(f"Too many files for {txt}.")
    return meta[0]


def load_and_plot_splits(cfg, splits, group, operator):
    """Create hexmap with tiles belonging to different meta data / files."""
    if len(splits) > 6:
        raise ValueError("Too many inputs for {group}. Max: 6.")
    labels = cfg.get("labels", list(splits.keys()))
    values, regions, meta = [], [], {}
    # extract regional average from each cube
    for split, meta in splits.items():
        meta = ensure_single_meta(meta, f"{group}/{split}")
        cube = iris.load_cube(meta["filename"])
        collapsed = ut.regional_stats(cfg, cube, operator)
        values.append(collapsed.data)
        regions = collapsed.coord("shape_id").points
    basename = cfg.get("basename", f"hexmap_regions_{group}")
    fmeta = meta.copy()
    fmeta["group"] = group
    filename = ut.get_plot_filename(cfg, basename, meta)
    hexmap(cfg, regions, values, labels=labels, filename=filename)


def load_and_plot_stats(cfg, metas, group, split, statistics):
    """Create hexmap with tiles for different operators for the same data."""
    meta = ensure_single_meta(metas, f"{group}")
    cube = iris.load_cube(meta["filename"])
    values, regions = [], []
    for operator in statistics:
        collapsed = ut.regional_stats(cfg, cube, operator)
        values.append(collapsed.data)
        regions = collapsed.coord("shape_id").points
    hexmap(
        cfg,
        regions,
        values,
        labels=statistics,
        suffix=f"_{group}_{split}_statistics",
    )


def set_defaults(cfg):
    """Ensure all config parameters are set."""
    cfg.setdefault("group_by", "dataset")
    cfg.setdefault("statistics", ["mean"])
    cfg.setdefault("cmap", "YlOrRd")
    cfg.setdefault("vmin", None)
    cfg.setdefault("vmax", None)
    cfg.setdefault("cbar", True)
    cfg.setdefault("labels", None)
    cfg.setdefault("show_values", False)
    cfg.setdefault("cb_label", "Decadal change of SPEI")
    cfg.setdefault("strip_plot", False)
    cfg.setdefault("split_by", "exp")
    cfg.setdefault("split_by_statistic", False)
    cfg.setdefault("plot_mmm", True)
    cfg.setdefault("exclude_regions", [])
    cfg.setdefault("regions", [])
    cfg.setdefault("filename", "{group}_{split}_{operator}.png")
    cfg.setdefault(
        "select_metadata",
        {"short_name": "spei", "diffmap_metric": "diff"},
    )


def main(cfg):
    set_defaults(cfg)
    # select metadata
    metas = cfg["input_data"].values()
    metas = select_metadata(metas, **cfg.get("select_metadata", {}))
    groups = group_metadata(metas, cfg["group_by"])
    statistics = cfg["statistics"]
    for group, metas in groups.items():
        # plot splits (segments in hex) for each group (figures)
        splits = group_metadata(metas, cfg["split_by"])
        if cfg.get("split_by_statistic"):
            for split, split_metas in splits.items():
                load_and_plot_stats(cfg, split_metas, group, split, statistics)
        else:
            for operator in statistics:
                load_and_plot_splits(cfg, splits, group, operator)


if __name__ == "__main__":
    with e.run_diagnostic() as config:
        main(config)
