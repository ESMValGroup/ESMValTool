"""Diagnostic script that creates a 2D histogram with 1D histograms in the margins.

It reproduces the same figure as
```
    seaborn.jointplot(
       kind="hist",
       marginal_kws={
          "stat": "probability",
       },
       stat="probability",
    )
```
but is much faster and has a lower memory footprint.
"""

import logging
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any, Literal

import dask.array as da
import iris.coords
import iris.cube
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn

from esmvaltool.diag_scripts.shared import (
    run_diagnostic,
    save_data,
    save_figure,
    select_metadata,
)

logger = logging.getLogger(Path(__file__).stem)


def get_provenance_record(
    caption: str,
    ancestor_files: Sequence[str],
) -> dict[str, str | Sequence[str]]:
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        "caption": caption,
        "plot_types": ["histogram"],
        "authors": [
            "lauer_axel",
            "bock_lisa",
            "schlund_manuel",
            "andela_bouwe",
        ],
        "ancestors": ancestor_files,
    }
    return record


def _compute_histograms(
    cube_x: iris.cube.Cube,
    cube_y: iris.cube.Cube,
    bins: int,
) -> tuple[
    np.ndarray,
    tuple[np.ndarray, np.ndarray],
    tuple[np.ndarray, np.ndarray],
]:
    """Compute the 2D and 1D histogram data."""
    array_x = cube_x.lazy_data().flatten()
    array_y = cube_y.lazy_data().flatten()
    # Could use preprocessor function to ensure same masks
    select = ~(da.ma.getmaskarray(array_x) | da.ma.getmaskarray(array_y))
    data_x = da.compress(select, da.ma.getdata(array_x))
    data_y = da.compress(select, da.ma.getdata(array_y))

    min_x, max_x, min_y, max_y = da.compute(
        data_x.min(), data_x.max(), data_y.min(), data_y.max()
    )

    # Compute the data for the plots
    hist_x, edges_x = da.histogram(data_x, bins=bins, range=(min_x, max_x))
    hist_y, edges_y = da.histogram(data_y, bins=bins, range=(min_y, max_y))
    hist = da.histogram2d(data_x, data_y, bins=[edges_x, edges_y])[0]
    # Scale the data according to seaborn.jointplot(stat="probability")
    scale = hist.sum()
    hist_x = hist_x.astype(float) / scale
    hist_y = hist_y.astype(float) / scale
    hist = hist.astype(float) / scale
    hist = da.ma.masked_less_equal(hist, 0)
    hist, hist_x, hist_y = da.compute(hist, hist_x, hist_y)

    return hist, (hist_x, edges_x), (hist_y, edges_y)


def _inv_transform(ax: mpl.axes.Axes, axis: Literal["x", "y"]) -> Callable:
    return getattr(ax, f"{axis}axis").get_transform().inverted().transform


def _plot_marginal(
    ax: mpl.axes.Axes,
    hist: np.ndarray,
    edges: np.ndarray,
    axis: Literal["x", "y"],
    **kwargs,
) -> None:
    """Plot the 1D histograms."""
    method = ax.barh if axis == "y" else ax.bar
    inv = _inv_transform(ax, axis)
    centers = inv(0.5 * (edges[:-1] + edges[1:]))
    widths = inv(np.diff(edges))
    method(centers, hist, widths, **kwargs)


def seaborn_histogram_jointplot(
    data: tuple[
        np.ndarray,
        tuple[np.ndarray, np.ndarray],
        tuple[np.ndarray, np.ndarray],
    ],
    label_x: str,
    label_y: str,
    suptitle: str | None = None,
    seaborn_settings: dict | None = None,
    joint_grid_kws: dict | None = None,
    joint_kws: dict | None = None,
    marginal_kws: dict | None = None,
    cbar: bool = True,
    cbar_kws: dict | None = None,
) -> None:
    """Plot the histogram data like `seaborn.jointplot` would plot it."""
    if seaborn_settings is not None:
        seaborn.set_theme(**seaborn_settings)

    joint_grid_kws = {} if joint_grid_kws is None else joint_grid_kws
    joint_kws = {} if joint_kws is None else joint_kws
    marginal_kws = {} if marginal_kws is None else marginal_kws
    cbar_kws = {} if cbar_kws is None else cbar_kws

    hist, (hist_x, edges_x), (hist_y, edges_y) = data

    x_width = np.diff(edges_x).mean()
    y_width = np.diff(edges_y).mean()
    xlim = (edges_x[0] - x_width, edges_x[-1] + x_width)
    ylim = (edges_y[0] - y_width, edges_y[-1] + y_width)

    grid = seaborn.JointGrid(xlim=xlim, ylim=ylim, **joint_grid_kws)
    grid.set_axis_labels(
        xlabel=label_x,
        ylabel=label_y,
    )
    if suptitle is not None:
        grid.figure.suptitle(suptitle)

    _plot_marginal(grid.ax_marg_y, hist_y, edges_y, axis="y", **marginal_kws)
    _plot_marginal(grid.ax_marg_x, hist_x, edges_x, axis="x", **marginal_kws)

    inv_edges_x = _inv_transform(grid.ax_joint, "x")(edges_x)
    inv_edges_y = _inv_transform(grid.ax_joint, "y")(edges_y)
    mesh = grid.ax_joint.pcolormesh(
        inv_edges_x, inv_edges_y, hist.T, **joint_kws
    )

    if cbar:
        grid.ax_joint.figure.colorbar(mesh, None, grid.ax_joint, **cbar_kws)
        # Reposition colorbar so it is to the right of marginals plot
        plt.subplots_adjust(left=0.13, right=0.83, top=0.9, bottom=0.1)
        # Get the current positions of the joint ax and the ax for the marginal x.
        pos_joint_ax = grid.ax_joint.get_position()
        pos_marg_x_ax = grid.ax_marg_x.get_position()
        # Reposition the joint ax so it has the same width as the marginal x ax.
        grid.ax_joint.set_position(
            (
                pos_joint_ax.x0,
                pos_joint_ax.y0,
                pos_marg_x_ax.width,
                pos_joint_ax.height,
            )
        )
        # Reposition the colorbar using new x positions and y positions of the joint ax.
        grid.figure.axes[-1].set_position(
            (0.85, pos_joint_ax.y0, 0.07, pos_joint_ax.height)
        )


def data_to_cubes(
    data: tuple[
        np.ndarray,
        tuple[np.ndarray, np.ndarray],
        tuple[np.ndarray, np.ndarray],
    ],
    cfg: dict[str, Any],
) -> tuple[
    iris.cube.Cube,
    iris.cube.Cube,
    iris.cube.Cube,
]:
    """Store the histogram data in iris cubes."""
    hist, (hist_x, edges_x), (hist_y, edges_y) = data
    hist_x_coord = iris.coords.DimCoord(
        var_name="bin_x",
        points=0.5 * (edges_x[:-1] + edges_x[1:]),
        bounds=np.array([edges_x[:-1], edges_x[1:]]).T,
    )
    hist_y_coord = iris.coords.DimCoord(
        var_name="bin_y",
        points=0.5 * (edges_y[:-1] + edges_y[1:]),
        bounds=np.array([edges_y[:-1], edges_y[1:]]).T,
    )
    cube_hist_x = iris.cube.Cube(
        var_name=f"{cfg['x']}_histogram",
        data=hist_x,
        dim_coords_and_dims=[(hist_x_coord, 0)],
    )
    cube_hist_y = iris.cube.Cube(
        var_name=f"{cfg['y']}_histogram",
        data=hist_y,
        dim_coords_and_dims=[(hist_y_coord, 0)],
    )
    cube_hist = iris.cube.Cube(
        var_name=f"{cfg['x']}_{cfg['y']}_histogram",
        data=hist,
        dim_coords_and_dims=[(hist_x_coord, 0), (hist_y_coord, 1)],
    )

    return cube_hist, cube_hist_x, cube_hist_y


def main(cfg: dict[str, Any]) -> None:
    """Create a 2D histogram with 1D histograms in the margins."""
    input_data = cfg["input_data"].values()

    filename_x = select_metadata(input_data, short_name=cfg["x"])[0][
        "filename"
    ]
    filename_y = select_metadata(input_data, short_name=cfg["y"])[0][
        "filename"
    ]
    cube_x = iris.load_cube(filename_x)
    cube_y = iris.load_cube(filename_y)

    # Compute the histograms
    data = _compute_histograms(cube_x, cube_y, cfg["bins"])

    # Create the figure
    seaborn_histogram_jointplot(
        data,
        label_x=f"{cube_x.var_name} [{cube_x.units}]",
        label_y=f"{cube_y.var_name} [{cube_y.units}]",
        suptitle=cfg["suptitle"],
        seaborn_settings=cfg.get("seaborn_settings"),
        joint_grid_kws=cfg.get("joint_grid_kws"),
        joint_kws=cfg.get("joint_kws"),
        cbar=cfg.get("cbar", True),
        cbar_kws=cfg.get("cbar_kws"),
        marginal_kws=cfg.get("marginal_kws"),
    )

    # Save results
    basename = cfg["plot_filename"]
    caption = f"Scatterplot of {cube_x.long_name} ({cfg['x']}) vs {cube_y.long_name} ({cfg['y']}) ({cfg['suptitle']})"
    provenance = get_provenance_record(caption, [filename_x, filename_y])
    save_figure(basename, provenance, cfg, dpi=cfg.get("dpi", 300))
    cube_hist, cube_hist_x, cube_hist_y = data_to_cubes(data, cfg)
    save_data(f"{basename}_{cube_hist.var_name}", provenance, cfg, cube_hist)
    provenance_x = get_provenance_record(caption, [filename_x])
    save_data(
        f"{basename}_{cube_hist_x.var_name}", provenance_x, cfg, cube_hist_x
    )
    provenance_y = get_provenance_record(caption, [filename_y])
    save_data(
        f"{basename}_{cube_hist_y.var_name}", provenance_y, cfg, cube_hist_y
    )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
