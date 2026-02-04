"""Diagnostic script that creates a 1D histogram for several datasets.

Description
-----------
This script computes 1D histograms for each dataset provided in the input
metadata. It then plots these histograms in a single figure for comparison.

Authors
-------
- Lisa Bock (DLR)

Cofiguration parameters through recipe
---------------------------------------
- bins: Number of bins for the histogram (default: 50).
- log_y: Boolean to set y-axis to logarithmic scale (default: False).
- ylim: List of two floats to set y-axis limits (default: None).
- suptitle: String for the figure's super title (default: None).
- bin_range: List of two floats specifying the min and max values for the histogram bins.
      Otherwise the min and max of the first dataset are used (default: None).

"""

import logging
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import dask.array as da
import iris.coords
import iris.cube
import matplotlib.pyplot as plt
import numpy as np
from iris.cube import CubeList

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
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
            "bock_lisa",
        ],
        "ancestors": ancestor_files,
    }
    return record


def set_defaults(cfg):
    """Set default values for most important config parameters."""
    cfg.setdefault("bins", 50)
    cfg.setdefault("bin_range", None)
    cfg.setdefault("log_y", False)
    cfg.setdefault("ylim", None)
    cfg.setdefault("suptitle", None)


def _compute_histograms(
    cfg: dict[str, Any],
    cube: iris.cube.Cube,
    bins: int,
) -> tuple[
    np.ndarray,
    np.ndarray,
]:
    """Compute the 1D histogram data."""
    array = cube.lazy_data().flatten()
    # Could use preprocessor function to ensure same masks
    select = ~(da.ma.getmaskarray(array))
    data = da.compress(select, da.ma.getdata(array))

    # Ensure min and max are set
    if cfg["bin_range"] is None:
        cfg["bin_range"] = [
            float(data.min().compute()),
            float(data.max().compute()),
        ]

    # Compute the data for the plots
    hist, edges = da.histogram(
        data, bins=bins, range=(cfg["bin_range"][0], cfg["bin_range"][1])
    )
    hist, edges = da.compute(hist, edges)

    # Calculate bin widths
    bin_widths = np.diff(edges)
    # Scale the data to get density (area under histogram = 1)
    scale = hist.sum() * bin_widths.sum()  # Total area of histogram
    hist = hist.astype(float) / scale
    hist = da.ma.masked_less_equal(hist, 0)
    hist = da.compute(hist)[0]

    return (hist, edges)


def plot_hist(
    cfg: dict[str, Any],
    datasets: Sequence[tuple[np.ndarray, np.ndarray]],
    label_x: str,
    labels: Sequence[str] | None = None,
    colors: Sequence[str] | None = None,
    log_y: bool = False,
    **kwargs,
) -> None:
    """Plot 1D histograms from the given datasets."""
    # normalize shapes and extract arrays
    hlist = []
    edges_ref = None
    for hist, edges in datasets:
        hlist.append(np.asarray(hist))
        if edges_ref is None:
            edges_ref = np.asarray(edges)
        else:
            if not np.allclose(edges_ref, np.asarray(edges)):
                raise ValueError("All datasets must share identical edges")

    centers = 0.5 * (edges_ref[:-1] + edges_ref[1:])
    widths = np.diff(edges_ref)

    fig, ax = plt.subplots()
    for i, hist in enumerate(hlist):
        h = np.where(hist > 0, hist, 0.0) if log_y else hist
        ax.bar(
            centers,
            h,
            width=widths,
            alpha=0.5,
            color=(None if colors is None else colors[i]),
            label=(None if labels is None else labels[i]),
        )

    ax.set_xlabel(label_x)
    ax.set_xlim(edges_ref[0], edges_ref[-1])
    ax.set_ylabel("Density")
    if labels is not None:
        ax.legend()
    if cfg["log_y"]:
        ax.set_yscale("log")
    if cfg["ylim"] is not None:
        ax.set_ylim(cfg["ylim"])
    if cfg["suptitle"] is not None:
        fig.suptitle(cfg["suptitle"])
    plt.tight_layout()


def data_to_cube(
    data: Sequence[tuple[np.ndarray, np.ndarray]],
    datasets: Sequence[str],
    group_name: str,
) -> iris.cube.Cube:
    """Store the histogram data in one iris cube."""
    cube_list = CubeList()
    for i, (hist, edges) in enumerate(data):
        if i == 0:
            hist_coord = iris.coords.DimCoord(
                var_name="bin",
                points=0.5 * (edges[:-1] + edges[1:]),
                bounds=np.array([edges[:-1], edges[1:]]).T,
            )
        cube_list.append(
            iris.cube.Cube(
                var_name=f"histogram_{group_name}_{datasets[i]}",
                data=hist,
                dim_coords_and_dims=[(hist_coord, 0)],
            )
        )

    return cube_list.merge()


def main(cfg: dict[str, Any]) -> None:
    """Create a 2D histogram with 1D histograms in the margins."""
    set_defaults(cfg)

    input_data = cfg["input_data"].values()

    groups = group_metadata(input_data, "variable_group", sort="dataset")

    all_data = []
    all_datasets = []
    all_filenames = []

    for group_name in groups:
        logger.info("Processing variable %s", group_name)
        for attributes in groups[group_name]:
            logger.info("Processing dataset %s", attributes["dataset"])
            input_file = attributes["filename"]

            cube = iris.load_cube(input_file)

            # Compute the histograms
            data = _compute_histograms(cfg, cube, cfg["bins"])

            all_data.append(data)
            all_datasets.append(attributes["dataset"])
            all_filenames.append(input_file)

        # Create the figure
        plot_hist(
            cfg,
            all_data,
            label_x=f"{cube.var_name} [{cube.units}]",
            labels=all_datasets,
            log_y=True,
        )

        # Save results
        basename = cfg["plot_filename"] + f"_{group_name}"
        caption = f"Scatterplot of {cube.long_name}"
        provenance = get_provenance_record(caption, all_filenames)
        save_figure(basename, provenance, cfg, dpi=cfg.get("dpi", 300))
        cube_hist = data_to_cube(all_data, all_datasets, group_name)
        save_data(f"{basename}", provenance, cfg, cube_hist)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
