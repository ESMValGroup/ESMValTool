"""Diagnostic script that creates a 1D histogram plot for several datasets.

Description
-----------
This script plots 1D histograms for each dataset provided in the input
metadata. It uses precomputed histogram data from the preprocessor "histogram".

Authors
-------
- Lisa Bock (DLR)

Cofiguration parameters through recipe
---------------------------------------
- log_y: Boolean to set y-axis to logarithmic scale (default: False).
- ylim: List of two floats to set y-axis limits (default: None).
- suptitle: String for the figure's super title (default: None).

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
    cfg.setdefault("log_y", False)
    cfg.setdefault("ylim", None)
    cfg.setdefault("suptitle", None)
    cfg.setdefault("plot_filename", "histogram_diag")


def plot_hist(
    cfg: dict[str, Any],
    datasets: Sequence[tuple[np.ndarray, np.ndarray]],
    xlabel: str,
    labels: Sequence[str] | None = None,
    colors: Sequence[str] | None = None,
) -> None:
    """Plot 1D histograms from the given datasets."""
    # normalize shapes and extract arrays
    hlist = []
    centers_ref = None
    for hist, centers in datasets:
        hlist.append(np.asarray(hist))
        if centers_ref is None:
            centers_ref = np.asarray(centers)

    widths = np.diff(centers_ref, prepend=centers_ref[0])

    fig, ax = plt.subplots()
    for i, hist in enumerate(hlist):
        h = np.where(hist > 0, hist, 0.0) if cfg["log_y"] else hist
        ax.bar(
            centers_ref,
            h,
            width=widths,
            alpha=0.5,
            color=(None if colors is None else colors[i]),
            label=(None if labels is None else labels[i]),
        )

    ax.set_xlabel(xlabel)
    ax.set_xlim(
        centers_ref[0] - 0.5 * widths[0], centers_ref[-1] + 0.5 * widths[-1],
    )
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
    for i, (hist, centers) in enumerate(data):
        if i == 0:
            hist_coord = iris.coords.DimCoord(
                var_name="bin",
                points=centers,
            )
        cube_list.append(
            iris.cube.Cube(
                var_name=f"histogram_{group_name}_{datasets[i]}",
                data=hist,
                dim_coords_and_dims=[(hist_coord, 0)],
            ),
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
            hist = da.compute(cube.data)[0]
            centers = cube.coord("lwe_precipitation_rate").points
            data = (hist, centers)

            all_data.append(data)
            all_datasets.append(attributes["dataset"])
            all_filenames.append(input_file)

        # Create the figure
        if cfg["xlabel"] is not None:
            xlabel = cfg["xlabel"]
        else:
            xlabel = f"{cube.var_name} [{cube.units}]"
        plot_hist(cfg, all_data, xlabel=xlabel, labels=all_datasets)

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
