"""Subantarctic Front position (russell18jgr figure 3b).

Python port of russell18jgr-fig3b.ncl.

Following Orsi et al. (1995), the Subantarctic Front is the poleward
location of the 4C (277.15 K) isotherm at the model level closest to
(and shallower than) 400 m depth.

- Takes the time average of thetao and makes sure it is in Kelvin.
- Extracts the level closest to and above 400 m.
- Restricts to 80S-35S and sorts by longitude.
- Extracts the 277.15 K isoline, discards stray segments and repeats
  it over -360..360 degrees longitude.
- Overlays the front of every dataset in one panel.
"""

import logging
import os
import sys

import iris
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    run_diagnostic,
)

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import russell_common as rc

logger = logging.getLogger(os.path.basename(__file__))

ISOTHERM = 277.15  # 4 degrees Celsius in Kelvin
TARGET_DEPTH = 400.0
LAT_RANGE = (-80.0, -35.0)


def front_position(cube):
    """Return (x, y) of the front for one dataset."""
    cube = rc.time_mean(cube)
    data = rc.to_kelvin(cube.data)

    depth = cube.coord(axis="Z").points
    idx = rc.closest_index(TARGET_DEPTH, depth)
    if depth[idx] >= TARGET_DEPTH:
        idx = max(idx - 1, 0)
    field = data[idx]

    lat = rc.coord_1d_or_2d(cube, "latitude")
    lon = rc.coord_1d_or_2d(cube, "longitude")
    if lat.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon, lat)
    else:
        lat2d, lon2d = lat, lon

    rows = (lat2d[:, 0] >= LAT_RANGE[0]) & (lat2d[:, 0] <= LAT_RANGE[1])
    order = np.argsort(lon2d[0, :])
    field = field[rows][:, order]
    lat2d = lat2d[rows][:, order]
    lon2d = lon2d[rows][:, order]

    x, y = rc.extract_isoline(lon2d, lat2d, field, ISOTHERM)
    return rc.replicate_isoline(x, y)


def main(cfg):
    """Run the diagnostic."""
    input_data = [
        m for m in cfg["input_data"].values() if m["short_name"] == "thetao"
    ]
    years = rc.year_range_str(input_data)

    fig, axes = plt.subplots(figsize=(10, 8))
    ancestors, nc_files = [], []

    for meta in input_data:
        cube = rc.load_cube(meta["filename"], meta["short_name"])
        try:
            x, y = front_position(cube)
        except ValueError as exc:
            logger.warning("Skipping %s: %s", meta["dataset"], exc)
            continue
        style = rc.style_for(meta["dataset"], cfg.get("styleset", "CMIP5"))
        axes.plot(
            x,
            y,
            label=meta["dataset"],
            color=style["color"],
            linestyle=style["dash"],
            linewidth=style["thick"],
        )

        out = iris.cube.Cube(
            np.vstack([x, y]),
            var_name="position_of_subantarctic_front",
            long_name=(
                "row 0: longitude and row 1: latitude of the 277.15K isoline"
            ),
        )
        nc_name = get_diagnostic_filename(
            f"russell18jgr_fig3b_subantarctic-front-position_"
            f"{meta['dataset']}_{meta['start_year']}-"
            f"{meta['end_year']}",
            cfg,
        )
        iris.save(out, nc_name)
        nc_files.append(nc_name)
        ancestors.append(meta["filename"])

    axes.set_xlim(-60, 300)
    axes.set_ylim(-70, -40)
    axes.set_xticks(np.arange(-60, 301, 30))
    axes.set_yticks(np.arange(-70, -39, 2))
    axes.set_title(" Subantarctic Fronts", fontsize=12)
    axes.text(
        0.0,
        1.03,
        "Russell et al -2018 - Figure 3 b",
        transform=axes.transAxes,
        fontsize=10,
    )
    axes.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=7,
        frameon=False,
    )
    fig.tight_layout()

    plot_file = get_plot_filename(
        f"Russell18jgr_fig3_Subantarctic-Fronts_{years}", cfg
    )
    fig.savefig(plot_file, bbox_inches="tight", dpi=200)
    plt.close(fig)
    logger.info("Wrote %s", plot_file)

    record = rc.provenance_record("Russell et al 2018 figure 3b", ancestors)
    with ProvenanceLogger(cfg) as prov:
        for filename in nc_files + [plot_file]:
            prov.log(filename, record)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
