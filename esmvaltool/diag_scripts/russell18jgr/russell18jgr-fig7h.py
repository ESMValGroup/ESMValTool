"""Zonal mean CO2 flux (russell18jgr figure 7h).

Python port of russell18jgr-fig7h.ncl.

- Expects time-averaged, land-masked fgco2 from the preprocessor.
- Converts kg m-2 s-1 to gC m-2 yr-1 (factor -3.1536e10; the sign flips
  the flux to the sea-to-air convention of the paper).
- Plots the zonal average against latitude for all datasets.
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

UNIT_FACTOR = -31536000000.0  # kg m-2 s-1  ->  gC m-2 yr-1 (sea-to-air)


def main(cfg):
    """Run the diagnostic."""
    input_data = [
        m for m in cfg["input_data"].values() if m["short_name"] == "fgco2"
    ]
    years = rc.year_range_str(input_data)

    fig, axes = plt.subplots(figsize=(9, 6.5))
    ancestors, nc_files = [], []

    for meta in input_data:
        cube = rc.load_cube(meta["filename"], meta["short_name"])
        cube = rc.time_mean(cube)
        data = np.ma.masked_invalid(cube.data) * UNIT_FACTOR
        zonal = rc.zonal_mean(data)
        lat = rc.lat_1d(cube)
        style = rc.style_for(meta["dataset"], cfg.get("styleset", "CMIP5"))
        axes.plot(
            lat,
            zonal,
            label=meta["dataset"],
            color=style["color"],
            linestyle=style["dash"],
            linewidth=style["thick"],
        )

        out_cube = iris.cube.Cube(
            zonal,
            var_name=meta["short_name"],
            units="g m-2 yr-1",
            long_name="zonal mean CO2 flux (sea to air)",
        )
        rc.add_latitude_coord(out_cube, lat)
        nc_name = get_diagnostic_filename(
            f"russell18jgr_fig-7h_fgco2_{meta['dataset']}_"
            f"{meta['start_year']}-{meta['end_year']}",
            cfg,
        )
        iris.save(out_cube, nc_name)
        nc_files.append(nc_name)
        ancestors.append(meta["filename"])

    axes.set_xlim(-80, -30)
    axes.set_ylim(-60, 40)
    axes.set_xticks(np.arange(-80, -29, 5))
    axes.set_yticks(np.arange(-60, 41, 10))
    axes.axhline(0, color="grey", linestyle="--", linewidth=0.8)
    axes.set_xlabel("Latitude")
    axes.set_title("Zonal-mean Flux", fontsize=12)
    axes.text(
        0.0,
        1.05,
        "Russell et al -2018 - Figure 7 h",
        transform=axes.transAxes,
        fontsize=10,
    )
    axes.text(
        1.0,
        1.05,
        "Units - gC/ (m$^2$ * yr)",
        ha="right",
        transform=axes.transAxes,
        fontsize=9,
    )
    axes.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=7,
        frameon=False,
    )
    fig.tight_layout()

    plot_file = get_plot_filename(f"Russell_figure7h_fgco2_{years}", cfg)
    fig.savefig(plot_file, bbox_inches="tight", dpi=200)
    plt.close(fig)
    logger.info("Wrote %s", plot_file)

    record = rc.provenance_record(
        "Russell et al 2018 figure 7h", ancestors, plot_types=("zonal",)
    )
    with ProvenanceLogger(cfg) as prov:
        for filename in nc_files + [plot_file]:
            prov.log(filename, record)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
