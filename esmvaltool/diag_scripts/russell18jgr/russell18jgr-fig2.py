"""Zonal and annual mean zonal wind stress (russell18jgr figure 2).

Python port of russell18jgr-fig2.ncl.

- Expects time-averaged tauu/tauuo from the ESMValTool preprocessor.
- Takes the zonal average and plots it against latitude for every
  dataset in one panel, using the CMIP5 style set.
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


def main(cfg):
    """Run the diagnostic."""
    input_data = list(cfg["input_data"].values())
    years = rc.year_range_str(input_data)

    fig, axes = plt.subplots(figsize=(9, 6))
    ancestors, nc_files = [], []

    for meta in input_data:
        cube = rc.load_cube(meta["filename"], meta["short_name"])
        cube = rc.time_mean(cube)
        zonal = rc.zonal_mean(cube.data)
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
            np.ma.masked_invalid(zonal),
            var_name=meta["short_name"],
            units=str(cube.units),
            long_name="zonal mean surface eastward wind stress",
        )
        out_cube.add_dim_coord(
            iris.coords.DimCoord(
                lat, standard_name="latitude", units="degrees_north"
            ),
            0,
        )
        nc_name = get_diagnostic_filename(
            f"russell18jgr_fig2_{meta['short_name']}_{meta['dataset']}_"
            f"{meta['start_year']}-{meta['end_year']}",
            cfg,
        )
        iris.save(out_cube, nc_name)
        nc_files.append(nc_name)
        ancestors.append(meta["filename"])

    axes.set_xlim(-80, -30)
    axes.set_ylim(-0.1, 0.25)
    axes.set_xticks(np.arange(-80, -29, 5))
    axes.set_yticks(np.arange(-0.1, 0.26, 0.1))
    axes.axhline(0, color="grey", linestyle="--", linewidth=0.8)
    axes.set_xlabel("Latitude")
    axes.set_ylabel("Surface eastward wind stress")
    axes.set_title("Russell et al -2018 - Figure 2", loc="left", fontsize=11)
    axes.set_title("Units - (Pa)", loc="right", fontsize=9)
    axes.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=7,
        frameon=False,
    )
    fig.tight_layout()

    plot_file = get_plot_filename(f"russell18jgr_fig2_{years}", cfg)
    fig.savefig(plot_file, bbox_inches="tight", dpi=200)
    plt.close(fig)
    logger.info("Wrote %s", plot_file)

    record = rc.provenance_record("Russell et al 2018 figure 2", ancestors)
    with ProvenanceLogger(cfg) as prov:
        for filename in nc_files + [plot_file]:
            prov.log(filename, record)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
