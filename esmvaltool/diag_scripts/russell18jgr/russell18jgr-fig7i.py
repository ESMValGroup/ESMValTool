"""Cumulative integrated CO2 flux 90S-30S (russell18jgr figure 7i).

Python port of russell18jgr-fig7i.ncl.

- Expects time-averaged, land-masked fgco2 from the preprocessor.
- Multiplies the flux with the cell area (areacello fx file for native
  ocean grids, or the manual regular-grid area calculation used by the
  NCL script) and accumulates the sum from the south pole northwards.
- Plots the cumulative integral (PgC/yr) against latitude.
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

UNIT_FACTOR = -0.000031536  # kg s-1  ->  PgC yr-1 (sea-to-air)


def main(cfg):
    """Run the diagnostic."""
    metadata = list(cfg["input_data"].values())
    fgco2_data = rc.select_by_short_name(metadata, "fgco2")
    area_data = rc.select_by_short_name(metadata, "areacello")
    years = rc.year_range_str(fgco2_data)

    fig, axes = plt.subplots(figsize=(9, 6.5))
    ancestors, nc_files = [], []

    for meta in fgco2_data:
        cube = rc.load_cube(meta["filename"], meta["short_name"])
        cube = rc.time_mean(cube)
        data = np.ma.masked_invalid(cube.data)
        lat = rc.lat_1d(cube)

        if rc.has_regular_grid(cube):
            area = rc.ncl_cell_area(lat, rc.lon_1d(cube), data.shape)
        else:
            area_meta = rc.match_by_dataset(area_data, meta["dataset"])
            if area_meta is None:
                raise ValueError(
                    f"areacello file for {meta['dataset']} not found; "
                    "please add the dataset to the areacello section "
                    "of the recipe."
                )
            area = np.ma.masked_invalid(
                rc.load_cube(area_meta["filename"], "areacello").data
            )

        carbon_flux = data * UNIT_FACTOR * area
        per_lat = carbon_flux.sum(axis=-1)
        cumulative = np.ma.cumsum(per_lat)

        style = rc.style_for(meta["dataset"], cfg.get("styleset", "CMIP5"))
        axes.plot(
            lat,
            cumulative,
            label=meta["dataset"],
            color=style["color"],
            linestyle=style["dash"],
            linewidth=style["thick"],
        )

        out_cube = iris.cube.Cube(
            cumulative,
            var_name="fgco2",
            units="Pg yr-1",
            long_name="cumulative integrated carbon flux from 90S",
        )
        out_cube.add_dim_coord(
            iris.coords.DimCoord(
                np.asarray(lat, dtype=float),
                standard_name="latitude",
                units="degrees_north",
            ),
            0,
        )
        nc_name = get_diagnostic_filename(
            f"russell_figure-7i_fgco2_{meta['dataset']}_"
            f"{meta['start_year']}-{meta['end_year']}",
            cfg,
        )
        iris.save(out_cube, nc_name)
        nc_files.append(nc_name)
        ancestors.append(meta["filename"])

    axes.set_xlim(-80, -30)
    axes.set_ylim(-1.2, 1.0)
    axes.set_xticks(np.arange(-80, -29, 5))
    axes.set_yticks(np.arange(-1.2, 1.01, 0.2))
    axes.axhline(0, color="grey", linestyle="--", linewidth=0.8)
    axes.set_title("Integrated Flux", fontsize=12)
    axes.text(
        0.0,
        1.05,
        "Russell et al -2018 - Figure 7i ",
        transform=axes.transAxes,
        fontsize=10,
    )
    axes.text(
        1.0,
        1.05,
        "Units - ( PgC/yr )",
        ha="right",
        transform=axes.transAxes,
        fontsize=9,
    )
    axes.set_xlabel("Latitude")
    axes.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=7,
        frameon=False,
    )
    fig.tight_layout()

    plot_file = get_plot_filename(f"Russell_figure7i_fgco2_{years}", cfg)
    fig.savefig(plot_file, bbox_inches="tight", dpi=200)
    plt.close(fig)
    logger.info("Wrote %s", plot_file)

    record = rc.provenance_record(
        "Russell et al 2018 figure 7i", ancestors, plot_types=("zonal",)
    )
    with ProvenanceLogger(cfg) as prov:
        for filename in nc_files + [plot_file]:
            prov.log(filename, record)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
