"""Annual cycle of Southern Ocean sea ice area (russell18jgr figure 5g).

Python port of russell18jgr-fig5g.ncl.

- Builds the monthly climatology of sic.
- Multiplies with the cell area (areacello fx file for native ocean
  grids; for regular lat-lon grids the area is computed manually, as
  some models provide sic on a different grid than areacello).
- Sums all cells south of the equator and plots the annual cycle in
  10^12 m^2.
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

MONTH_LABELS = [
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
]


def main(cfg):
    """Run the diagnostic."""
    metadata = list(cfg["input_data"].values())
    sic_data = [m for m in metadata if m["short_name"] in ("sic", "siconc")]
    area_data = rc.select_by_short_name(metadata, "areacello")
    var0 = sic_data[0]["short_name"]
    years = rc.year_range_str(sic_data)

    fig, axes = plt.subplots(figsize=(9, 6.5))
    months = np.arange(12)
    ancestors, nc_files = [], []

    for meta in sic_data:
        cube = rc.load_cube(meta["filename"], meta["short_name"])
        clim = rc.monthly_climatology(cube)
        sic = rc.sea_ice_percent(clim)

        lat = rc.coord_1d_or_2d(clim, "latitude")
        if rc.has_regular_grid(clim):
            area = rc.ncl_cell_area(lat, rc.lon_1d(clim), sic.shape[1:])
            sh_rows = lat < 0.0
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
            sh_rows = lat[:, 0] < 0.0

        ice_area = sic * area[np.newaxis, :, :]
        # 1e14 = 1e12 (unit scale) * 100 (percent -> fraction)
        annual_cycle = (
            ice_area[:, sh_rows, :].reshape(12, -1).sum(axis=1) / 1.0e14
        )

        style = rc.style_for(meta["dataset"], cfg.get("styleset", "CMIP5"))
        axes.plot(
            months,
            annual_cycle,
            label=meta["dataset"],
            color=style["color"],
            linestyle=style["dash"],
            linewidth=style["thick"],
            marker="D",
            markersize=3,
            markerfacecolor="none",
        )

        out_cube = iris.cube.Cube(
            annual_cycle,
            var_name=var0,
            units="1e12 m2",
            long_name="southern hemisphere area under sea ice",
        )
        out_cube.add_dim_coord(
            iris.coords.DimCoord(np.arange(1, 13), var_name="month"), 0
        )
        nc_name = get_diagnostic_filename(
            f"russell18jgr-fig5g_{var0}_{meta['dataset']}_"
            f"{meta['start_year']}-{meta['end_year']}",
            cfg,
        )
        iris.save(out_cube, nc_name)
        nc_files.append(nc_name)
        ancestors.append(meta["filename"])

    axes.set_ylim(0, 24)
    axes.set_xticks(months)
    axes.set_xticklabels(MONTH_LABELS)
    axes.set_xlabel("months")
    axes.set_ylabel("Area under sea ice ( 10$^{12}$ m$^2$ )")
    axes.set_title("Russell et al -2018 - Figure 5 g", loc="left", fontsize=11)
    axes.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=7,
        frameon=False,
    )
    fig.tight_layout()

    plot_file = get_plot_filename(f"russell18jgr-fig5g_{var0}_{years}", cfg)
    fig.savefig(plot_file, bbox_inches="tight", dpi=200)
    plt.close(fig)
    logger.info("Wrote %s", plot_file)

    record = rc.provenance_record(
        "Russell et al 2018 figure 5g", ancestors, plot_types=("times",)
    )
    with ProvenanceLogger(cfg) as prov:
        for filename in nc_files + [plot_file]:
            prov.log(filename, record)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
