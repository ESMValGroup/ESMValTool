"""Westerly band width vs. carbon uptake (russell18jgr figure 9b).

Python port of russell18jgr-fig9b.ncl.

Scatter plot of the width of the Southern Hemisphere westerly wind
band against the annual-mean integrated carbon uptake south of 30S
(fgco2, in Pg C/yr), with the line of best fit.
"""

import logging
import os
import sys

import numpy as np

from esmvaltool.diag_scripts.shared import run_diagnostic

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import russell_common as rc
import russell_fig9_shared as f9

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    """Run the diagnostic."""
    metadata = list(cfg["input_data"].values())
    fgco2_data = rc.select_by_short_name(metadata, "fgco2")
    area_data = rc.select_by_short_name(metadata, "areacello")
    tauu_data = [m for m in metadata if m["short_name"] in ("tauu", "tauuo")]
    years = rc.year_range_str(fgco2_data)

    datasets, widths, fluxes, metas, ancestors = [], [], [], [], []
    for tauu_meta in tauu_data:
        dataset = tauu_meta["dataset"]
        fgco2_meta = rc.match_by_dataset(fgco2_data, dataset)
        if fgco2_meta is None:
            logger.warning("No fgco2 for %s, skipping", dataset)
            continue
        widths.append(f9.band_width(tauu_meta))
        fluxes.append(
            f9.integrated_flux(fgco2_meta, area_data, f9.CARBON_FACTOR)
        )
        datasets.append(dataset)
        metas.append(fgco2_meta)
        ancestors.append([tauu_meta["filename"], fgco2_meta["filename"]])

    xlim = (
        np.round(min(widths) * 2) / 2.0 - 0.5,
        np.round(max(widths) * 2) / 2.0 + 0.5,
    )
    plot_file, regline = f9.scatter_with_regression(
        cfg,
        widths,
        fluxes,
        datasets,
        "Latitudinal width of Southern Hemisphere Westerly Band",
        "Southern ocean carbon uptake (Pg/yr)",
        "Russell et at 2018 - fig 9b ",
        xlim,
        0.5,
        f"russell18jgr-fig9b_{years}",
    )

    f9.save_pairs(
        cfg,
        "russell18jgr_fig9b",
        "carbon-flux_lat-width",
        "total carbon flux and lat width of southern westerly band",
        datasets,
        metas,
        [(fluxes[i], widths[i]) for i in range(len(datasets))],
        regline,
        plot_file,
        "Russell et al 2018 figure 9b",
        ancestors,
    )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
