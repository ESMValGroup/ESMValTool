"""Westerly band width vs. heat uptake (russell18jgr figure 9a).

Python port of russell18jgr-fig9a.ncl.

Scatter plot of the width of the Southern Hemisphere westerly wind
band against the annual-mean integrated heat uptake south of 30S
(hfds, in PW; negative is heat lost from the ocean), with the line of
best fit.
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
    hfds_data = rc.select_by_short_name(metadata, "hfds")
    area_data = rc.select_by_short_name(metadata, "areacello")
    tauu_data = [m for m in metadata
                 if m["short_name"] in ("tauu", "tauuo")]
    years = rc.year_range_str(tauu_data)

    datasets, widths, fluxes, metas, ancestors = [], [], [], [], []
    for tauu_meta in tauu_data:
        dataset = tauu_meta["dataset"]
        hfds_meta = rc.match_by_dataset(hfds_data, dataset)
        if hfds_meta is None:
            logger.warning("No hfds for %s, skipping", dataset)
            continue
        widths.append(f9.band_width(tauu_meta))
        fluxes.append(f9.heat_flux(hfds_meta, area_data))
        datasets.append(dataset)
        metas.append(tauu_meta)
        ancestors.append([tauu_meta["filename"], hfds_meta["filename"]])

    xlim = (np.floor(min(widths) - 0.5), np.ceil(max(widths) + 0.5))
    plot_file, regline = f9.scatter_with_regression(
        cfg, widths, fluxes, datasets,
        "Latitudinal width of Southern Hemisphere Westerly Band",
        "Southern ocean heat uptake (PW)",
        " Russell et al 2018 - Figure 9a ",
        xlim, 1.0, f"russell18jgr-fig9a_{years}")

    f9.save_pairs(
        cfg, "russell18jgr_fig9a", "heat-flux_lat-width",
        "total heat flux and lat width of southern westerly band",
        datasets, metas,
        [(fluxes[i], widths[i]) for i in range(len(datasets))],
        regline, plot_file, "Russell et al 2018 figure 9a", ancestors)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
