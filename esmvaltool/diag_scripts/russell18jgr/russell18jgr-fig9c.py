"""Heat uptake vs. carbon uptake south of 30S (russell18jgr figure 9c).

Python port of russell18jgr-fig9c.ncl.

Scatter plot of the net heat uptake south of 30S (PW) against the
annual-mean integrated carbon uptake south of 30S (Pg C/yr), with the
line of best fit.
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
    hfds_data = rc.select_by_short_name(metadata, "hfds")
    area_data = rc.select_by_short_name(metadata, "areacello")
    years = rc.year_range_str(fgco2_data)

    datasets, heat, carbon, metas, ancestors = [], [], [], [], []
    for fgco2_meta in fgco2_data:
        dataset = fgco2_meta["dataset"]
        hfds_meta = rc.match_by_dataset(hfds_data, dataset)
        if hfds_meta is None:
            logger.warning("No hfds for %s, skipping", dataset)
            continue
        carbon.append(
            f9.integrated_flux(fgco2_meta, area_data, f9.CARBON_FACTOR)
        )
        heat.append(f9.heat_flux(hfds_meta, area_data))
        datasets.append(dataset)
        metas.append(fgco2_meta)
        ancestors.append([fgco2_meta["filename"], hfds_meta["filename"]])

    xlim = (
        np.floor((min(heat) - 0.2) * 10) / 10,
        np.ceil((max(heat) + 0.2) * 10) / 10,
    )
    plot_file, regline = f9.scatter_with_regression(
        cfg,
        heat,
        carbon,
        datasets,
        "Southern ocean heat uptake (PW)",
        "Southern ocean carbon uptake (Pg/yr)",
        " Russell et al 2018 - Figure 9c ",
        xlim,
        None,
        f"russell18jgr-fig9c_{years}",
    )

    f9.save_pairs(
        cfg,
        "russell18jgr_fig9c",
        "heat-flux_carbon-flux",
        "total heat and carbon flux south of 30S",
        datasets,
        metas,
        [(heat[i], carbon[i]) for i in range(len(datasets))],
        regline,
        plot_file,
        "Russell et al 2018 figure 9c",
        ancestors,
    )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
