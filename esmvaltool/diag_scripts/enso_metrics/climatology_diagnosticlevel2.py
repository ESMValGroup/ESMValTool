"""Plot map comparisons of background climatology."""

import logging
import os
from pprint import pformat

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
import numpy as np
from esmvalcore.preprocessor import convert_units

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
)

logger = logging.getLogger(os.path.basename(__file__))


def plotmaps_level2(input_data, grp):
    """Create map plots for pair of input data."""
    fig = plt.figure(figsize=(18, 6))
    proj = ccrs.Orthographic(central_longitude=210.0)
    for plt_pos, dataset in enumerate(input_data, start=121):
        logger.info(
            "dataset: %s - %s",
            dataset["dataset"],
            dataset["long_name"],
        )
        cube, cbar_label = load_seacycle_stdev(dataset)

        ax1 = plt.subplot(plt_pos, projection=proj)
        ax1.add_feature(cfeature.LAND, facecolor="gray")
        ax1.coastlines()
        # set levels for tauu map
        if grp == "tauu_bias":
            cf1 = iplt.contourf(
                cube,
                cmap="coolwarm",
                extend="both",
                levels=np.arange(-100, 100, 20),
            )
        else:
            cf1 = iplt.contourf(cube, cmap="coolwarm", extend="both")

        ax1.set_extent([130, 290, -20, 20], crs=ccrs.PlateCarree())
        ax1.set_title(dataset["dataset"])

        # Add gridlines for latitude and longitude
        gl1 = ax1.gridlines(draw_labels=True, linestyle="--")
        gl1.top_labels = False
        gl1.right_labels = False

    # Add a single colorbar at the bottom
    cax = plt.axes([0.15, 0.08, 0.7, 0.05])
    cbar = fig.colorbar(cf1, cax=cax, orientation="horizontal", extend="both")
    cbar.set_label(cbar_label)

    return fig


def load_seacycle_stdev(dataset):
    """Load, seasonal cycle std dev if required."""
    var_units = {"tos": "degC", "pr": "mm/day", "tauu": "1e-3 N/m2"}
    sname = dataset["short_name"]
    cube = iris.load_cube(dataset["filename"])
    # convert units for different variables
    cube = convert_units(cube, units=var_units[sname])

    diag_label = sname.upper()
    if len(cube.coords("month_number")) == 1:
        cube.coord("month_number").guess_bounds()
        cube = cube.collapsed("month_number", iris.analysis.STD_DEV)
        diag_label = f"{sname.upper()} std"

    cbar_label = f"{diag_label} {var_units[sname]}"

    return cube, cbar_label


def provenance_record(var_grp, ancestor_files):
    """Create a provenance record describing the diagnostic plot."""
    caption = {
        "pr_bias": (
            "Time-mean precipitation bias in the equatorial Pacific, "
            "primarily highlighting the double intertropical convergence "
            "zone (ITCZ) bias."
        ),
        "pr_seacycle": (
            "Bias in the amplitude of the mean seasonal cycle of "
            "precipitation in the equatorial Pacific."
        ),
        "sst_bias": (
            "Time-mean sea surface temperature bias in the equatorial Pacific."
        ),
        "sst_seacycle": (
            "Bias in the amplitude of the mean seasonal cycle of "
            "sea surface temperature in the equatorial Pacific."
        ),
        "tauu_bias": "Time-mean zonal wind stress bias in the "
        "equatorial Pacific.",
        "tauu_seacycle": (
            "Bias in the amplitude of the mean seasonal cycle of "
            "zonal wind stress in the equatorial Pacific."
        ),
    }
    record = {
        "caption": caption[var_grp],
        "authors": [
            "chun_felicity",
            "beucher_romain",
        ],
        "references": [
            "planton2021",
        ],
        "ancestors": ancestor_files,
    }
    return record


def save_plotdata(plotdata, group, pairs, cfg):
    """Save both obs and model plotted data."""
    for i, cube in enumerate(plotdata):
        data_prov = provenance_record(group, [pairs[i]["filename"]])
        datafile = [
            pairs[i]["dataset"],
            pairs[i]["short_name"],
            pairs[i]["preprocessor"],
        ]
        save_data("_".join(datafile), data_prov, cfg, cube)


def main(cfg):
    """Run basic climatology diagnostic level 2 for each input dataset."""
    input_data = cfg["input_data"].values()

    # group by variables
    variable_groups = group_metadata(
        input_data,
        "variable_group",
        sort="project",
    )
    # for each select obs and iterate others, obs last
    for grp, var_attr in variable_groups.items():
        logger.info("%s : %d, %s", grp, len(var_attr), pformat(var_attr))
        prov = provenance_record(grp, list(cfg["input_data"].keys()))
        for metadata in var_attr:
            # create pairs, add obs first to list
            pairs = [var_attr[-1]]
            logger.info("iterate though datasets\n %s", pformat(metadata))
            if metadata["project"] == "CMIP6":
                pairs.append(metadata)
                fig = plotmaps_level2(pairs, grp)
                filename = "_".join(
                    [
                        metadata["dataset"],
                        metadata["short_name"],
                        metadata["preprocessor"],
                    ],
                )
                save_figure(
                    filename,
                    prov,
                    cfg,
                    figure=fig,
                    dpi=300,
                )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
