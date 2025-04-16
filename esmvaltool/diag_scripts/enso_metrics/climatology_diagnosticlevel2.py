"""Plot map comparisons of background climatology."""

import logging
import os
from pprint import pformat

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
from esmvalcore.preprocessor import convert_units

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_figure,
)

logger = logging.getLogger(os.path.basename(__file__))


def plotmaps_level2(input_data):
    """Create map plots for pair of input data."""
    var_units = {"tos": "degC", "pr": "mm/day", "tauu": "1e-3 N/m2"}
    fig = plt.figure(figsize=(18, 6))
    proj = ccrs.Orthographic(central_longitude=210.0)

    for plt_pos, dataset in enumerate(input_data, start=121):
        sname = dataset["short_name"]

        logger.info(
            "dataset: %s - %s", dataset["dataset"], dataset["long_name"],
        )

        cube = iris.load_cube(dataset["filename"])
        # convert units for different variables
        cube = convert_units(cube, units=var_units[sname])
        diag_label = sname.upper()
        if len(cube.coords("month_number")) == 1:
            cube = sea_cycle_stdev(cube)
            diag_label = f"{sname.upper()} std"

        cbar_label = f"{diag_label} {var_units[sname]}"
        ax1 = plt.subplot(plt_pos, projection=proj)
        ax1.add_feature(cfeature.LAND, facecolor="gray")
        ax1.coastlines()
        cf1 = iplt.contourf(cube, cmap="coolwarm")

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


def sea_cycle_stdev(cube):
    """Process seasonal cycle standard deviation."""
    cube.coord("month_number").guess_bounds()
    cube = cube.collapsed("month_number", iris.analysis.STD_DEV)

    return cube


def main(cfg):
    """Run basic climatology diagnostic level 2 for each input dataset."""
    provenance_record = {
        "caption": "ENSO metrics comparison maps",
        "authors": [
            "chun_felicity",
            "beucher_romain",
        ],
        "references": [""],
        "ancestors": list(cfg["input_data"].keys()),
    }
    input_data = cfg["input_data"].values()

    # group by variables
    variable_groups = group_metadata(
        input_data, "variable_group", sort="project",
    )
    # for each select obs and iterate others, obs last
    for grp, var_attr in variable_groups.items():
        # create pairs
        logger.info("%s : %d, %s", grp, len(var_attr), pformat(var_attr))
        obs_data = var_attr[-1]

        for metadata in var_attr:
            logger.info("iterate though datasets\n %s", pformat(metadata))
            pairs = [obs_data]
            if metadata["project"] == "CMIP6":
                pairs.append(metadata)
                fig = plotmaps_level2(pairs)
                filename = "_".join(
                    [
                        metadata["dataset"],
                        metadata["short_name"],
                        metadata["preprocessor"],
                    ],
                )
                save_figure(
                    filename, provenance_record, cfg, figure=fig, dpi=300,
                )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
