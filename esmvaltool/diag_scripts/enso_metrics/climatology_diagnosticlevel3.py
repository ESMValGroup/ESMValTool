"""Plot monthly comparisons of background climatology."""

import logging
import os
from pprint import pformat

import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
import numpy as np
from esmvalcore.preprocessor import (
    convert_units,
    extract_region,
    meridional_statistics,
    zonal_statistics,
)

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
)

logger = logging.getLogger(os.path.basename(__file__))


def plotmaps_level3(input_data, itcz=False):
    """Create month structure plots for pair of input data.

    Args:
    ----
    input_data (List) : List of input datasets to plot
    itcz (bool) : Boolean flag to trigger longitude formatting

    """
    fig = plt.figure(figsize=(14, 8))
    colmap = {
        "pr": "YlGn",
        "tos": "inferno",
        "ts": "inferno",
        "tauu": "RdBu_r",
    }

    levels = {
        "pr": np.arange(0, 12, 1),
        "tos": np.arange(20, 31, 1),
        "ts": np.arange(20, 31, 1),
        "tauu": np.arange(-80, 90, 10),
    }
    for plt_pos, dataset in enumerate(input_data, start=121):
        logger.info(
            "dataset: %s - %s",
            dataset["dataset"],
            dataset["long_name"],
        )
        cube, cbar_label, x_label = load_seacycle_stat(dataset, itcz)

        ax1 = plt.subplot(plt_pos)
        cf1 = iplt.contourf(
            cube,
            coords=[x_label, "month_number"],
            levels=levels[dataset["short_name"]],
            extend="both",
            cmap=colmap[dataset["short_name"]],
        )
        ax1.set_ylim(1, 12)
        ax1.set_yticks(ticks=np.arange(1, 13, 4), labels=["Jan", "May", "Sep"])
        ax1.set_title(dataset["dataset"])
        ax1.set_ylabel("months")
        if not itcz:
            ax1.xaxis.set_major_formatter(plt.FuncFormatter(format_longitude))
        ax1.set_xlabel(x_label)

    plt.subplots_adjust(bottom=0.2)
    # Add a single colorbar at the bottom
    cax = plt.axes([0.15, 0.07, 0.7, 0.05])
    cbar = fig.colorbar(cf1, cax=cax, orientation="horizontal", extend="both")
    cbar.set_label(cbar_label)

    return fig


def load_seacycle_stat(dataset, itcz=False):
    """Load, seasonal cycle std dev if required.

    Args:
    ----
    dataset (dict)) : Dictionary of dataset metadata from recipe
    itcz (bool) : Boolean flag to trigger zonal statistics

    """
    var_units = {
        "tos": "degC",
        "ts": "degC",
        "pr": "mm/day",
        "tauu": "1e-3 N/m2",
    }
    sname = dataset["short_name"]
    cube = iris.load_cube(dataset["filename"])
    # convert units for different variables
    cube = convert_units(cube, units=var_units[sname])

    diag_label = sname.upper()
    if itcz:
        nino3_latext_region = {
            "start_longitude": 210.0,
            "end_longitude": 270.0,
            "start_latitude": -15.0,
            "end_latitude": 15.0,
        }
        x_label = "latitude"
        cube = extract_region(cube, **nino3_latext_region)
        cube = zonal_statistics(cube, "mean")
    else:
        eq_region = {
            "start_longitude": 160.0,
            "end_longitude": 270.0,
            "start_latitude": -5.0,
            "end_latitude": 5.0,
        }
        x_label = "longitude"
        cube = extract_region(cube, **eq_region)
        cube = meridional_statistics(cube, "mean")

    cbar_label = f"{diag_label} {var_units[sname]}"

    return cube, cbar_label, x_label


def format_longitude(x, _pos):
    """Format longitude values for plotting."""
    if x > 180:
        return f"{int(360 - x)}°W"
    if x == 180:
        return f"{int(x)}°"
    return f"{int(x)}°E"


def provenance_record(var_grp, ancestor_files):
    """Create a provenance record describing the diagnostic plot."""
    caption = {
        "doubleITCZ_seacycle": (
            "Meridional structure of the mean seasonal cycle of "
            "precipitation (PR) in the eastern Pacific (150-90°W averaged)."
        ),
        "pr_seacycle": (
            "Zonal structure of the mean seasonal cycle of "
            "precipitation (PR) in the eastern Pacific (5°S-5°N averaged)."
        ),
        "sst_seacycle": (
            "Zonal structure of the mean seasonal cycle of "
            "sea surface temperature (SST) in the eastern Pacific (5°S-5°N averaged)."
        ),
        "tauu_seacycle": "Zonal structure of the mean seasonal cycle of "
        "zonal wind stress (TAUX) in the eastern Pacific (5°S-5°N averaged).",
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
    """Run basic climatology diagnostic level 3 for each input dataset."""
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
        if grp in ["pr_seacycle", "sst_seacycle", "tauu_seacycle"]:
            for metadata in var_attr:
                # create pairs, add obs first to list
                pairs = [var_attr[-1]]
                logger.info("iterate though datasets\n %s", pformat(metadata))
                if metadata["project"] == "CMIP6":
                    pairs.append(metadata)
                    fig = plotmaps_level3(pairs, itcz=False)
                    # save_plotdata(data_cubes, grp, pairs, cfg)
                    filename = "_".join(
                        [
                            metadata["dataset"],
                            grp,
                            "level3",
                        ],
                    )
                    prov = provenance_record(
                        grp,
                        list(cfg["input_data"].keys()),
                    )
                    save_figure(
                        filename,
                        prov,
                        cfg,
                        figure=fig,
                        dpi=300,
                    )
                    if grp == "pr_seacycle":
                        # replace pr with doubleITCZ
                        grp_itcz = "doubleITCZ_seacycle"
                        fig = plotmaps_level3(pairs, itcz=True)
                        filename = "_".join(
                            [
                                metadata["dataset"],
                                grp_itcz,
                                "level3",
                            ],
                        )
                        prov = provenance_record(
                            grp_itcz,
                            list(cfg["input_data"].keys()),
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
