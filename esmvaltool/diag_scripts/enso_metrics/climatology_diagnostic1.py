"""Compute basic climatology metrics for ENSO."""

import logging
import os
from pprint import pformat

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np
from esmvalcore.preprocessor import (
    convert_units,
    meridional_statistics,
    zonal_statistics,
)

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
)

logger = logging.getLogger(os.path.basename(__file__))


def plot_level1(input_data, cfg, provenance):
    """Create plots for pair of input data."""
    plt.clf()
    obs_data, model_data = None, None
    figure = plt.figure(figsize=(10, 6), dpi=300)

    for dataset in input_data:
        logger.info(
            "dataset: %s - %s",
            dataset["dataset"],
            dataset["long_name"],
        )
        cube, title, ylabel = load_process_seacycle(dataset)
        # save data
        save_data(
            f"{dataset['dataset']}_{dataset['variable_group']}",
            provenance,
            cfg,
            cube,
        )

        if dataset["project"] == "CMIP6":
            qplt.plot(cube, label=dataset["dataset"])
            model_data = cube.data
        else:
            qplt.plot(cube, label=f"ref: {dataset['dataset']}", color="black")
            obs_data = cube.data

    rmse = np.sqrt(np.mean((obs_data - model_data) ** 2))
    filename = [input_data[1]["dataset"], input_data[1]["variable_group"]]

    plt.title(title)
    plt.legend()
    plt.grid(linestyle="--")
    plt.ylabel(ylabel)

    plt.text(
        0.5,
        0.95,
        f"RMSE: {rmse:.2f} {cube.units}",
        fontsize=12,
        ha="center",
        transform=plt.gca().transAxes,
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
    )

    if input_data[1]["preprocessor"].startswith("ITCZ"):
        plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(format_lat))
    else:
        plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(format_lon))

    return figure, filename, rmse


def load_process_seacycle(dataset):
    """Load, standard deviation for seasonal cycle if req."""
    var_units = {"tos": "degC", "pr": "mm/day", "tauu": "1e-3 N/m2"}
    sname = dataset["short_name"]
    cube = iris.load_cube(dataset["filename"])
    # convert units for different variables
    cube = convert_units(cube, units=var_units[sname])

    title = f"Mean {dataset['long_name']}"
    ylabel = f"{sname.upper()} ({cube.units})"
    if len(cube.coords("month_number")) == 1:
        cube.coord("month_number").guess_bounds()
        cube = cube.collapsed("month_number", iris.analysis.STD_DEV)
        # ITCZ or zonal
        if dataset["preprocessor"].startswith("ITCZ"):
            cube = zonal_statistics(cube, "mean")
        else:
            cube = meridional_statistics(cube, "mean")

        ylabel = f"{sname.upper()} std ({cube.units})"
        title = f"{dataset['long_name']} seasonal cycle"

    return cube, title, ylabel


def format_lat(x_val, _) -> str:
    """Format latitude for plot axis."""
    if x_val < 0:
        return f"{abs(x_val):.0f}°S"
    if x_val > 0:
        return f"{x_val:.0f}°N"

    return "0°"


def format_lon(val, _) -> str:
    """Format longitude for plot axis."""
    if val > 180:
        return f"{(360 - val):.0f}°W"
    if val == 180:
        return f"{val:.0f}°"

    return f"{val:.0f}°E"


def provenance_record(var_grp, ancestor_files):
    """Create a provenance record describing the diagnostic."""
    caption = {
        "pr_double": (
            "Meridional bias in the time-mean precipitation structure "
            + "across the eastern Pacific (averaged between 150-90°W), "
            + "primarily illustrating the double intertropical convergence "
            + "zone (ITCZ) bias."
        ),
        "eq_pr_bias": (
            "Zonal bias in the time-mean precipitation structure across "
            + "the equatorial Pacific (averaged between 5°S-5°N), "
            + "illustrating the increased precipitation in the eastern "
            + "Pacific and decreased precipitation in the western Pacific."
        ),
        "eq_sst_bias": (
            "Zonal bias in the sea surface temperature structure across "
            + "the equatorial Pacific (averaged between 5°S-5°N), primarily "
            + "illustrating the cold tongue bias (typically warmer near "
            + "South America and cooler further west)."
        ),
        "eq_tauu_bias": (
            "Zonal bias in the structure of zonal wind stress across "
            + "the equatorial Pacific (averaged between 5°S-5°N), primarily "
            + "highlighting the trade winds bias (typically weaker "
            + "circulation in the central Pacific and stronger in the "
            + "western Pacific)."
        ),
        "pr_double_seacycle": (
            "Meridional bias in the amplitude of the mean seasonal "
            + "precipitation cycle in the eastern Pacific "
            + "(averaged between 150-90°W). "
        ),
        "eq_pr_seacycle": (
            "Zonal bias in the amplitude of the mean seasonal cycle of "
            + "precipitation in the equatorial Pacific "
            + "(averaged between 5°S-5°N)."
        ),
        "eq_sst_seacycle": (
            "Zonal bias in the amplitude of the mean seasonal cycle of sea "
            + "surface temperature in the equatorial Pacific "
            + "(averaged between 5°S-5°N)."
        ),
        "eq_tauu_seacycle": (
            "Zonal bias in the amplitude of the mean seasonal cycle of "
            + "zonal wind stress in the equatorial Pacific "
            + "(averaged between 5°S-5°N)."
        ),
        "values": "List of metric values.",
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


def main(cfg):
    """Compute sea ice area for each input dataset."""
    input_data = cfg["input_data"].values()
    metricfile = get_diagnostic_filename("matrix", cfg, extension="csv")
    # group by variables
    variable_groups = group_metadata(
        input_data,
        "variable_group",
        sort="project",
    )
    # for each select obs and iterate others, obs last
    for grp, var_attr in variable_groups.items():
        logger.info("%s : %d, %s", grp, len(var_attr), pformat(var_attr))
        pairs = [var_attr[-1]]  # obs to list
        prov = provenance_record(grp, list(cfg["input_data"].keys()))
        for metadata in var_attr:
            logger.info("iterate though datasets\n %s", pformat(metadata))
            if metadata["project"] == "CMIP6":
                pairs.append(metadata)
                fig, filename, rmse = plot_level1(pairs, cfg, prov)

                save_figure(
                    "_".join(filename),
                    prov,
                    cfg,
                    figure=fig,
                    dpi=300,
                )
                with open(metricfile, "a+", encoding="utf-8") as fileo:
                    fileo.write(f"{filename[0]},{filename[1]},{rmse}\n")
    # write provenance for csv metrics
    prov = provenance_record("values", list(cfg["input_data"].keys()))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(metricfile, prov)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
