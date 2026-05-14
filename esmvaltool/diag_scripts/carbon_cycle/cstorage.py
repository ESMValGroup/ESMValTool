#!/usr/bin/env python
"""Diagnostic script to plot figure 9.42a of IPCC AR5 chapter 9.

Description
-----------
Calculate and plot trends in CO2 Seasonal cycle amplitude

Author
------
Bettina Gier (Univ. of Bremen, Germany)

Project
-------
Eval4CMIP

Configuration options in recipe
-------------------------------
save : dict, optional
    Keyword arguments for the `fig.saveplot()` function.

"""

import logging
import os

import iris
import iris.quickplot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.lines import Line2D

# import iris.plot as iplt
from esmvaltool.diag_scripts.shared import (
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    plot,
    run_diagnostic,
    select_metadata,
)

logger = logging.getLogger(os.path.basename(__file__))


def n_models(dataset):
    n_cycle_models = [
        "ACCESS-ESM1-5",
        "BNU-ESM",
        "CESM1-BGC",
        "NorESM1-ME",
        "UKESM1-0-LL",
        "NorESM2-LM",
        "EC-Earth3-Veg",
        "EC-Earth3-CC",
        "CESM2",
        "CESM2-WACCM",
        "SAM0-UNICON",
        "MIROC-ES2L",
        "MPI-ESM1-2-LR",
    ]
    return dataset in n_cycle_models


def plot_scatter(data, cfg, proj, exp, model_chars):
    # One plot for each CMIP Project+Experiment for different highlighting
    # for project in proj:
    #    for experiment in exp:
    #        if any(d['exp'] == experiment for d in proj[project]):
    plot_path = get_plot_filename("carbon_storage_global", cfg)
    print(plot_path)
    fig, axx = plt.subplots(figsize=(10, 7))
    leg_added = []
    add_legend_elements = []
    for dataset in data:
        data_parts = str.split(dataset, "__")
        if "OBS" in data_parts[1].upper():
            print("in obs")
            style_name = data_parts[0]
        else:
            style_name = data_parts[1] + "_" + data_parts[2]
        # if ((data_parts[1]==project and data_parts[2]==experiment)
        #    or ('OBS' in data_parts[1].upper())):
        #    style_name = data_parts[0]
        #    lab = style_name
        # else:
        #    style_name = data_parts[1] + "_" + data_parts[2]
        if style_name in leg_added:
            lab = "_nolegend_"
        else:
            if data_parts[2] == "historical":
                exp = "c"
            else:
                exp = "e"
            if "OBS" not in data_parts[1].upper():
                lab = data_parts[1] + exp
        #        lab = data_parts[1] + " " + data_parts[2]
        #        leg_added.append(style_name)
        # print(dataset)
        # print(style_name)
        style = plot.get_dataset_style(style_name, "cmip6_gier21")
        if (data_parts[0] in model_chars) and (
            model_chars[data_parts[0]]["n_cycle"]
        ):
            fill_style = "full"
            style["facecolor"] = style["color"]
            alpha = 0.5
        else:
            alpha = 0.5
            fill_style = "none"
        if (
            "mmm" in data_parts[0] and "OBS" not in data_parts[1].upper()
        ):  # have all MMMs marked by Diamonds in the diff colors! and increase size
            style["mark"] = "X"
            markersize = 14
            # if data_parts[2] == 'historical':
            #    exper = 'c'
            # else:
            #    exper = 'e'
            if (data_parts[0] in model_chars) and (
                model_chars[data_parts[0]]["n_cycle"]
            ):
                suffix = "Ncycle MMM"
            else:
                suffix = "non-Ncycle MMM"
            lab = data_parts[1] + exp + " " + suffix
            leg_added.append(lab)
        elif "OBS" in data_parts[1].upper():
            # continue
            markersize = 14
            alpha = 1
            lab = style_name
            leg_added.append(lab)
            print("in obs2")
            print(lab)
            print(style_name)
        else:
            markersize = 10
            if style_name not in leg_added:
                add_legend_elements.append(
                    Line2D(
                        [0],
                        [0],
                        marker=style["mark"],
                        color=style["color"],
                        markerfacecolor=style["facecolor"],
                        fillstyle="full",
                        linestyle="none",
                        markersize=markersize,
                        markeredgewidth=2.0,
                        alpha=alpha,
                        label=lab + " Ncycle",
                    )
                )
                add_legend_elements.append(
                    Line2D(
                        [0],
                        [0],
                        marker=style["mark"],
                        color=style["color"],
                        markerfacecolor=style["facecolor"],
                        fillstyle="none",
                        linestyle="none",
                        markersize=markersize,
                        markeredgewidth=2.0,
                        alpha=alpha,
                        label=lab + " non-Ncycle",
                    )
                )
            lab = "_nolegend_"  # Only add my custom one
            leg_added.append(style_name)
        # if ((data_parts[0] in model_chars)
        #    and (model_chars[data_parts[0]]['n_cycle'])):
        #    fill_style = "full"
        #    style['facecolor'] = style['color']
        #    alpha = 0.5
        # else:
        #    alpha = 1
        #    fill_style = "none"
        # print(leg_added)
        axx.plot(
            data[dataset]["cTotal"],
            data[dataset]["cVeg"],
            marker=style["mark"],
            color=style["color"],
            markerfacecolor=style["facecolor"],
            fillstyle=fill_style,
            linestyle="none",
            markersize=markersize,
            markeredgewidth=2.0,
            alpha=alpha,
            label=lab,
        )
    axx.set_xlabel("Soil Carbon (+ Litter) [PgC]")
    axx.set_ylabel("Vegetation Carbon [PgC]")
    handles, labels = axx.get_legend_handles_labels()
    add_legend_elements.extend(
        handles
    )  # Make custom legend entries appear first, then premade ones
    lgd = axx.legend(
        handles=add_legend_elements,
        bbox_to_anchor=(1, 1.01),
        loc="upper left",
        fontsize=12,
    )
    # axx.legend(bbox_to_anchor=(1,1), loc="upper left")
    plt.tight_layout()
    fig.savefig(
        plot_path, bbox_extra_artists=(lgd,), dpi=300, bbox_inches="tight"
    )
    plt.close()


def main(cfg):
    """Run the diagnostic."""

    # Get model carbon cycle characteristics from yaml file
    base_dir = os.path.dirname(os.path.realpath(__file__))
    file = base_dir + "/carbon_model_characteristics.yml"
    with open(file) as infile:
        model_chars = yaml.safe_load(infile)

    # Select data from the different variables
    data = {}
    data["cVeg"] = select_metadata(
        cfg["input_data"].values(), short_name="cVeg"
    )
    data["cSoil"] = select_metadata(
        cfg["input_data"].values(), short_name="cSoil"
    )
    # print(data['cSoil'])
    data["cSoil"] = data["cSoil"] + select_metadata(
        cfg["input_data"].values(), short_name="cSoilAbove1m"
    )

    # print(data['cSoil'])
    data["cLitter"] = select_metadata(
        cfg["input_data"].values(), short_name="cLitter"
    )
    plot_data = {}
    save_data = {}

    input_data = cfg["input_data"].values()
    projects = group_metadata(input_data, "project")
    cmip_projects = [proj for proj in projects.keys() if "CMIP" in proj]
    exps = group_metadata(input_data, "exp")
    experiments = list(exps.keys())
    mmm_cubes = {}
    for proj in cmip_projects:
        for exp in experiments:
            str1 = "__".join(["Ncycle_mmm", proj, exp])
            mmm_cubes[str1] = []
            str2 = "__".join(["non-Ncycle_mmm", proj, exp])
            mmm_cubes[str2] = []

    # Read data and collect
    for var in data:
        logger.info("Processing %s", var)
        for dataset in data[var]:
            name = dataset["dataset"]
            if (
                name == "HadGEM2-CC" or name == "HadGEM2-ES"
            ):  # this model has dec 2005 in the rcp files
                dataset["exp"] = "historical"
            long_name = (
                name + "__" + dataset["project"] + "__" + dataset["exp"]
            )
            logger.info("Processing %s", name)
            if (
                "reference_dataset" in dataset
                and name == dataset["reference_dataset"]
            ):
                long_name = "__".join(
                    ["OBS", dataset["project"], dataset["exp"]]
                )
                # print(name)
                # print(long_name)
            if long_name not in plot_data.keys():
                plot_data[long_name] = {}
                if var == "cVeg":
                    save_data[
                        (long_name.replace("__", " ")).replace("_", r"\_")
                    ] = {}
            cube_data = iris.load_cube(dataset["filename"])
            plot_data[long_name][var] = cube_data.data
            if var == "cVeg":
                save_data[(long_name.replace("__", " ")).replace("_", r"\_")][
                    var
                ] = cube_data.data
            if dataset["project"] in cmip_projects:
                if not (
                    ("exclude_mmm" in cfg) and (name in cfg["exclude_mmm"])
                ):
                    if (name in model_chars) and (
                        model_chars[name]["n_cycle"]
                    ):
                        str = "__".join(
                            ["Ncycle_mmm", dataset["project"], dataset["exp"]]
                        )
                    else:
                        str = "__".join(
                            [
                                "non-Ncycle_mmm",
                                dataset["project"],
                                dataset["exp"],
                            ]
                        )
                    mmm_cubes[str].append(long_name)
            if name == "HWSD":
                plot_data[long_name][var] = (
                    1300  # FIX WHY THIS ONE ONLY GOES WRONG! cube_data.data = inf??
                )
                save_data[long_name.replace("__", " ").replace("_", r"\_")][
                    var
                ] = 1300.0

    # Add cSoil and cLitter
    for dataset in plot_data:
        if "cLitter" in plot_data[dataset]:
            plot_data[dataset]["cTotal"] = (
                plot_data[dataset]["cSoil"] + plot_data[dataset]["cLitter"]
            )
            save_data[(dataset.replace("__", " ")).replace("_", r"\_")][
                "cTotal"
            ] = plot_data[dataset]["cSoil"] + plot_data[dataset]["cLitter"]
        else:
            plot_data[dataset]["cTotal"] = plot_data[dataset]["cSoil"]
            save_data[(dataset.replace("__", " ")).replace("_", r"\_")][
                "cTotal"
            ] = plot_data[dataset]["cSoil"]

    # add n-cycle mms to plot_data
    for long_name in mmm_cubes:
        # Subset plot_data according to this mmm
        cvegs = []
        ctotals = []
        for dataset in mmm_cubes[long_name]:
            cvegs.append(plot_data[dataset]["cVeg"])
            ctotals.append(plot_data[dataset]["cTotal"])
        if len(cvegs) > 1:
            plot_data[long_name] = {}
            plot_data[long_name]["cVeg"] = np.mean(cvegs)
            plot_data[long_name]["cTotal"] = np.mean(ctotals)

            save_data[(long_name.replace("__", " ")).replace("_", r"\_")] = {}
            save_data[(long_name.replace("__", " ")).replace("_", r"\_")][
                "cVeg"
            ] = np.mean(cvegs)
            save_data[(long_name.replace("__", " ")).replace("_", r"\_")][
                "cTotal"
            ] = np.mean(ctotals)

    # Plot
    plot_scatter(plot_data, cfg, projects, experiments, model_chars)

    # Save data
    path_scatter_latex = get_diagnostic_filename(
        "data_scatter_csoil_cveg", cfg, extension="tex"
    )

    # logger.info(save_data)
    dfs = pd.DataFrame.from_dict(
        save_data, orient="index", columns=["cVeg", "cTotal"], dtype="float"
    )

    dfs.to_latex(
        path_scatter_latex,
        float_format="%.0f",
        longtable=True,
        caption="Carbon Storage Totals [PgC]",
        label="tab:csoil-cveg",
    )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
