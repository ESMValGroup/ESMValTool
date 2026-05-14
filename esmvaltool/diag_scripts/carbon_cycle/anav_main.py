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
import pickle
from copy import deepcopy

import iris
import iris.pandas
import iris.quickplot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from esmvalcore.preprocessor import (
    annual_statistics,
    climate_statistics,
    linear_trend,
    linear_trend_stderr,
    multi_model_statistics,
)
from matplotlib.lines import Line2D

# import iris.plot as iplt
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    extract_variables,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    plot,
    run_diagnostic,
)
from esmvaltool.diag_scripts.shared.io import iris_save

logger = logging.getLogger(os.path.basename(__file__))


def plot_scatter(data, data_with_error, model_chars, cfg):
    # Prepare attributes
    timespan, var, var_unit, projects, cmip_projects, experiments = (
        get_plot_attributes(cfg)
    )
    if var != "lai":
        trend_unit = "PgC / y$^2$"
    else:
        trend_unit = "1"

    # One plot for each CMIP Project+Experiment for different highlighting
    # for project in cmip_projects:
    #    for experiment in experiments:
    #        if any(d['exp']==experiment for d in projects[project]):
    #            plot_path = get_plot_filename('scatter_%s_%s_%s_%s_%s'
    #                                          %(var, timespan, cfg['region'],
    #                                            project, experiment),
    #                                          cfg)

    basename = "scatter_" + var + "_" + timespan + "_" + cfg["region"]
    plot_path = get_plot_filename(
        "scatter_%s_%s_%s" % (var, timespan, cfg["region"]), cfg
    )
    pickle_path = os.path.join(
        cfg["plot_dir"],
        f"{basename}.pickle",
    )
    fig, axx = plt.subplots(figsize=(10, 7))
    leg_added = []
    add_legend_elements = []
    # Get dataset style
    for dataset in data:
        data_parts = str.split(dataset, "__")
        # if ((data_parts[1]==project and data_parts[2]==experiment)
        #    or ('OBS' in data_parts[1].upper())):
        #    style_name = data_parts[0]
        #    lab = style_name
        # else:
        style_name = data_parts[1] + "_" + data_parts[2]
        # lab = '_nolegend_' #By default don't add points to legend
        if style_name in leg_added:
            lab = "_nolegend_"
        else:
            if data_parts[2] == "historical":
                exp = "c"
            else:
                exp = "e"
            lab = data_parts[1] + exp
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
            if (data_parts[0] in model_chars) and (
                model_chars[data_parts[0]]["n_cycle"]
            ):
                suffix = "Ncycle MMM"
            else:
                suffix = "non-Ncycle MMM"
            lab = data_parts[1] + exp + " " + suffix
            leg_added.append(lab)
            alpha = 1
        elif "OBS" in data_parts[1].upper():
            continue
            markersize = 14
            alpha = 1
        else:
            markersize = 10  # Add custom markers, both N-cycle and non-Ncycle
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
        axx.plot(
            data[dataset]["mean"],
            data[dataset]["trend"],
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

    for dataset in data_with_error:
        style = plot.get_dataset_style(dataset, "cmip6_gier21")
        axx.errorbar(
            data_with_error[dataset]["mean"][0],
            data_with_error[dataset]["trend"][0],
            xerr=data_with_error[dataset]["mean"][1],
            yerr=data_with_error[dataset]["trend"][1],
            color=style["color"],
            marker=style["mark"],
            linestyle="none",
            markersize=14,
            label=dataset,
            capsize=10,
        )
    axx.set_xlabel(
        timespan + " Mean " + var.upper() + " [%s]" % var_unit, fontsize=12
    )
    axx.set_ylabel(timespan + " Linear Trend  [%s]" % trend_unit, fontsize=12)
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
    # fig_legend.legend(axx,get_legend_handles_labels())
    axx.set_title(region_label(cfg["region"]))
    pickle.dump(fig, open(pickle_path, "wb"))
    plt.tight_layout()
    fig.savefig(
        plot_path, bbox_extra_artists=(lgd,), dpi=300, bbox_inches="tight"
    )
    plt.close()


def get_style(dataset, model_chars, mmm, refs=[]):
    data_parts = str.split(dataset, "__")
    if dataset in refs:
        style_name = data_parts[0]
        lab = data_parts[0]  # dataset
    elif "mmm" in dataset and mmm:
        style_name = "_".join(data_parts[-2:])
        lab = get_label(dataset, mmm)
    else:
        style_name = data_parts[0]  # "_".join(data_parts[-2:])
        lab = get_label(dataset, mmm)
    style = plot.get_dataset_style(style_name, "cmip6_gier21")
    if "mmm" in dataset and not mmm:
        style["color"] = "red"
    style["label"] = lab
    if dataset not in refs:
        dparts_split = str.split(data_parts[0], "_")
        if model_chars[dparts_split[0]]["n_cycle"]:
            style["thick"] = 5
            if mmm:
                style["dash"] = "--"
            style["hatch"] = "----"
            style["alpha"] = 0.5
            style["fill_style"] = "full"
        else:
            if mmm:
                style["dash"] = "-"
            style["thick"] = 3
            style["hatch"] = "||||"
            style["alpha"] = 1
            style["fill_style"] = "none"
    else:
        style["hatch"] = None
    style["facecolor"] = style["color"]
    return style


def get_label(name, mmm):
    data_parts = str.split(name, "__")
    lab_parts = []
    if "mmm" in data_parts[0]:
        mmm_split = str.split(data_parts[0], "_")
        lab_parts.append(mmm_split[0])
        lab_parts.append("MMM")
        if mmm:
            if data_parts[2] == "historical":
                lab_parts.append(data_parts[1] + "c")
            elif (
                data_parts[2] == "esm-hist" or data_parts[2] == "esmHistorical"
            ):
                lab_parts.append(data_parts[1] + "e")
            else:
                lab_parts.append(data_parts[1])
                lab_parts.append(data_parts[2])
    else:
        lab_parts.append(data_parts[0])

    return " ".join(lab_parts)


def region_label(region):
    if region == "global":
        return "Global"
    elif region == "nh":
        return r"Northern Hemisphere ($20\degree$N - $90\degree$N)"
    elif region == "sh":
        return r"Southern Hemisphere ($20\degree$S - $90\degree$S)"
    elif region == "trop":
        return r"Tropics ($20\degree$S - $20\degree$N)"
    else:
        return region


def plot_project_cycle(data, model_chars, cfg):
    timespan, var, var_unit, projects, cmip_projects, experiments = (
        get_plot_attributes(cfg)
    )
    mmm = True
    refs = [dataset for dataset in list(data.keys()) if "CMIP" not in dataset]
    style = {
        dataset: get_style(dataset, model_chars, mmm, refs) for dataset in data
    }
    basename = "cycle_" + var + "_" + timespan + "_" + cfg["region"]
    plot_path = get_plot_filename(basename, cfg)
    netcdf_path = get_diagnostic_filename(basename, cfg)

    fig, axx = plt.subplots(figsize=(10, 7))
    for dataset in data:
        cube = data[dataset]["mean"]
        axx.plot(
            cube.coord("month_number").points,
            cube.data,
            color=style[dataset]["color"],
            linestyle=style[dataset]["dash"],
            label=style[dataset]["label"],
            linewidth=style[dataset]["thick"],
        )
        if dataset not in refs:
            axx.fill_between(
                cube.coord("month_number").points,
                cube.data - data[dataset]["std"].data,
                cube.data + data[dataset]["std"].data,
                alpha=0.2,
                color=style[dataset]["color"],
                hatch=style[dataset]["hatch"],
            )

    axx.axhline(color="lightgrey", linestyle="--")
    axx.set_xlabel("Month", fontsize=12)
    axx.set_ylabel(
        timespan + " " + var.upper() + "  [%s]" % var_unit, fontsize=12
    )
    axx.legend(
        bbox_to_anchor=(1, 1.01), loc="upper left", fontsize=12, handlelength=4
    )
    axx.set_xlim(1, 12)
    axx.set_title(region_label(cfg["region"]), fontsize=14)
    plt.tight_layout()
    fig.savefig(plot_path, dpi=300)
    plt.close()

    provenance_record = get_provenance_record(
        "dummy_caption", ["mean", "stddev"], ["seas"], cfg["region"]
    )
    var_attrs = extract_variables(cfg)
    # data_cubes = iris.cube.CubeList(data)
    # io.save_1d_data(data_cubes, netcdf_path, 'month_number',
    #                var_attrs[list(var_attrs.keys())[0]])
    # Write provenance
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)
        provenance_logger.log(plot_path, provenance_record)


def plot_cycle(data, model_chars, cfg):
    # Prepare attributes
    timespan, var, var_unit, projects, cmip_projects, experiments = (
        get_plot_attributes(cfg)
    )
    mmm = False
    refs = [dataset for dataset in list(data.keys()) if "CMIP" not in dataset]
    style = {
        dataset: get_style(dataset, model_chars, mmm, refs) for dataset in data
    }
    for project in cmip_projects:
        for experiment in experiments:
            if any(d["exp"] == experiment for d in projects[project]):
                plot_path = get_plot_filename(
                    "cycle_%s_%s_%s_%s_%s"
                    % (var, timespan, cfg["region"], project, experiment),
                    cfg,
                )
                fig, axx = plt.subplots(figsize=(10, 7))
                # Get dataset style
                for dataset in data:
                    data_parts = str.split(dataset, "__")
                    if (
                        data_parts[1] == project
                        and data_parts[2] == experiment
                    ) or ("OBS" in data_parts[1].upper()):
                        cube = data[dataset]
                        axx.plot(
                            cube.coord("month_number").points,
                            cube.data,
                            color=style[dataset]["color"],
                            linestyle=style[dataset]["dash"],
                            label=style[dataset]["label"],
                            linewidth=style[dataset]["thick"],
                        )

                axx.axhline(color="lightgrey", linestyle="--")
                axx.set_xlabel("Month", fontsize=12)
                axx.set_ylabel(
                    timespan + " " + var.upper() + "  [%s]" % var_unit,
                    fontsize=12,
                )
                axx.legend(
                    bbox_to_anchor=(1, 1.01),
                    loc="upper left",
                    fontsize=12,
                    handlelength=4,
                )
                axx.set_xlim(1, 12)
                axx.set_title(region_label(cfg["region"]))
                plt.tight_layout()
                fig.savefig(plot_path, dpi=300)
                plt.close()

                provenance_record = get_provenance_record(
                    "dummy_caption",
                    ["mean", "stddev"],
                    ["seas"],
                    cfg["region"],
                )

                # Write provenance
                # with ProvenanceLogger(cfg) as provenance_logger:
                #    provenance_logger.log(netcdf_path, provenance_record)
                #    provenance_logger.log(plot_path, provenance_record)


def plot_errorbar(data, model_chars, cfg):
    # Prepare attributes
    timespan, var, var_unit, _, _, _ = get_plot_attributes(cfg)

    # One plot for each region
    plot_path = get_plot_filename(
        "errorbar_%s_%s_%s" % (var, timespan, cfg["region"]), cfg
    )
    fig, axx = plt.subplots(figsize=(14, 7))

    datasets = data
    datasets = dict(sorted(datasets.items(), key=lambda x: x[1]["mean"]))
    lab, lcolors, hlight = [], [], []
    for i, dataset in enumerate(list(datasets.keys())):
        data_parts = str.split(dataset, "__")
        if "OBS" in data_parts[1].upper():
            style = plot.get_dataset_style(data_parts[0], "cmip6_gier21")
            # style = plot.get_dataset_style('OBS', 'cmip6_gier21')
            msize, csize = 15, 10
            axx.axhline(y=data[dataset]["mean"], color="black", linestyle="--")
            axx.fill_between(
                np.arange(-1, len(datasets) + 1),
                data[dataset]["mean"] - data[dataset]["std"],
                data[dataset]["mean"] + data[dataset]["std"],
                color="lightgrey",
                alpha=0.5,
            )
            lcolors.append(style["color"])
            hlight.append(i)
        else:
            style = plot.get_dataset_style(
                data_parts[1] + "_" + data_parts[2], "cmip6_gier21"
            )
            msize = 8
            if "mmm" in data_parts[0]:
                csize = 10
                lcolors.append(style["color"])
                hlight.append(i)
            else:
                csize = 3
                lcolors.append("black")
            if (data_parts[0] in model_chars) and (
                model_chars[data_parts[0]]["n_cycle"]
            ):
                style["mark"] = "o"
            else:
                style["mark"] = "D"

        lab.append(data_parts[0])
        # Plot
        axx.errorbar(
            i,
            data[dataset]["mean"],
            yerr=data[dataset]["std"],
            marker=style["mark"],
            markersize=msize,
            color=style["color"],
            capsize=csize,
        )
    axx.axhline(color="darkgrey", linestyle="--")
    axx.set_ylabel(
        timespan + " " + var.upper() + "  [%s]" % var_unit, fontsize=12
    )
    axx.set_xticks(range(len(lab)), lab, fontsize=12, rotation=45, ha="right")
    # axx.set_xticklabels(lab, rotation=45, ha = "right", fontsize=10)
    [t.set_color(i) for (i, t) in zip(lcolors, axx.xaxis.get_ticklabels())]
    [axx.get_xticklabels()[t].set_weight("bold") for t in hlight]
    axx.set_xlim(-1, len(datasets))
    axx.set_title(region_label(cfg["region"]))

    # Make custom legend
    legend_elements = [
        Line2D([0], [0], color="#62B1F5", lw=4, label="CMIP5 historical"),
        Line2D([0], [0], color="#2A5072", lw=4, label="CMIP5 esm hist"),
        Line2D([0], [0], color="#EB523F", lw=4, label="CMIP6 historical"),
        Line2D([0], [0], color="#7B2E25", lw=4, label="CMIP6 esm hist"),
        Line2D(
            [0],
            [0],
            marker="*",
            markerfacecolor="black",
            label="OBS",
            color="w",
            markersize=18,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            markerfacecolor="black",
            label="Nitrogen cycle",
            color="w",
            markersize=12,
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            markerfacecolor="black",
            label="no N cycle",
            color="w",
            markersize=12,
        ),
    ]
    axx.legend(handles=legend_elements, loc="upper left", ncol=2, fontsize=12)

    plt.tight_layout()
    fig.savefig(plot_path, dpi=300)
    plt.close()


def get_plot_attributes(cfg):
    input_data = cfg["input_data"].values()
    start_years = group_metadata(input_data, "start_year")
    start_year = min(start_years.keys())
    end_years = group_metadata(input_data, "end_year")
    end_year = max(end_years.keys())
    timespan = str(start_year) + "-" + str(end_year)  #'1986-2005'
    vars = extract_variables(cfg)
    var = list(vars.keys())[0]  # 'nbp'
    if var == "lai":
        var_unit = "1"
    else:
        var_unit = "PgC/y"

    projects = group_metadata(input_data, "project")
    cmip_projects = [proj for proj in projects.keys() if "CMIP" in proj]
    exps = group_metadata(input_data, "exp")
    experiments = list(exps.keys())
    return timespan, var, var_unit, projects, cmip_projects, experiments


def calc_region_stat(cube, var):
    clim_mean = climate_statistics(cube, operator="mean")
    clim_std = climate_statistics(cube, operator="std_dev")
    trend = linear_trend(cube)
    if var != "lai":
        trend.convert_units("Pg m-2 yr-2")
    else:
        trend.convert_units("yr-1")
    trend_stderr = linear_trend_stderr(cube)
    if var != "lai":
        trend_stderr.convert_units("Pg m-2 yr-2")
    else:
        trend_stderr.convert_units("yr-1")
    return clim_mean.data, clim_std.data, trend.data, trend_stderr.data


def calc_region_cycle(cube):
    clim_cycle = climate_statistics(cube, operator="mean", period="monthly")
    return clim_cycle


def custom_std(cubelist):
    zstd_cube = deepcopy(cubelist[0])
    test = np.zeros(((len(cubelist)), len(zstd_cube.data)))
    for ii, cube in enumerate(cubelist):
        test[ii, :] = cubelist[ii].data
    meow = np.std(test, axis=0)
    zstd_cube.data = meow
    return zstd_cube


def get_provenance_record(caption, statistics, plot_type, region):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        "authors": ["gier_bettina"],
        "caption": caption,
        "domains": [region],
        "plot_type": plot_type,
        "realms": ["land"],
        "references": ["gier24bg"],
        "statistics": statistics,
        "themes": ["carbon"],
    }
    return record

    # provenance_record = get_provenance_record(
    #    attributes, ancestor_files=[input_file])


def main(cfg):
    """Run the diagnostic."""
    # Get model carbon cycle characteristics from yaml file
    base_dir = os.path.dirname(os.path.realpath(__file__))
    file = base_dir + "/carbon_model_characteristics.yml"
    with open(file) as infile:
        model_chars = yaml.safe_load(infile)

    # New dicts for plot data
    data_scatter = {}
    data_scatter_save = {}
    data_scatter_errorbar = {}
    obs_cubes = {"mean": [], "trend": []}
    data_cycle = {}
    data_cycle_pm = {}
    data_errorbar = {}

    # Make custom multi-model means by grouping N-cycle and non N-cycle models
    mmm_cubes = {}
    _, var, _, _, cmip_projects, experiments = get_plot_attributes(cfg)
    for proj in cmip_projects:
        for exp in experiments:
            str1 = "__".join(["Ncycle_mmm", proj, exp])
            mmm_cubes[str1] = []
            str2 = "__".join(["non-Ncycle_mmm", proj, exp])
            mmm_cubes[str2] = []

    input_data = cfg["input_data"].values()
    for data in input_data:
        name = data["dataset"]
        if (
            name == "HadGEM2-CC" or name == "HadGEM2-ES"
        ):  # this model has dec 2005 in the rcp files
            data["exp"] = "historical"
        long_name = name + "__" + data["project"] + "__" + data["exp"]
        logger.info("Processing %s", name)
        cube = iris.load_cube(data["filename"])
        if var != "lai":
            cube.convert_units("Pg m-2 yr-1")
        if data["project"] in cmip_projects:
            if not (("exclude_mmm" in cfg) and (name in cfg["exclude_mmm"])):
                if (name in model_chars) and (model_chars[name]["n_cycle"]):
                    str = "__".join(
                        ["Ncycle_mmm", data["project"], data["exp"]]
                    )
                else:
                    str = "__".join(
                        ["non-Ncycle_mmm", data["project"], data["exp"]]
                    )
                mmm_cubes[str].append(cube)

        # Add dataset entries to plot data dictionaries
        data_scatter[long_name] = {}
        save_longname = long_name.replace("__", " ").replace("_", r"\_")
        data_scatter_save[save_longname] = {}
        data_errorbar[long_name] = {}

        if name == "GCP":
            ann_cube = cube.collapsed(
                ["longitude", "latitude"], iris.analysis.MEAN
            )
            ann_cube.data = (
                ann_cube.data * 148300000000000.0
            )  # multiply by area
            iris.coord_categorisation.add_year(ann_cube, "time")
        else:
            ann_cube = annual_statistics(cube)

        mean, mean_std, trend, trend_std = calc_region_stat(ann_cube, var)

        if "OBS" in data["project"]:
            obs_cubes["mean"].append(mean)
            if name not in [
                "MTE",
                "FLUXCOM",
            ]:  # Supposed to leave these out since trend non-existent
                # TODO: turn this into recipe var 'no trend' for modularity
                obs_cubes["trend"].append(trend)

        data_scatter[long_name]["mean"] = mean
        data_scatter[long_name]["trend"] = trend
        data_scatter_save[save_longname]["mean"] = mean
        data_scatter_save[save_longname]["trend"] = trend
        data_scatter_save[save_longname]["mean std"] = mean_std
        data_scatter_save[save_longname]["trend stderr"] = trend_std

        data_errorbar[long_name]["mean"] = mean
        data_errorbar[long_name]["std"] = mean_std

        if name != "GCP":
            cycle = calc_region_cycle(cube)
            data_cycle[long_name] = cycle
            if "CMIP" not in long_name:
                data_cycle_pm[long_name] = {"mean": cycle}

    # print(obs_cubes)
    data_scatter_errorbar["OBS"] = {
        "mean": [np.mean(obs_cubes["mean"]), np.std(obs_cubes["mean"])],
        "trend": [np.mean(obs_cubes["trend"]), np.std(obs_cubes["trend"])],
    }
    # print(data_scatter_errorbar)

    for long_name in mmm_cubes:
        # print(long_name)
        if len(mmm_cubes[long_name]) > 1:
            cube_s = multi_model_statistics(
                products=iris.cube.CubeList(mmm_cubes[long_name]),
                span="overlap",
                statistics=["mean", "std_dev"],
                keep_input_datasets=False,
            )
            cube = cube_s["mean"]

            # Add dataset entries to plot data dictionaries
            save_longname = long_name.replace("__", " ").replace("_", r"\_")
            data_scatter[long_name] = {}
            data_scatter_save[save_longname] = {}
            data_errorbar[long_name] = {}

            ann_cube = annual_statistics(cube)
            mean, mean_std, trend, trend_std = calc_region_stat(ann_cube, var)

            data_scatter[long_name]["mean"] = mean
            data_scatter[long_name]["trend"] = trend

            data_scatter_save[save_longname]["mean"] = mean
            data_scatter_save[save_longname]["trend"] = trend
            data_scatter_save[save_longname]["mean std"] = mean_std
            data_scatter_save[save_longname]["trend stderr"] = trend_std

            data_errorbar[long_name]["mean"] = mean
            data_errorbar[long_name]["std"] = mean_std

            cycle = calc_region_cycle(cube)
            data_cycle[long_name] = cycle

            new_list = []
            # custom std needed because of LAI? would go from std of 0.94 to 0.05
            for cubee in mmm_cubes[long_name]:
                new_list.append(calc_region_cycle(cubee))
            data_cycle_pm[long_name] = {
                "mean": cycle,
                "std": custom_std(new_list),
            }

    # Plot
    plot_scatter(data_scatter, data_scatter_errorbar, model_chars, cfg)

    plot_cycle(data_cycle, model_chars, cfg)
    plot_project_cycle(data_cycle_pm, model_chars, cfg)

    plot_errorbar(data_errorbar, model_chars, cfg)

    # save data
    # do cycle later
    path_scatter = get_diagnostic_filename(
        "data_scatter_" + cfg["region"], cfg, extension="csv"
    )
    path_scatter_latex = get_diagnostic_filename(
        var + "_data_" + cfg["region"], cfg, extension="tex"
    )
    path_cycle_pm = get_diagnostic_filename("data_cycle_" + cfg["region"], cfg)
    path_errorbar = get_diagnostic_filename(
        "data_errorbar_" + cfg["region"], cfg, extension="csv"
    )
    # logger.info(data_scatter[region])

    dfs = pd.DataFrame.from_dict(
        data_scatter_save,
        orient="index",
        columns=["mean", "mean std", "trend", "trend stderr"],
        dtype="float",
    )
    dfs.index.name = "Models"
    # print(dfs)
    dfs.to_csv(path_scatter)
    dfs.to_latex(
        path_scatter_latex,
        float_format="%.4f",
        longtable=True,
        caption=var + " " + cfg["region"],
        label="tab:" + var + "_" + cfg["region"],
    )

    # dfs.sort_index()
    # dfs_cube = iris.pandas.as_cubes(dfs, axis=)
    # print(dfs_cube)

    # cycle: -> cubelist, each dataset gets mean, std where possible as vars?
    # add datasetname attribute + statistics std/mean
    cycle_list = []
    for dataset in data_cycle_pm:
        mean = data_cycle_pm[dataset]["mean"]
        mean.attributes["dataset"] = dataset
        mean.attributes["statistics"] = "mean"
        cycle_list.append(mean)
        if "std" in data_cycle_pm[dataset]:
            std = data_cycle_pm[dataset]["std"]
            std.attributes["dataset"] = dataset
            std.attributes["statistics"] = "std"
            cycle_list.append(std)
    cube_list = iris.cube.CubeList(cycle_list)
    # I should probably do save_1d_data instead, but var_attrs are diff?
    iris_save(cube_list, path_cycle_pm)

    dfe = pd.DataFrame.from_dict(
        data_errorbar, orient="index", columns=["mean", "std"]
    )
    # print(dfe)
    dfe.to_csv(path_errorbar)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
