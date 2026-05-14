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
from copy import deepcopy

import dask.array as da
import esmvalcore.preprocessor
import iris
import iris.plot as iplt
import iris.quickplot
import matplotlib.pyplot as plt
import numpy as np
import shapely.vectorized as shp_vect
import yaml
from cartopy.io import shapereader
from esmvalcore.preprocessor import (
    area_statistics,
    climate_statistics,
    mask_landsea,
    multi_model_statistics,
    zonal_statistics,
)

# import iris.plot as iplt
from esmvaltool.diag_scripts.shared import (
    extract_variables,
    get_plot_filename,
    group_metadata,
    plot,
    run_diagnostic,
    select_metadata,
)

logger = logging.getLogger(os.path.basename(__file__))


def _get_ne_land_mask_cube(n_lats=45, n_lons=90):
    """Get Natural Earth land mask."""
    ne_dir = os.path.join(
        os.path.dirname(os.path.realpath(esmvalcore.preprocessor.__file__)),
        "ne_masks",
    )
    ne_file = os.path.join(ne_dir, "ne_10m_land.shp")
    reader = shapereader.Reader(ne_file)
    geometries = list(reader.geometries())

    # Setup grid
    lat_coord = iris.coords.DimCoord(
        np.linspace(-90.0, 90.0, n_lats),
        var_name="lat",
        standard_name="latitude",
        long_name="latitude",
        units="degrees",
    )
    lon_coord = iris.coords.DimCoord(
        np.linspace(-180, 180.0, n_lons),
        var_name="lon",
        standard_name="longitude",
        long_name="longitude",
        units="degrees",
    )
    (lats, lons) = np.meshgrid(lat_coord.points, lon_coord.points)

    # Setup mask (1: land, 0: sea)
    mask = np.full(lats.shape, False, dtype=bool)
    for geometry in geometries:
        mask |= shp_vect.contains(geometry, lons, lats)
    land_mask = np.swapaxes(np.where(mask, 1, 0), 0, 1)

    # Setup cube
    cube = iris.cube.Cube(
        land_mask,
        var_name="land_mask",
        long_name="Land mask (1: land, 0: sea)",
        units="no_unit",
        dim_coords_and_dims=[(lat_coord, 0), (lon_coord, 1)],
    )

    return cube


def get_var_plot_atts(var, data):
    if var == "cSoil":
        varname = "Soil Carbon"
        crange_min = 0  # np.min([np.min(data[dataset].data) for dataset in data]) #5e-11
        crange_max = (
            5  # np.max([np.max(data[dataset].data) for dataset in data])
        )
        crange_anom = 6  # 4e-11#6
    elif var == "cVeg":
        varname = "Vegetation Carbon"
        crange_min = 0  # np.min([np.min(data[dataset].data) for dataset in data]) #5e-11
        crange_max = (
            5  # np.max([np.max(data[dataset].data) for dataset in data])
        )
        crange_anom = 6  # 2e-11 #6
    elif var == "nbp":
        varname = "NBP"
        crange_min = -0.3
        crange_max = (-1) * crange_min
        crange_anom = 0.2  # 6e-2 #5e-13*10e-15
    elif var == "gpp":
        varname = "GPP"
        crange_min = 0.0
        crange_max = 3.6
        crange_anom = 2.3
    elif var == "lai":
        varname = "Leaf Area Index"
        crange_min = 0.0
        crange_max = 6.0
        crange_anom = 9.0
    else:
        varname = var
        crange_min = np.min(
            [np.min(np.abs(data[dataset].data)) for dataset in data]
        )  # 5e-11
        crange_max = np.max(
            [np.max(np.abs(data[dataset].data)) for dataset in data]
        )
    return varname, crange_min, crange_max, crange_anom


def precision_round(number, digits=1):
    try:
        power = f"{number:e}".split("e")[1]
        rounded = round(number, -(int(power) - digits))
    except:
        rounded = []
        for num in number:
            power = f"{num:e}".split("e")[1]
            rounded.append(round(num, -(int(power) - digits)))
    return rounded


def without_keys(d, keys):
    return {x: d[x] for x in d if x not in keys}


def calc_hatching(data, reference, margin):
    hatching = deepcopy(reference)
    hatching.data = np.where(
        np.logical_or(
            (np.sign(data.data) == np.sign(reference.data)),
            (np.abs(data.data) + np.abs(reference.data) < margin),
        ),
        1,
        0,
    )
    mask = np.logical_or(
        np.logical_or((hatching.data == 0), reference.data.mask),
        data.data.mask,
    )
    hatching.data = da.ma.masked_array(hatching.data, mask=mask)
    return hatching


def combine_hatching(h1, h2):
    hatching = deepcopy(h1)
    hatching.data = np.where(h1.data == h2.data, h1.data, 0)
    mask = np.logical_or(
        (hatching.data == 0), np.logical_or(h1.data.mask, h2.data.mask)
    )
    hatching.data = da.ma.masked_array(hatching.data, mask=mask)
    return hatching


def calc_std_hatching(data, reference):
    hatching = deepcopy(reference)
    hatching.data = np.abs(data["mean"].data - reference.data)
    # hatching.data = np.where((data['mean'].data - data['std'].data - reference.data) < 0, 1, 0)
    hatching.data = np.where(hatching.data < data["std"].data, 1, 0)
    # hatching.data = np.where(np.logical_and((data['mean'].data-data['std'].data - reference.data) < 0,
    #                                        (data['mean'].data+data['std'].data) >= reference.data),
    #                         1, 0)
    # mask = np.logical_or(data['mean'].data.mask, reference.data.mask)
    mask = np.logical_or(hatching.data == 0, reference.data.mask)
    # np.where(hatching.data == 0, True, False)
    hatching.data = da.ma.masked_array(hatching.data, mask=mask)
    return hatching


def plot_maps(ref, data, plot_path, var, unit, cfg):
    # Calculate the number of rows and columns for the panel plot
    ncol = int((len(data) + len(ref)) / 2 + 0.5)
    nrow = 2
    fontsize = 15

    fig, axx = plt.subplots(figsize=(10, 12))
    plt.gcf().subplots_adjust(
        hspace=0.05,
        wspace=0.05,
        top=0.95,
        bottom=0.05,
        left=0.075,
        right=0.925,
    )

    varname, crange_min, crange_max, crange_anom = get_var_plot_atts(var, data)

    # Add reference
    cmap = "viridis"
    for ii, ref_dataset in enumerate(list(ref.keys())):
        plt.subplot(ncol, nrow, ii + 1)
        mesh = iplt.contourf(
            ref[ref_dataset],
            np.linspace(crange_min, crange_max, 13),
            cmap=cmap,
            extend="both",
        )
        plt.gca().coastlines()
        # plt.title(ref_dataset)
        plt.title(get_label(ref_dataset), loc="left", fontsize=fontsize)
        if var == "lai":
            area_operator = "mean"
        else:
            area_operator = "sum"
        map_mean = area_statistics(ref[ref_dataset], area_operator)
        plt.title("%.2g" % map_mean.data, loc="right", fontsize=fontsize)

    for ii, dataset in enumerate(list(data.keys())):
        if len(ref) > 1:
            ref_mean = ref["Reference Mean"]
        else:
            ref_mean = ref[list(ref.keys())[0]]
        hatch = calc_std_hatching(data[dataset], ref_mean)
        # Calc Anomaly
        anomaly = data[dataset]["mean"] - ref_mean
        if cfg["relative_difference"]:
            abs_cube = ref_mean
            abs_cube.data = np.abs(abs_cube.data)
            anomaly = anomaly / abs_cube

        plt.subplot(
            ncol, nrow, len(ref) + 1 + len(ref) % nrow + ii
        )  # (nrow - len(ref)%nrow)
        cf = iplt.contourf(
            anomaly,
            np.linspace(-crange_anom, crange_anom, 13),
            cmap="coolwarm",
            extend="both",
        )  # ,
        # vmin = -crange, vmax = crange)#, contour_levels)
        iplt.contourf(
            hatch,
            colors="none",
            # hatches=['+++', None])
            hatches=[r" \ \ \ \ \ ", None],
        )
        # Add coastlines
        plt.gca().coastlines()
        # Add year as title
        plt.title(get_label(dataset), loc="left", fontsize=fontsize)

        map_mean = area_statistics(anomaly, "mean")
        #    anomaly.collapsed(['longitude', 'latitude'], iris.analysis.MEAN))
        # map_mean=np.mean(anomaly.data)
        plt.title("%.2g" % map_mean.data, loc="right", fontsize=fontsize)

    plt.tight_layout()
    # Make room for a combined colorbars
    plt.subplots_adjust(bottom=0.15)
    plt.subplots_adjust(left=0.08)
    colorbar_axes = plt.gcf().add_axes([0.35, 0.1, 0.4, 0.025])
    colorbar = plt.colorbar(cf, colorbar_axes, orientation="horizontal")  # ,
    # ticks = np.linspace(-crange_anom, crange_anom, 13))
    if cfg["relative_difference"]:
        colorbar.set_label(
            "Relative " + varname + " Bias [1]", fontsize=fontsize
        )
    else:
        colorbar.set_label(varname + " Bias [" + unit + "]", fontsize=fontsize)
    colorbar.ax.tick_params(labelsize=12)
    colorbar_axes1 = plt.gcf().add_axes([0.5, 0.827, 0.015, 0.145])
    colorbar1 = plt.colorbar(mesh, colorbar_axes1, orientation="vertical")
    colorbar1.set_label(varname + " [" + unit + "]", fontsize=fontsize)

    colorbar1.ax.tick_params(labelsize=12)

    fig.savefig(plot_path, dpi=300)
    plt.close()


def plot_maps_seperate(ref, data, var, unit, cfg):
    ncol = int(len(ref) / 2 + 0.5)  # int(len(ref)/2 + 1)
    nrow = 2
    plot_path = get_plot_filename(var + "_maps_ref", cfg)
    fig, axx = plt.subplots(figsize=(10, 3 * ncol))
    fontsize = 13
    avg_rnd = 4
    # Calculate the number of rows and columns for the panel plot

    # nrow = int(np.ceil(np.sqrt(len(data)+1)))
    # ncol = int(np.ceil((len(data)+1) / nrow))
    plt.gcf().subplots_adjust(
        hspace=0.05,
        wspace=0.05,
        top=0.95,
        bottom=0.05,
        left=0.075,
        right=0.925,
    )

    varname, crange_min, crange_max, crange_anom = get_var_plot_atts(var, data)
    if var in ["gpp", "nbp"]:
        zunit = "PgC yr-1"
    elif var == "lai":
        zunit = ""
    else:
        zunit = "PgC"

    # Calculate hatching: Same sign or both values close to zero
    margin = 0.5e-13
    # (np.max(ref["Reference Mean"].data) - np.min(ref["Reference Mean"].data)) / 5.
    if len(ref) > 1:
        ref_names = list(
            (without_keys(ref, ["Reference Mean", "Reference Range"])).keys()
        )
        hatching = calc_hatching(ref[ref_names[0]], ref[ref_names[1]], margin)
        if len(ref_names) == 2:
            ref_anomaly = ref[ref_names[0]] - ref[ref_names[1]]

    # Add reference
    if len(ref) > 1:
        if var == "nbp":
            cmap = "BrBG"
        else:
            cmap = "viridis"
        for ii, ref_dataset in enumerate(list(ref.keys())):
            if ref_dataset == "Reference Range":
                continue
            elif ref_dataset == "Reference Mean":
                plt.subplot(ncol, nrow, len(ref) - 1 + len(ref) % nrow)
            else:
                plt.subplot(ncol, nrow, ii + 1)
            mesh = iplt.contourf(
                ref[ref_dataset],
                np.linspace(crange_min, crange_max, 13),
                cmap=cmap,
            )
            if var == "nbp":
                iplt.contourf(hatching, colors="none", hatches=["/////", None])
            plt.gca().coastlines()
            plt.title(ref_dataset, loc="left", fontsize=fontsize)

            area_operator = "mean"
            # if var == 'lai':
            #    area_operator = 'mean'
            # else:
            #    area_operator = 'sum'
            map_mean = area_statistics(ref[ref_dataset], area_operator)
            plt.title(
                np.round(map_mean.data, avg_rnd),
                loc="right",
                fontsize=fontsize,
            )
            # plt.title('%.2e'%map_mean.data, loc='right', fontsize=fontsize)
        if var == "nbp":
            ref_range = [0, 0.3]
        elif var == "gpp":
            ref_range = [0, 0.9]
        elif var == "lai":
            ref_range = [0, 1.8]
        plt.subplot(ncol, nrow, len(ref) + len(ref) % nrow)  # 1+
        anom = iplt.contourf(
            ref["Reference Range"],
            np.linspace(ref_range[0], ref_range[1], 13),
            cmap="plasma",
            extend="max",
        )
        plt.gca().coastlines()
        plt.title("Reference Range", loc="left", fontsize=fontsize)
        # map_mean = ref[ref_dataset].collapsed(['longitude', 'latitude'], iris.analysis.MEAN)
        map_mean = area_statistics(ref[ref_dataset], "mean")
        plt.title(
            np.round(map_mean.data, avg_rnd), loc="right", fontsize=fontsize
        )
        # plt.title('%.2g'%map_mean.data, loc='right', fontsize=fontsize)
    else:
        for ii, ref_dataset in enumerate(list(ref.keys())):
            plt.subplot(ncol, nrow, ii + 1)
            mesh = iplt.contourf(
                ref[ref_dataset],
                np.linspace(crange_min, crange_max, 13),
                cmap="viridis",
            )
            plt.gca().coastlines()
            plt.title(ref_dataset, loc="left", fontsize=fontsize)
            area_operator = "mean"
            # if var == 'lai':
            #    area_operator = 'mean'
            # else:
            #    area_operator = 'sum'
            map_mean = area_statistics(ref[ref_dataset], area_operator)
            plt.title(
                np.round(map_mean.data, avg_rnd),
                loc="right",
                fontsize=fontsize,
            )
            # plt.title('%.2g'%map_mean.data, loc='right', fontsize=fontsize)

    plt.subplots_adjust(bottom=0.15)
    plt.subplots_adjust(left=0.08)

    # if len(ref) == 3:
    if len(ref) > 1:
        colorbar_axes = plt.gcf().add_axes([0.079, 0.1, 0.413, 0.025])
        colorbar_axes_anom = plt.gcf().add_axes([0.513, 0.1, 0.413, 0.025])
        colorbar_anom = plt.colorbar(
            anom, colorbar_axes_anom, orientation="horizontal"
        )

        colorbar_anom.ax.tick_params(labelsize=12)
        colorbar_anom.set_label(
            "Reference Range" + " [" + unit + "]", fontsize=fontsize
        )
    else:
        colorbar_axes = plt.gcf().add_axes([0.35, 0.1, 0.4, 0.025])

    colorbar = plt.colorbar(
        mesh,
        colorbar_axes,
        orientation="horizontal",
        ticks=np.linspace(crange_min, crange_max, 7),
    )
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_label(varname + " [" + unit + "]", fontsize=fontsize)

    fig.savefig(plot_path, dpi=300)
    plt.close()

    plot_path = get_plot_filename(var + "_maps_means", cfg)
    fig, axx = plt.subplots(figsize=(8, 10))
    # Calculate the number of rows and columns for the panel plot
    ncol = int(len(data.keys()) / 2 + 1)
    nrow = 2

    for ii, dataset in enumerate(list(data.keys())):
        if len(ref) > 1:
            ref_mean = ref["Reference Mean"]
        else:
            # hatch = calc_hatching(data[dataset], ref[list(ref.keys())[0]], margin)
            ref_mean = ref[list(ref.keys())[0]]
        hatch = calc_std_hatching(data[dataset], ref_mean)
        # Calc Anomaly
        anomaly = data[dataset]["mean"] - ref_mean
        if cfg["relative_difference"]:
            abs_cube = ref_mean
            abs_cube.data = np.abs(abs_cube.data)
            anomaly = anomaly / abs_cube
        plt.subplot(ncol, nrow, ii + 1)  # 2 if project_mean = False
        cf = iplt.contourf(
            anomaly,
            precision_round(np.linspace(-crange_anom, crange_anom, 13)),
            cmap="coolwarm",
            extend="both",
        )  # ,
        # vmin = -crange, vmax = crange)#, contour_levels)
        iplt.contourf(hatch, colors="none", hatches=[r" \ \ \ \ \ ", None])
        if var == "nbp":
            iplt.contourf(hatching, colors="none", hatches=["/////", None])
        # mesh = iplt.pcolormesh(data[dataset], cmap='coolwarm')
        # Add coastlines
        plt.gca().coastlines()
        # Add datasetname as title
        plt.title(get_label(dataset), loc="left", fontsize=fontsize)

        # map_mean = anomaly.collapsed(['longitude', 'latitude'], iris.analysis.MEAN)
        map_mean = area_statistics(anomaly, "mean")
        # map_mean = np.mean(anomaly.data)
        plt.title(
            np.round(map_mean.data, avg_rnd), loc="right", fontsize=fontsize
        )
        # plt.title('%.2e'%map_mean.data, loc='right', fontsize=fontsize)

    plt.tight_layout()
    # Make room for a combined colorbars
    colorbar_axes = plt.gcf().add_axes([0.35, 0.15, 0.4, 0.025])
    colorbar = plt.colorbar(cf, colorbar_axes, orientation="horizontal")
    colorbar.ax.tick_params(labelsize=12)
    # ,
    # ticks = np.linspace(-crange_anom, crange_anom, 7))
    if cfg["relative_difference"]:
        colorbar.set_label(
            "Relative " + varname + " Bias [" + unit + "]", fontsize=fontsize
        )
    else:
        colorbar.set_label(varname + " Bias [" + unit + "]", fontsize=fontsize)

    fig.savefig(plot_path, dpi=300)
    plt.close()


def plot_models(ref, data, var, unit, cfg, model_chars):
    for projs in data:
        if len(data[projs]) > 0:
            plot_path = get_plot_filename(
                var + "_" + projs + "_maps_models", cfg
            )
            ncol = int(np.ceil(np.sqrt(len(data[projs])) - 0.5))
            nrow = int(np.ceil((len(data[projs]) + 1) / ncol))

            fig, axx = plt.subplots(figsize=(5 * ncol, 3 * nrow))
            # Calculate the number of rows and columns for the panel plot

            varname, crange_min, crange_max, crange_anom = get_var_plot_atts(
                var, data[projs]
            )

            if var == "nbp":
                cmap = "BrBG"
            elif var == "cSoil":
                cmap = "viridis"
                crange_min = 0
                crange_max = 10e-11
            else:
                cmap = "viridis"
            for ii, dataset in enumerate(list(data[projs].keys())):
                plt.subplot(nrow, ncol, ii + 1)
                cf = iplt.contourf(
                    data[projs][dataset],
                    np.linspace(crange_min, crange_max, 13),
                    cmap=cmap,
                    extend="both",
                )
                # Add coastlines
                plt.gca().coastlines()
                # Add name as title
                if dataset in model_chars:
                    if model_chars[dataset]["n_cycle"]:
                        color = "green"
                    else:
                        color = "black"
                else:
                    color = "black"
                # plt.title(dataset, color = color)
                plt.title(get_label(dataset), loc="left", color=color)

                if var == "lai":
                    area_operator = "mean"
                else:
                    area_operator = "sum"
                map_mean = area_statistics(data[projs][dataset], area_operator)
                plt.title("%.2g" % map_mean.data, loc="right")

            plt.tight_layout()
            # Make room for a combined colorbars
            plt.subplots_adjust(bottom=0.15)
            plt.subplots_adjust(left=0.08)
            colorbar_axes = plt.gcf().add_axes([0.35, 0.1, 0.4, 0.025])
            colorbar = plt.colorbar(
                cf, colorbar_axes, orientation="horizontal"
            )  # ,
            # ticks = np.linspace(crange_min, crange_max, 13))
            colorbar.set_label(varname + " [" + unit + "]")

            fig.savefig(plot_path, dpi=300)
            plt.close()


def plot_zmeans(data, plot_path, var, unit, style, stds=[0]):
    varname = var
    fig, axx = plt.subplots(figsize=(10, 6))
    for dataset in data:
        cube = data[dataset]
        axx.plot(
            cube.coord("latitude").points,
            cube.data,
            label=style[dataset]["label"],
            color=style[dataset]["color"],
            linestyle=style[dataset]["dash"],
            linewidth=style[dataset]["thick"],
        )
        if dataset in stds:
            axx.fill_between(
                cube.coord("latitude").points,
                cube.data - stds[dataset].data,
                cube.data + stds[dataset].data,
                alpha=0.2,
                color=style[dataset]["color"],
                hatch=style[dataset]["hatch"],
            )
    axx.axhline(color="lightgrey", linestyle="--")
    axx.set_xlabel("Latitude [°]", fontsize=12)
    if var == "cSoil":
        axx.set_xlim(30, 83)
    else:
        axx.set_xlim(-50, 83)
    axx.set_ylabel(varname.upper() + " [" + unit + "]", fontsize=12)
    axx.legend(
        bbox_to_anchor=(-0.10, -0.09),
        loc="upper left",
        ncol=3,
        handlelength=3.5,
        fontsize=12,
    )
    #     axx.legend(loc="best", ncol = 2, handlelength = 4)
    plt.tight_layout()
    fig.savefig(plot_path, dpi=300)
    plt.close()


def plot_zonal_means(data, refs, var, unit, cfg, model_chars, stds=[0]):
    """Plot zonal means of data"""
    if type(refs) == dict:
        for proj in data:
            new_data = {**data[proj], **refs}
            plot_path = get_plot_filename(var + "_zonal_means_" + proj, cfg)
            style = {
                dataset: get_style(dataset, model_chars, refs)
                for dataset in new_data
            }
            plot_zmeans(new_data, plot_path, var, unit, style)
    else:
        plot_path = get_plot_filename(var + "_zonal_means", cfg)
        style = {
            dataset: get_style(dataset, model_chars, refs) for dataset in data
        }
        plot_zmeans(data, plot_path, var, unit, style, stds)


def get_label(datasetname):
    data_parts = str.split(datasetname, "_")
    if len(data_parts) != 4:
        label = datasetname
    elif "CMIP" in data_parts[-2]:
        experiment = data_parts.pop()
        if experiment == "historical":
            data_parts[-1] = data_parts[-1] + "c"
        elif experiment in ["esmHistorical", "esm-hist"]:
            data_parts[-1] = data_parts[-1] + "e"
        else:
            data_parts[-1] = data_parts[-1] + experiment
        label = data_parts[-1] + " " + data_parts[0] + " " + data_parts[1]
    else:
        label = datasetname

    return label


def get_style(dataset, model_chars, refs=[]):
    data_parts = str.split(dataset, "_")
    if dataset in refs:
        style_name = dataset
        lab = dataset
    else:
        style_name = "_".join(data_parts[-2:])
        lab = " ".join(data_parts)
    style = plot.get_dataset_style(style_name, "cmip6_gier21")
    style["label"] = get_label(dataset)
    # style['label'] = lab
    if dataset not in refs:
        if model_chars[data_parts[0]]["n_cycle"]:
            style["thick"] = 3
            style["dash"] = "--"
            style["hatch"] = "----"
        else:
            style["thick"] = 2
            style["hatch"] = "||||"
    else:
        style["hatch"] = None
    return style


def mark_reference_datasets(cfg, input_data):
    """Mark reference datasets."""
    input_data = deepcopy(input_data)
    for dataset in input_data:
        if dataset["dataset"] in cfg["reference_datasets"]:
            dataset["ref"] = True
        else:
            dataset["ref"] = False
    return input_data


def custom_std(cubelist):
    zstd_cube = deepcopy(cubelist[0])
    test = np.zeros(((len(cubelist)), len(zstd_cube.data)))
    for ii, cube in enumerate(cubelist):
        test[ii, :] = cubelist[ii].data
    meow = np.std(test, axis=0)
    zstd_cube.data = meow
    return zstd_cube


def main(cfg):
    """Run the diagnostic."""
    # Get model carbon cycle characteristics from yaml file
    base_dir = os.path.dirname(os.path.realpath(__file__))
    file = base_dir + "/carbon_model_characteristics.yml"
    with open(file) as infile:
        model_chars = yaml.safe_load(infile)

    # Read variable name
    vars = extract_variables(cfg)
    var = list(vars.keys())[0]

    input_data = list(cfg["input_data"].values())
    input_data = mark_reference_datasets(cfg, input_data)

    ref = list(group_metadata(input_data, "reference_dataset").keys())[0]
    collect_data = {}
    collect_data_model = {}
    zmeans = {}
    zmeans_ref = {}
    zstds = {}
    zmeans_models = {}

    if var == "lai":
        zoperator = "mean"
    else:
        zoperator = "sum"

    project_mean = True
    if project_mean == True:
        projects = group_metadata(input_data, "project")
        cmip_projects = [proj for proj in projects.keys() if "CMIP" in proj]
        exps = group_metadata(input_data, "exp")
        experiments = list(exps.keys())
        mmm_cubes = {}
        mmm_names = {}
        for proj in cmip_projects:
            meow = projects[proj]
            proj_exp = group_metadata(meow, "exp")
            for exp in experiments:
                if exp in proj_exp:
                    str1 = "_".join(["Ncycle_MMM", proj, exp])
                    mmm_cubes[str1] = []
                    mmm_names[str1] = []
                    str2 = "_".join(["non-Ncycle_MMM", proj, exp])
                    mmm_cubes[str2] = []
                    mmm_names[str2] = []
                    collect_data_model["_".join([proj, exp])] = {}
                    zmeans_models["_".join([proj, exp])] = {}

    # Load reference data first
    refs = {}
    ref_datasets = select_metadata(input_data, ref=True)

    # Prepare knowledge of allowed missing values
    cfg.setdefault("missing_threshold", 0.85)
    land_mask = _get_ne_land_mask_cube(n_lats=90, n_lons=180)
    land_zmean = zonal_statistics(land_mask, operator="sum")

    for dataset in ref_datasets:
        ref_name = dataset["dataset"]
        ref_cube = iris.load_cube(dataset["filename"])
        ref_cube.convert_units(dataset["plot_units"].replace("C", ""))
        # ref_cube = mask_landsea(ref_cube, 'sea', threshold=100.)
        if ref_cube.ndim > 2:
            ref_cube = climate_statistics(ref_cube, operator="mean")
            ref_cube.remove_coord("time")
            try:
                ref_cube.remove_coord("month_number")
            except:
                pass
            if var == "gpp":
                ref_cube.standard_name = (
                    "gross_primary_productivity_of_biomass_expressed_as_carbon"
                )
                ref_cube.long_name = "Carbon Mass Flux out of Atmosphere Due to Gross Primary Production on Land [kgC m-2 s-1]"
        refs[ref_name] = ref_cube
        climate_mean = deepcopy(ref_cube)
        if zoperator == "sum":
            weights = iris.analysis.cartography.area_weights(climate_mean)
            climate_mean.data *= weights.data
            if var in ["gpp", "nbp"]:
                new_unit = "Pg m-2 yr-1"
            else:
                new_unit = "Pg m-2"
            climate_mean.convert_units(new_unit)
        zmeans_ref[ref_name] = zonal_statistics(climate_mean, zoperator)

        # Do not display zonal mean for ref if too many missing values per latitude
        ref_cube_mv = deepcopy(ref_cube)
        ref_cube_mv.data = np.where(ref_cube.data.mask, 0, 1)
        land_zmean_ref = zonal_statistics(ref_cube_mv, operator="sum")

        zmeans_ref[ref_name].data.mask = np.where(
            land_zmean_ref.data < land_zmean.data * cfg["missing_threshold"],
            True,
            False,
        )

    if len(refs) > 1:
        ref_cubelist = iris.cube.CubeList(list(refs.values()))
        ref_data_mean = multi_model_statistics(
            products=ref_cubelist,
            span="overlap",
            statistics=["mean", "min", "max"],
            keep_input_datasets=False,
        )
        refs["Reference Mean"] = ref_data_mean["mean"]
        refs["Reference Range"] = ref_data_mean["max"] - ref_data_mean["min"]

    # # Get unit for conversion, but iris doesn't know PgC
    unit = select_metadata(input_data, dataset=ref)[0]["plot_units"]
    # if var=='lai':
    zunit = unit
    # elif var in ['cVeg', 'cSoil']:
    #    zunit = "PgC"
    # elif var in ['gpp', 'nbp']:
    #    zunit = "PgC yr-1"

    # refs= {ref: ref_data_mean}
    # Read data and calculate anomaly
    for data in input_data:
        name = data["dataset"]
        if (
            name == "HadGEM2-CC" or name == "HadGEM2-ES"
        ):  # this model has dec 2005 in the rcp files
            data["exp"] = "historical"
        if not data["ref"]:
            cube_data = iris.load_cube(data["filename"])
            cube_data.convert_units(data["plot_units"].replace("C", ""))
            # Regular landsea masking doesn't work.. let's combine it with a search for a value of 0!

            mask_data = deepcopy(cube_data)
            mask_data = mask_landsea(mask_data, "sea", threshold=100.0)
            cube_data.data.mask = np.where(
                np.logical_and(
                    cube_data.data == 0, mask_data.data.mask == True
                ),
                True,
                cube_data.data.mask,
            )

            # if not name=='SAM0-UNICON': #that model does some weird shit- TODO! but both mask and data seem ok
            # cube_data = mask_landsea(cube_data, 'sea', threshold=100.)
            if project_mean == True:
                if data["project"] in cmip_projects:
                    if not (
                        ("exclude_mmm" in cfg) and (name in cfg["exclude_mmm"])
                    ):
                        if (name in model_chars) and (
                            model_chars[name]["n_cycle"]
                        ):
                            str = "_".join(
                                ["Ncycle_MMM", data["project"], data["exp"]]
                            )
                        else:
                            str = "_".join(
                                [
                                    "non-Ncycle_MMM",
                                    data["project"],
                                    data["exp"],
                                ]
                            )
                        mmm_cubes[str].append(cube_data)
                        mmm_names[str].append(name)
            proj, exp = data["project"], data["exp"]
            clim_mean = climate_statistics(cube_data, operator="mean")
            clim_mean.remove_coord("time")
            collect_data_model["_".join([proj, exp])][name] = (
                clim_mean  # - ref_data_mean
            )
            climate_mean = deepcopy(clim_mean)
            if zoperator == "sum":
                weights = iris.analysis.cartography.area_weights(climate_mean)
                climate_mean.data *= weights.data
                if var in ["gpp", "nbp"]:
                    new_unit = "Pg m-2 yr-1"
                    zunit = "PgC yr-1"
                else:
                    new_unit = "Pg m-2"
                    zunit = "PgC"
                climate_mean.convert_units(new_unit)
            zmeans_models["_".join([proj, exp])][name] = zonal_statistics(
                climate_mean, zoperator
            )

    if project_mean == True:
        for long_name in mmm_cubes:
            if len(mmm_cubes[long_name]) > 1:
                cube_s = multi_model_statistics(
                    products=iris.cube.CubeList(mmm_cubes[long_name]),
                    span="overlap",
                    statistics=["mean", "std_dev"],
                    keep_input_datasets=False,
                )
                clim_mean = climate_statistics(cube_s["mean"], operator="mean")
                std_mean = climate_statistics(
                    cube_s["std_dev"], operator="mean"
                )
                collect_data[long_name] = {}
                collect_data[long_name]["mean"] = clim_mean  # col_data
                collect_data[long_name]["std"] = std_mean
                collect_data_model["_".join(long_name.split("_")[2:])][
                    long_name
                ] = clim_mean
                climate_mean = deepcopy(clim_mean)
                if zoperator == "sum":
                    weights = iris.analysis.cartography.area_weights(
                        climate_mean
                    )
                    climate_mean.data *= weights.data
                    if var in ["gpp", "nbp"]:
                        new_unit = "Pg m-2 yr-1"
                        zunit = "PgC yr-1"
                    else:
                        new_unit = "Pg m-2"
                        zunit = "PgC"
                    climate_mean.convert_units(new_unit)
                zmeans[long_name] = zonal_statistics(climate_mean, zoperator)
                zstd_cubelist = [
                    zmeans_models["_".join(long_name.split("_")[2:])][dataset]
                    for dataset in mmm_names[long_name]
                ]
                cube_zstd = multi_model_statistics(
                    products=iris.cube.CubeList(zstd_cubelist),
                    span="overlap",
                    statistics=["mean", "std_dev"],
                    keep_input_datasets=False,
                )
                zstds[long_name] = cube_zstd[
                    "std_dev"
                ]  # custom_std(zstd_cubelist)
                zmeans_models["_".join(long_name.split("_")[2:])][
                    long_name
                ] = zonal_statistics(climate_mean, zoperator)

    # Plot

    plt.rcParams["hatch.linewidth"] = 0.2

    plot_path = get_plot_filename(var + "_maps", cfg)
    plot_maps(refs, collect_data, plot_path, var, unit, cfg)
    plot_maps_seperate(refs, collect_data, var, unit, cfg)
    plot_models(refs, collect_data_model, var, unit, cfg, model_chars)
    plot_zonal_means(
        {**zmeans, **zmeans_ref},
        list(
            (without_keys(refs, ["Reference Mean", "Reference Range"])).keys()
        ),
        var,
        zunit,
        cfg,
        model_chars,
        zstds,
    )
    plot_zonal_means(zmeans_models, zmeans_ref, var, zunit, cfg, model_chars)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
