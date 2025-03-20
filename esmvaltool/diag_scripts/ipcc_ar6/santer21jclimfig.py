#!/usr/bin/env python
"""Diagnostic script to plot figure 3_12 of IPCC AR 6 chapter 3.

Description
-----------
Trends in CMIP5 and CMIP6 compared to observational and reanalysis data
based on methods used in Santer et al. (2007) and Santer et al. (2021).
https://www.pnas.org/doi/10.1073/pnas.0702872104
https://doi.org/10.1175/JCLI-D-20-0768.1

Author
------
Katja Weigel (IUP, Uni Bremen, and DLR, Germany)

Configuration options in recipe
-------------------------------
sample_obs: optional, filter all data sets (netCDF file with 0 and 1 for
            used grid). The data sets must be interpolated to the same lat/lon
            grid as the filter.
            The filter must cover at least the time period used for the data.

##############################################################################

"""

import logging
import os
from collections import OrderedDict
from pprint import pformat

import iris
import iris.coord_categorisation as cat
import matplotlib.pyplot as plt
import numpy as np
from cf_units import Unit
from scipy import stats

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    extract_variables,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    plot,
    run_diagnostic,
    select_metadata,
    variables_available,
)

logger = logging.getLogger(os.path.basename(__file__))


def _apply_filter(cfg, cube):
    """Apply filter from RSS Anomalies to all data and calculates mean."""
    if "sample_obs" in cfg:
        filt = iris.load(
            os.path.join(cfg["auxiliary_data_dir"], cfg["sample_obs"])
        )[0]

        for iii, dim_str in enumerate(["time", "latitude", "longitude"]):
            filt = _fit_dim(filt, cube, dim_str, iii)

        cube.data = cube.data * filt.data

    coords = ("longitude", "latitude")
    for coord in coords:
        if not cube.coord(coord).has_bounds():
            cube.coord(coord).guess_bounds()
    cube_grid_areas = iris.analysis.cartography.area_weights(cube)
    cube = cube.collapsed(coords, iris.analysis.MEAN, weights=cube_grid_areas)

    return cube


def _fit_dim(filt, cube, dim_str, dim_nr):
    """Adjust array for given dimension."""
    dimfil = filt.coord(dim_str)
    dimcu = cube.coord(dim_str)

    if dim_str == "time":
        start = dimcu.units.num2date(dimcu.points[0])
        end = dimcu.units.num2date(dimcu.points[-1])
        sdim = dimfil.nearest_neighbour_index(dimfil.units.date2num(start))
        edim = dimfil.nearest_neighbour_index(dimfil.units.date2num(end))
    else:
        start = dimcu.points[0]
        end = dimcu.points[-1]
        sdim = dimfil.nearest_neighbour_index(start)
        edim = dimfil.nearest_neighbour_index(end)

    if dim_nr == 0:
        filt = filt[sdim : edim + 1, :, :]
    elif dim_nr == 1:
        filt = filt[:, sdim : edim + 1, :]
    elif dim_nr == 2:
        filt = filt[:, :, sdim : edim + 1]

    return filt


def _calc_trend(cube_anom):
    """Calculate trends."""
    reg_var = stats.linregress(
        np.linspace(0, len(cube_anom.data) - 1, len(cube_anom.data)),
        cube_anom.data,
    )

    return reg_var.slope


def _calculate_anomalies(cube):
    """Remove annual cycle from time series."""
    c_n = cube.aggregated_by("month_number", iris.analysis.MEAN)
    month_number = cube.coord("month_number").points
    startind = np.where(month_number == 1)
    endind = np.where(month_number == 12)
    data_in = cube.data
    data_new = np.full(data_in.shape, 0.5)
    for iii in range((startind[0])[0], (endind[0])[-1], 12):
        data_new[iii : iii + 12] = (
            ((cube.data[iii : iii + 12] - c_n.data) / c_n.data) * 100.0 * 120.0
        )

    cube.data = data_new
    cube.units = Unit("percent")

    return cube


def _check_full_data(dataset_path, cube):
    """Check if cube covers time series from start year to end year."""
    cstart_year = cube.coord("year").points[0]
    cend_year = cube.coord("year").points[-1]
    start_year = int(dataset_path["start_year"])
    end_year = int(dataset_path["end_year"])

    check = 0
    if start_year == cstart_year and end_year == cend_year:
        check = 1
    else:
        logger.info(
            "Time series incomplete: "
            "Start(Recipe, Data): %i, %i; "
            "End(Recipe, Data): %i, %i"
            "Dataset %s excluded",
            start_year,
            cstart_year,
            end_year,
            cend_year,
            dataset_path,
        )

    return check


def _get_hem_letter(lat):
    """Get S or N from Latitude."""
    shem = ""
    if lat < 0:
        shem = "S"
    if lat > 0:
        shem = "N"

    return shem


def _get_plotobscol(iii):
    """Get color based on index."""
    if iii > 12:
        obscoli = float(iii) - 12.25
    elif iii > 8:
        obscoli = float(iii) - 8.75
    elif iii > 4:
        obscoli = float(iii) - 4.5
    else:
        obscoli = float(iii)

    plotobscol = (
        1.0 - 0.25 * obscoli,
        0.25 * obscoli,
        0.5 + obscoli * 0.1,
    )
    return plotobscol


def _get_valid_datasets(input_data):
    """Get valid datasets list and number for each model."""
    number_of_subdata = OrderedDict()
    available_dataset = list(group_metadata(input_data, "dataset"))
    valid_datasets = []
    period = {}
    period["start_year"] = []
    period["end_year"] = []
    period["span"] = []
    period["slat"] = []
    period["elat"] = []
    for dataset in available_dataset:
        meta = select_metadata(input_data, dataset=dataset)
        number_of_subdata[dataset] = float(len(meta))
        for dataset_path in meta:
            cube = iris.load(dataset_path["filename"])[0]
            cat.add_year(cube, "time", name="year")
            if not _check_full_data(dataset_path, cube):
                number_of_subdata[dataset] = number_of_subdata[dataset] - 1
            else:
                valid_datasets.append(dataset_path["filename"])
                period["start_year"].append(int(dataset_path["start_year"]))
                period["end_year"].append(int(dataset_path["end_year"]))
                period["span"].append(
                    int(dataset_path["end_year"])
                    - int(dataset_path["start_year"])
                    + 1
                )
                period["slat"].append(round(cube.coord("latitude").points[0]))
                period["elat"].append(round(cube.coord("latitude").points[-1]))
    if min(period["span"]) == max(period["span"]):
        period["common_span"] = str(min(period["span"]))
    if min(period["start_year"]) == max(period["start_year"]) and min(
        period["end_year"]
    ) == max(period["end_year"]):
        period["common_period"] = (
            str(min(period["start_year"]))
            + " - "
            + str(min(period["end_year"]))
        )

    if min(period["slat"]) == max(period["slat"]) and min(
        period["elat"]
    ) == max(period["elat"]):
        shem = _get_hem_letter(min(period["slat"]))
        ehem = _get_hem_letter(min(period["elat"]))

        period["common_lats"] = (
            str(abs(min(period["slat"])))
            + "°"
            + shem
            + " - "
            + str(abs(min(period["elat"])))
            + "°"
            + ehem
        )

    return valid_datasets, number_of_subdata, period


def _make_list_dict(res_ar, trends, var, print_long_name):
    """Make List diectionary."""
    list_dict = {}
    list_dict["data"] = [res_ar["xhist"]]
    list_dict["name"] = [
        {
            "var_name": var + "_trends_bins",
            "long_name": print_long_name + " Trend bins",
            "units": "percent",
        }
    ]
    if trends["cmip5"]:
        list_dict["data"].append(res_ar["artrend_c5"])
        list_dict["name"].append(
            {
                "var_name": var + "_trends_cmip5",
                "long_name": print_long_name + " Trends CMIP5",
                "units": "percent",
            }
        )
        list_dict["data"].append(res_ar["weights_c5"])
        list_dict["name"].append(
            {
                "var_name": "data_set_weights",
                "long_name": "Weights for each data set.",
                "units": "1",
            }
        )
        list_dict["data"].append(res_ar["kde1_c5"](res_ar["xhist"]))
        list_dict["name"].append(
            {
                "var_name": var + "_trend_distribution_cmip5",
                "long_name": print_long_name
                + " Trends "
                + "distribution CMIP5",
                "units": "1",
            }
        )
    if trends["cmip6"]:
        list_dict["data"].append(res_ar["artrend_c6"])
        list_dict["name"].append(
            {
                "var_name": var + "_trends_cmip6",
                "long_name": print_long_name + " Trends CMIP6",
                "units": "percent",
            }
        )
        list_dict["data"].append(res_ar["weights_c6"])
        list_dict["name"].append(
            {
                "var_name": "data_set_weights",
                "long_name": "Weights for each data set.",
                "units": "1",
            }
        )
        list_dict["data"].append(res_ar["kde1_c6"](res_ar["xhist"]))
        list_dict["name"].append(
            {
                "var_name": var + "_trend_distribution_cmip6",
                "long_name": print_long_name
                + " Trends "
                + "distribution CMIP6",
                "units": "1",
            }
        )

    if trends["obs"]:
        for obsname in trends["obs"].keys():
            list_dict["data"].append(trends["obs"][obsname])
            list_dict["name"].append(
                {
                    "var_name": var + "_trend_" + obsname,
                    "long_name": print_long_name + " Trend " + obsname,
                    "units": "percent",
                }
            )

    return list_dict


def _plot_extratrends(cfg, extratrends, trends, period, axx_lim):
    """Plot trends for ensembles."""
    res_ar = {"artrend": {}, "kde1": {}}
    res_ar["xval"] = axx_lim["xval"]
    res_ar["xhist"] = axx_lim["xhist"]
    fig, axx = plt.subplots(figsize=(8, 6))
    names = {}

    names["valid_datasets"] = []

    for xtrmdl in cfg["add_model_dist"]:
        alias = list(extratrends[xtrmdl].keys())[0]
        names["valid_datasets"].append(
            select_metadata(cfg["input_data"].values(), alias=alias)[0][
                "filename"
            ]
        )
        if alias in trends["cmip6"].keys():
            style = plot.get_dataset_style(xtrmdl, style_file="cmip6")
        elif alias in trends["cmip5"].keys():
            style = plot.get_dataset_style(xtrmdl, style_file="cmip5")
        else:
            style = {"facecolor": (0, 0, 1, 0.2), "color": (0, 0, 1, 1.0)}

        res_ar["artrend"][xtrmdl] = np.fromiter(
            extratrends[xtrmdl].values(), dtype=float
        )
        res_ar["kde1"][xtrmdl] = stats.gaussian_kde(
            res_ar["artrend"][xtrmdl], bw_method="scott"
        )
        hbin, bins1, patches = axx.hist(
            res_ar["artrend"][xtrmdl],
            bins=res_ar["xhist"],
            density=True,
            edgecolor=style["color"],
            facecolor=style["facecolor"],
        )
        del bins1, patches
        axx.plot(
            res_ar["xval"],
            res_ar["kde1"][xtrmdl](res_ar["xval"]),
            color=style["color"],
            linewidth=3,
            label=xtrmdl,
        )
    axx_lim["maxh"] = (1.0 + 0.5 * axx_lim["factor"]) * np.max(hbin)
    names["caption"] = _plot_obs(trends, axx, axx_lim)
    names["caption"] = (
        _plot_settings(cfg, axx, period, axx_lim) + names["caption"]
    )
    fig.tight_layout()
    fig.savefig(get_plot_filename("fig2", cfg), dpi=300)
    plt.close()

    names["var"] = " ".join(list(extract_variables(cfg).keys()))
    names["caption"] = (
        "Probability density function of the decadal trend in "
        + "the "
        + extract_variables(cfg)[names["var"]]["long_name"]
        + names["caption"]
    )

    names["provenance_record"] = get_provenance_record(
        names["valid_datasets"], names["caption"], ["trend", "other"], ["reg"]
    )

    names["diagnostic_file"] = get_diagnostic_filename("fig2", cfg)

    logger.info("Saving analysis results to %s", names["diagnostic_file"])

    iris.save(
        cube_to_save_vars(_write_list_dict(cfg, trends, res_ar)),
        target=names["diagnostic_file"],
    )

    logger.info(
        "Recording provenance of %s:\n%s",
        names["diagnostic_file"],
        pformat(names["provenance_record"]),
    )
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(
            names["diagnostic_file"], names["provenance_record"]
        )


def _plot_obs(trends, axx, axx_lim):
    """Plot observational or reanalysis data as vertical line."""
    obs_str = ""
    # IPCC colors for obs from
    # https://github.com/IPCC-WG1/colormaps/blob/master/
    # categorical_colors_rgb_0-255/dark_cat.txt
    if trends["obs"]:
        obs_str = " Vertical lines show the trend for"
        for iii, obsname in enumerate(trends["obs"].keys()):
            if obsname == "RSS":
                plotobscol = (
                    221 / 255.0,
                    84 / 255.0,
                    46 / 255.0,
                )
            elif obsname == "CRU":
                plotobscol = (0.8, 0.4, 0.1, 1)
            elif obsname == "ERA5":
                plotobscol = (
                    128 / 255.0,
                    54 / 255.0,
                    168 / 255.0,
                )
            elif obsname == "MERRA2":
                plotobscol = (0.8, 0.1, 0.4, 1)
            elif obsname == "NCEP-NCAR-R1":
                plotobscol = (1, 0.6, 0, 1)
            else:
                plotobscol = _get_plotobscol(iii)

            axx.vlines(
                trends["obs"][obsname],
                0,
                axx_lim["maxh"],
                colors=plotobscol,
                linewidth=3,
                label=obsname,
            )
        obs_str = obs_str + "."
    return obs_str


def _get_ax_limits(cfg, trends):
    """Automatically find plot bounds, if 'histmin' and 'histmax' not given.

    If no values are given, the limits for the x-axis are then determined
    by the width of the histograms (plus offset).
    """
    axx_lim = {"factor": 0.3}

    if "histmin" in cfg.keys() and "histmax" in cfg.keys():
        histmin = cfg["histmin"]
        xmin = histmin
        histmax = cfg["histmax"]
        xmax = histmax
    else:
        if "histmin" in cfg.keys():
            histmin = cfg["histmin"]
            xmin = histmin
            histmax = -1000
            for label in ["cmip5", "cmip6", "obs"]:
                if trends[label]:
                    histmax = np.max(
                        np.array(
                            [
                                histmax,
                                np.max(
                                    np.fromiter(
                                        trends[label].values(), dtype=float
                                    )
                                ),
                            ]
                        )
                    )
            xmax = histmax + axx_lim["factor"] * abs(
                abs(histmax) - abs(histmin)
            )
        elif "histmax" in cfg.keys():
            histmax = cfg["histmax"]
            xmax = histmax
            histmin = 1000
            for label in ["cmip5", "cmip6", "obs"]:
                if trends[label]:
                    histmin = np.min(
                        np.array(
                            [
                                histmin,
                                np.min(
                                    np.fromiter(
                                        trends[label].values(), dtype=float
                                    )
                                ),
                            ]
                        )
                    )
            xmin = histmin - axx_lim["factor"] * abs(
                abs(histmax) - abs(histmin)
            )
        else:
            histmin = 1000
            histmax = -1000
            for label in ["cmip5", "cmip6", "obs"]:
                if trends[label]:
                    histmin = np.min(
                        np.array(
                            [
                                histmin,
                                np.min(
                                    np.fromiter(
                                        trends[label].values(), dtype=float
                                    )
                                ),
                            ]
                        )
                    )
                    histmax = np.max(
                        np.array(
                            [
                                histmax,
                                np.max(
                                    np.fromiter(
                                        trends[label].values(), dtype=float
                                    )
                                ),
                            ]
                        )
                    )
            xmin = histmin - axx_lim["factor"] * abs(
                abs(histmax) - abs(histmin)
            )
            xmax = histmax + axx_lim["factor"] * abs(
                abs(histmax) - abs(histmin)
            )

    # Saving values in dictionary axx_lim:
    axx_lim["xmin"] = xmin
    axx_lim["xmax"] = xmax
    axx_lim["histmin"] = histmin
    axx_lim["histmax"] = histmax
    axx_lim["xhist"] = np.linspace(histmin, histmax, 41)
    axx_lim["xval"] = np.linspace(xmin, xmax, 41)

    return axx_lim


def _plot_trends(cfg, trends, valid_datasets, period, axx_lim):
    """Plot probability density function of trends.

    The probability density function is estimated by using Gaussian kernel
    density estimation. The models are weighted by their respective number
    of realisations. Trends are depicted in a histogram, the probability
    density function as a curve.
    """
    res_ar = {"xval": axx_lim["xval"]}
    res_ar["xhist"] = axx_lim["xhist"]
    fig, axx = plt.subplots(figsize=(8, 6))
    maxh = -1.0
    names = {}

    # IPCC colors for CMIP5 and CMIP6 from
    # https://github.com/IPCC-WG1/colormaps/blob/master/
    # categorical_colors_rgb_0-255/cmip_cat.txt
    # CMIP5
    if trends["cmip5"]:
        maxh = _res_ar_hist("cmip5", trends, res_ar, maxh, axx)

    # CMIP6
    if trends["cmip6"]:
        maxh = _res_ar_hist("cmip6", trends, res_ar, maxh, axx)

    # Find the highest bin in the histograms for the y-axis limit
    axx_lim["maxh"] = maxh * (1.0 + 0.5 * axx_lim["factor"])
    names["caption"] = _plot_settings(cfg, axx, period, axx_lim) + _plot_obs(
        trends, axx, axx_lim
    )
    fig.tight_layout()
    fig.savefig(get_plot_filename("fig1", cfg), dpi=300)
    plt.close()

    names["var"] = " ".join(list(extract_variables(cfg).keys()))
    names["long_name"] = extract_variables(cfg)[names["var"]]["long_name"]
    names["caption"] = (
        "Probability density function of the decadal trend in "
        + "the "
        + names["long_name"]
        + names["caption"]
    )

    names["provenance_record"] = get_provenance_record(
        valid_datasets, names["caption"], ["trend", "other"], ["reg"]
    )

    names["diagnostic_file"] = get_diagnostic_filename("fig1", cfg)

    logger.info("Saving analysis results to %s", names["diagnostic_file"])

    list_dict = _make_list_dict(
        res_ar, trends, names["var"], names["long_name"]
    )

    iris.save(cube_to_save_vars(list_dict), target=names["diagnostic_file"])

    logger.info(
        "Recording provenance of %s:\n%s",
        names["diagnostic_file"],
        pformat(names["provenance_record"]),
    )
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(
            names["diagnostic_file"], names["provenance_record"]
        )


def _plot_settings(cfg, axx, period, axx_lim):
    """Define Settings for pdf figures."""
    for spine in ["top", "right"]:
        axx.spines[spine].set_visible(False)
    axx.legend(loc=1)

    if "ymax" in cfg.keys():
        ymax = cfg["ymax"]
    else:
        ymax = axx_lim["maxh"]
    axx.set_ylim([0, ymax])
    axx.set_ylabel("Probability density")
    add = "."
    if "common_period" in period.keys():
        add = " for " + period["common_period"] + "."
    elif "common_span" in period.keys():
        add = " for " + period["common_span"] + " years."
    if "common_lats" in period.keys():
        add = " between " + period["common_lats"] + add

    var = " ".join(list(extract_variables(cfg).keys()))
    print_long_name = extract_variables(cfg)[var]["long_name"]
    axx.set_title("Trends in " + print_long_name + add)
    axx.set_xlabel("Trend (%/decade)")

    axx.set_xlim(axx_lim["xmin"], axx_lim["xmax"])

    return add


def _res_ar_hist(cmipstr, trends, res_ar, maxh, axx):
    """Set and plot CMIP histograms."""
    if cmipstr == "cmip6":
        shortstr = "_c6"
        eco = (204 / 250.0, 35 / 255.0, 35 / 255.0, 1.0)
        fco = (204 / 250.0, 35 / 255.0, 35 / 255.0, 0.2)
        cco = (204 / 250.0, 35 / 255.0, 35 / 255.0, 1)
    elif cmipstr == "cmip5":
        shortstr = "_c5"
        eco = (37 / 250.0, 81 / 255.0, 204 / 255.0, 1.0)
        fco = (37 / 250.0, 81 / 255.0, 204 / 255.0, 0.2)
        cco = (37 / 250.0, 81 / 255.0, 204 / 255.0, 1)
    else:
        shortstr = ""
        eco = (107 / 250.0, 81 / 255.0, 204 / 255.0, 1.0)
        fco = (107 / 250.0, 81 / 255.0, 204 / 255.0, 0.2)
        cco = (107 / 250.0, 81 / 255.0, 204 / 255.0, 1)

    res_ar["artrend" + shortstr] = np.fromiter(
        trends[cmipstr].values(), dtype=float
    )
    res_ar["weights" + shortstr] = np.fromiter(
        trends[cmipstr + "weights"].values(), dtype=float
    )
    res_ar["kde1" + shortstr] = stats.gaussian_kde(
        res_ar["artrend" + shortstr],
        weights=res_ar["weights" + shortstr],
        bw_method="scott",
    )
    hbinx1, bins1, patches = axx.hist(
        res_ar["artrend" + shortstr],
        bins=res_ar["xhist"],
        density=True,
        weights=res_ar["weights" + shortstr],
        edgecolor=eco,
        facecolor=fco,
    )
    del bins1, patches
    axx.plot(
        res_ar["xval"],
        res_ar["kde1" + shortstr](res_ar["xval"]),
        color=cco,
        linewidth=3,
        label=cmipstr.upper,
    )
    maxh = np.max([maxh, np.max(hbinx1)])
    return maxh


def _set_extratrends_dict(cfg):
    """Set dictionary to plot pdf over model ensembles."""
    extratrends = {}
    for extramodel in cfg["add_model_dist"]:
        extratrends[extramodel] = OrderedDict()
    return extratrends


def _write_list_dict(cfg, trends, res_ar):
    """Collect data for provenance."""
    var = " ".join(list(extract_variables(cfg).keys()))
    print_long_name = extract_variables(cfg)[var]["long_name"]

    list_dict = {}
    list_dict["data"] = [res_ar["xhist"]]
    list_dict["name"] = [
        {
            "var_name": var + "_trends_bins",
            "long_name": print_long_name + " Trend bins",
            "units": "percent",
        }
    ]

    for extramodel in cfg["add_model_dist"]:
        list_dict["data"].append(res_ar["artrend"][extramodel])
        list_dict["name"].append(
            {
                "var_name": var + "_trends_" + extramodel,
                "long_name": print_long_name + " Trends " + extramodel,
                "units": "percent",
            }
        )
        list_dict["data"].append(res_ar["kde1"][extramodel](res_ar["xhist"]))
        list_dict["name"].append(
            {
                "var_name": var + "_trend_distribution_" + extramodel,
                "long_name": print_long_name
                + " Trends "
                + "distribution "
                + extramodel,
                "units": "1",
            }
        )

    if trends["obs"]:
        for obsname in trends["obs"].keys():
            list_dict["data"].append(trends["obs"][obsname])
            list_dict["name"].append(
                {
                    "var_name": var + "_trend_" + obsname,
                    "long_name": print_long_name + " Trend " + obsname,
                    "units": "percent",
                }
            )
    return list_dict


def cube_to_save_vars(list_dict):
    """Create cubes to prepare bar plot data for saving to netCDF."""
    for iii, var in enumerate(list_dict["data"]):
        if iii == 0:
            cubes = iris.cube.CubeList(
                [
                    iris.cube.Cube(
                        var,
                        var_name=list_dict["name"][iii]["var_name"],
                        long_name=list_dict["name"][iii]["long_name"],
                        units=list_dict["name"][iii]["units"],
                    )
                ]
            )
        else:
            cubes.append(
                iris.cube.Cube(
                    var,
                    var_name=list_dict["name"][iii]["var_name"],
                    long_name=list_dict["name"][iii]["long_name"],
                    units=list_dict["name"][iii]["units"],
                )
            )

    return cubes


def get_provenance_record(
    ancestor_files, caption, statistics, domains, plot_type="probability"
):
    """Get Provenance record."""
    record = {
        "caption": caption,
        "statistics": statistics,
        "domains": domains,
        "plot_type": plot_type,
        "realms": ["atmos"],
        "themes": ["atmDyn"],
        "authors": [
            "weigel_katja",
        ],
        "references": [
            "santer07jclim",
            "santer21jclim",
            "eyring21ipcc",
        ],
        "ancestors": ancestor_files,
    }
    return record


##############################################################################
# Setup diagnostic
##############################################################################


def main(cfg):
    """Run the diagnostic."""
    ##########################################################################
    # Read recipe data
    ##########################################################################

    # Dataset data containers
    input_data = cfg["input_data"].values()

    if not variables_available(cfg, ["prw"]):
        logger.warning("This diagnostic was written and tested only for prw.")

    ##########################################################################
    # Read data
    ##########################################################################

    # Create iris cube for each dataset and save annual means
    trends = {}
    trends["cmip5"] = OrderedDict()
    trends["cmip6"] = OrderedDict()
    trends["cmip5weights"] = OrderedDict()
    trends["cmip6weights"] = OrderedDict()
    trends["obs"] = {}

    if "add_model_dist" in cfg:
        extratrends = _set_extratrends_dict(cfg)

    valid_datasets, number_of_subdata, period = _get_valid_datasets(input_data)

    for dataset_path in input_data:
        project = dataset_path["project"]
        cube = _apply_filter(cfg, iris.load(dataset_path["filename"])[0])
        cat.add_month_number(cube, "time", name="month_number")
        cat.add_year(cube, "time", name="year")
        alias = dataset_path["alias"]
        dataset = dataset_path["dataset"]

        if not _check_full_data(dataset_path, cube):
            continue
        cube_anom = _calculate_anomalies(cube)

        if project == "CMIP6":
            trend = _calc_trend(cube_anom)
            trends["cmip6"][alias] = trend
            trends["cmip6weights"][alias] = 1.0 / number_of_subdata[dataset]
        elif project == "CMIP5":
            trend = _calc_trend(cube_anom)
            trends["cmip5"][alias] = trend
            trends["cmip5weights"][alias] = 1.0 / number_of_subdata[dataset]
        else:
            trend = _calc_trend(cube_anom)
            trends["obs"][dataset] = trend

        if "add_model_dist" in cfg:
            if dataset in cfg["add_model_dist"]:
                extratrends[dataset][alias] = trend

    axx_lim = _get_ax_limits(cfg, trends)
    _plot_trends(cfg, trends, valid_datasets, period, axx_lim)
    if "add_model_dist" in cfg:
        _plot_extratrends(cfg, extratrends, trends, period, axx_lim)

    ##########################################################################
    # Process data
    ##########################################################################


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
