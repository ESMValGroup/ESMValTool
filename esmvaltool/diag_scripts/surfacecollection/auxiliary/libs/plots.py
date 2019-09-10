#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 08:06:54 2019

@author: bmueller
"""
import logging
import os
import random
import string
import textwrap

import cartopy.crs as ccrs
import cf_units
import iris
import iris.plot as iplt
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm

from .utilities import adjust_minmax, cmap_switch, mean_std_txt
from .utilities import get_filenames, clean_filename, multiplot_gridspec

logger = logging.getLogger(os.path.basename(__file__))


def __plot_map__(data, gs, **kwargs):
    """
    basic 2D plotting routine
    -------------------------
    unified 2D map plotting routine
    """

    ax = plt.subplot(gs, projection=ccrs.Robinson())  # does not work
    iplt.pcolormesh(data, **kwargs)
    ax.coastlines(linewidth=0.75, color='navy')
    ax.gridlines(crs=ccrs.PlateCarree(), linestyle=':')

    return


def simple_plot(data, **kwargs):
    """
    simple 2D plotting routine
    --------------------------
    for development use only
    """

    def random_string(stringLength=10):
        """Generate a random string of fixed length """
        letters = string.ascii_lowercase
        return ''.join(random.choice(letters) for i in range(stringLength))

    figsize = (10, 5)
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 1)
    __plot_map__(data, gs[0, 0])
    plt.colorbar()
    fig.savefig(fname=("{}" + os.sep +
                       "simple_plot_{}.{}").format(kwargs["plotdir"],
                                                   random_string(), "png"))
    return


def time_series(data, **kwargs):
    """
    produces time_series plots
    --------------------------
    produces the figure
    """

    # predefined settings
    basic_figsize = (15, 5)

    # check for format
    fformat = kwargs.pop("fformat", "pdf")

    # check for plotdir
    plotdir = kwargs.pop("plotdir", ".")

    # get size informations
    num_elements = len(data["global"].metadata.attributes["elements"])
    num_plots = len(data)
    ylim = kwargs.pop("vminmax", [None, None])

    gs, figsize = multiplot_gridspec(num_elements, num_plots, basic_figsize)

    fig = plt.figure(figsize=figsize)

    col = plt.get_cmap("jet", num_elements)

    files = []

    for jindx, sp in enumerate(gs):
        if jindx < len(data):
            ax = plt.subplot(sp)
            plt.sca(ax)
            name = list(data.keys())[jindx]

            if np.ma.is_masked(data[name].data):
                if np.all(data[name].data.mask):
                    data[name].data = np.ma.MaskedArray(
                        (data[name].data.data * 0 + 1) * -1e16,
                        mask=np.logical_not(data[name].data.mask))

            for indx, cube in enumerate(data[name].slices('unified_time')):
                cube_label = data[name].metadata.attributes["elements"][
                    cube.coord("element").points[0]].split(".")[0]
                if jindx == 0:
                    files.append(cube_label)
                iplt.plot(cube, label=cube_label, color=col(indx))

            ax.autoscale(tight=True)

            ax.set(
                ylim=ylim,
                xlabel='Time',
                ylabel=textwrap.fill("{} [{}]".format(cube.long_name,
                                                      cube.units), 50),
                title=name
            )

            plt.grid(True)

            plt.tight_layout()

    all_axes = fig.get_axes()

    for indx, ax in enumerate(all_axes):
        if indx not in [num_plots - 2, num_plots - 1]:
            ax.set_xticklabels([])
            ax.set(xlabel='')
        if indx % 2:
            ax.set_yticklabels([])
            ax.set(ylabel='')

    # Add the legend
    plt.legend(bbox_to_anchor=(1, 0), loc="lower right",
               bbox_transform=fig.transFigure, ncol=2)

    fig.savefig(fname=("{}" + os.sep +
                       "time_series_{}.{}").format(plotdir,
                                                   "_".join(files), fformat))

    return


def __time_label_and_loc__(time_coord):
    """
    produces labels and locators for time axis
    ------------------------------------------
    returns locations and labels
    """

    locs = [tp for tp in time_coord.points if
            np.isclose(tp % 365.2425, 0,
                       atol=np.mean(
                           np.diff(time_coord.points)))]

    while len(locs) > 7:
        locs = [locs[lind] for (lind, _) in
                enumerate(locs) if not lind % 2]

    locs = locs + [time_coord.points[-1]]

    labels = cf_units.num2date(locs,
                               time_coord.units.name,
                               time_coord.units.calendar)
    for (idx, _) in enumerate(labels):
        labels[idx] = labels[idx].strftime('%Y-%m-%d')

    locs = [l - 365.2425 * 50 for l in locs]

    return locs, labels


def glob_temp_mean(data, **kwargs):
    """
    produces global_temporal_mean plots
    -----------------------------------
    produces the figure
    """
    # predefined settings
    numcols = 11
    figsize = (10, 5)
    xpos_add_txt = [.05, .80]

    # adjust min and max values
    if "vminmax" in kwargs:
        vminmax = kwargs.pop("vminmax")
        if all(np.isnan(vminmax)):
            vmin = None
            vmax = None
            ticks = None
        else:
            vmin = np.nanmin(vminmax)
            vmax = np.nanmax(vminmax)
            ticks = np.linspace(vmin, vmax, num=numcols)
    else:
        vmin = None
        vmax = None
        ticks = None

    # check for format
    fformat = kwargs.pop("fformat", "pdf")

    # check for plotdir
    plotdir = kwargs.pop("plotdir", ".")

    file = data.metadata.attributes["source_file"].split(os.sep)[-1]

    # adjust colorbar
    if "cmap" in kwargs:
        kwargs["cmap"] = cmap_switch(kwargs["cmap"], numcols)
    else:
        kwargs["cmap"] = None

    mst = mean_std_txt(data)

    # plot data
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 1)
    __plot_map__(data, gs[0, 0], vmin=vmin, vmax=vmax, **kwargs)
    fig.text(xpos_add_txt[0], .05, mst, ha='left', fontsize=8)
    fig.text(xpos_add_txt[1], .05, file, ha='right', fontsize=8)
    plt.colorbar(ticks=ticks)
    plt.title("Temporal Mean of {} [{}]".format(data.long_name,
                                                data.units),
              pad=10, fontsize=12)
    plt.tight_layout()
    fig.savefig(fname=("{}" + os.sep +
                       "glob_temp_mean_{}.{}").format(plotdir,
                                                      file.split(".")[0], fformat))

    return


def glob_temp_mean_absdiff(data, **kwargs):
    """
    produces global temporal mean absolute difference plots
    -------------------------------------------------------
    produces the figure
    """
    # predefined settings
    numcols = 10
    figsize = (10, 5)
    xpos_add_txt = [.05, .80]

    # adjust min and max values
    if "vminmax" in kwargs:
        vmax = np.abs(np.diff([np.nanmin(kwargs["vminmax"]),
                               np.nanmax(kwargs["vminmax"])]))
        vmin = -vmax
        kwargs.pop("vminmax")

    else:
        vmin = data.collapsed(["longitude", "latitude"], iris.analysis.MIN)
        vmax = data.collapsed(["longitude", "latitude"], iris.analysis.MAX)

        vmin, vmax = adjust_minmax([vmin.data, vmax.data], symmetric=True)

    ticks = np.linspace(vmin, vmax, num=numcols + 1)

    # check for format
    fformat = kwargs.pop("fformat", "pdf")

    # check for plotdir
    plotdir = kwargs.pop("plotdir", ".")

    files = get_filenames(data.metadata.attributes["source_file"])

    # adjust colorbar
    kwargs["cmap"] = plt.cm.get_cmap("coolwarm_r", numcols)
    logging.warning("colormap for differences predefined; cfg overwritten!")

    mst = mean_std_txt(data)

    # plot data
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 1)
    __plot_map__(data, gs[0, 0], vmin=vmin, vmax=vmax, **kwargs)
    fig.text(xpos_add_txt[0], .05, mst, ha='left', fontsize=8)
    fig.text(xpos_add_txt[1], .05, "\n".join(files), ha='right', fontsize=8)
    plt.colorbar(ticks=ticks)
    plt.title(textwrap.fill(
        "Temporal Mean absolute difference of {} [{}]".format(
            data.long_name,
            data.units), 60),
        pad=10, fontsize=12)
    plt.tight_layout()
    fig.savefig(fname=("{}" + os.sep +
                       "glob_temp_mean_absdiff_{}.{}").format(plotdir,
                                                              clean_filename("_".join(
                                                                  f.split(".")[0] for
                                                                  f in files)), fformat))

    return


def glob_temp_mean_reldiff(data, **kwargs):
    """
    produces global temporal mean relative difference plots
    -------------------------------------------------------
    produces the figure
    """
    # predefined settings
    numcols = 10
    figsize = (10, 5)
    xpos_add_txt = [.05, .80]
    extend = "neither"

    # adjust min and max values
    if "vminmax" in kwargs:
        vmin = -1
        vmax = 1
        logging.warning("ranges for relative differences predefined;" +
                        " cfg overwritten!")
        extend = "both"
        kwargs.pop("vminmax")

    else:
        vmin = data.collapsed(["longitude", "latitude"], iris.analysis.MIN)
        vmax = data.collapsed(["longitude", "latitude"], iris.analysis.MAX)

        vmin, vmax = adjust_minmax([vmin.data, vmax.data], symmetric=True)

    ticks = np.linspace(vmin, vmax, num=numcols + 1)

    # check for format
    fformat = kwargs.pop("fformat", "pdf")

    # check for plotdir
    plotdir = kwargs.pop("plotdir", ".")

    files = get_filenames(data.metadata.attributes["source_file"])

    logger.info(files)

    # adjust colorbar
    kwargs["cmap"] = plt.cm.get_cmap("coolwarm_r", numcols)
    logging.warning("colormap for differences predefined; cfg overwritten!")

    mst = mean_std_txt(data)

    # plot data
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 1)
    __plot_map__(data, gs[0, 0], vmin=vmin, vmax=vmax, **kwargs)
    fig.text(xpos_add_txt[0], .05, mst, ha='left', fontsize=8)
    fig.text(xpos_add_txt[1], .05, "\n".join(files), ha='right', fontsize=8)
    plt.colorbar(ticks=ticks, extend=extend)
    plt.title(textwrap.fill(
        "Temporal Mean relative difference of {} [{}]".format(
            data.long_name,
            data.units), 60),
        pad=10, fontsize=12)
    plt.tight_layout()
    fig.savefig(fname=("{}" + os.sep +
                       "glob_temp_mean_reldiff_{}.{}").format(plotdir,
                                                              clean_filename("_".join(
                                                                  f.split(".")[0] for
                                                                  f in files)), fformat))

    return


def percentiles(data, **kwargs):
    """
    produces percentiles map and line plots
    ---------------------------------------
    deligates figures
    """

    if isinstance(data, pd.DataFrame):
        __perc_lin_plot__(data, **kwargs)
    elif isinstance(data, dict):
        kwargs["percentile"] = list(data.keys())[0]
        r = (list(data.values())[0]).get("ref")
        for ind, cube in enumerate(list(data.values())[0].get("nonref")):
            kwargs["corr"] = (list(data.values())[0]).get("corr")[ind]
            nr = cube
            __perc_map_plot__(r, nr, **kwargs)
    else:
        logger.error("This is no plottable feature!")

    return


def __perc_lin_plot__(pdf, **kwargs):
    """
    produces percentiles line plot with all percentile correlations
    ---------------------------------------------------------------
    produces the figure
    """

    # predefined settings
    figsize = (10, 5)

    # check for format
    fformat = kwargs.pop("fformat", "pdf")

    # check for plotdir
    plotdir = kwargs.pop("plotdir", ".")

    pdf = pdf.reset_index()

    fig, ax = plt.subplots(figsize=figsize)

    ax.plot([], [], ' ', label="ref: {}".format(pdf.columns[1]))

    col = plt.get_cmap("jet", len(pdf.columns[2:]))

    for i, c in enumerate(pdf.columns[2:]):
        ax.plot("index", c, "D--", data=pdf, color=col(i), label=c)
    ax.set_ylim(0., 1.)
    ax.set_ylabel("correlation")
    ax.set_xlim(0., 1.)
    ax.set_xlabel("percentile")
    ax.set_xticks(pdf["index"])
    ax.xaxis.set_major_formatter(plt.FuncFormatter(__perc_format__))
    ax.grid(linestyle=":")
    ax.legend()

    plt.tight_layout()

    files = [pcol.split(".")[0] for pcol in pdf.columns[1:]]

    fig.savefig(fname=("{}" + os.sep +
                       "percentile_correlation_{}.{}").format(plotdir,
                                                              "_".join(files), fformat))
    return


def __perc_map_plot__(ref, nonref, **kwargs):
    """
    produces percentiles map plot with correlation info
    ---------------------------------------------------
    produces the figure
    """

    data = nonref

    # predefined settings
    numcols = 11
    figsize = (14, 5)
    xpos_add_txt = [.02, .48, .52, .98]

    # adjust min and max values
    if "vminmax" in kwargs:
        vminmax = kwargs.pop("vminmax")
        if all(np.isnan(vminmax)):
            vmin = None
            vmax = None
            ticks = None
        else:
            vmin = np.nanmin(vminmax)
            vmax = np.nanmax(vminmax)
            ticks = np.linspace(vmin, vmax, num=numcols)
    else:
        vmin = None
        vmax = None
        ticks = None

    # check for percentile
    percentile = kwargs.pop("percentile", None)

    # check for correlation values
    if "corr" in kwargs:
        rval = kwargs["corr"]["r"]
        # pval = kwargs["corr"]["p-value"]
        kwargs.pop("corr")
    else:
        rval = None

    # check for format
    fformat = kwargs.pop("fformat", "pdf")

    # check for plotdir
    plotdir = kwargs.pop("plotdir", ".")

    ref_file = ref.metadata.attributes["source_file"].split(os.sep)[-1]
    nonref_file = nonref.metadata.attributes["source_file"].split(os.sep)[-1]

    # adjust colorbar
    if "cmap" in kwargs:
        kwargs["cmap"] = cmap_switch(kwargs["cmap"], numcols)
    else:
        kwargs["cmap"] = None

    mst_ref = mean_std_txt(ref)
    mst_nonref = mean_std_txt(nonref)

    # plot data
    fig = plt.figure(figsize=figsize)
    fig.suptitle(r"{0:.0%}-Percentiles of {1} [{2}] (R = {3:3.2f})".format(
        float(percentile),
        data.long_name,
        data.units,
        rval,
    ),
        #                  pad=10,
        fontsize=12)
    gs = gridspec.GridSpec(4, 4,
                           width_ratios=[1, 5, 5, 1],
                           height_ratios=[1, 80, 1, 4])
    __plot_map__(ref, gs[1, 0:2], vmin=vmin, vmax=vmax, **kwargs)
    __plot_map__(nonref, gs[1, 2:4], vmin=vmin, vmax=vmax, **kwargs)
    fig.text(xpos_add_txt[0], .15, mst_ref, ha='left', fontsize=8)
    fig.text(xpos_add_txt[1], .15, ref_file, ha='right', fontsize=8)
    fig.text(xpos_add_txt[2], .15, mst_nonref, ha='left', fontsize=8)
    fig.text(xpos_add_txt[3], .15, nonref_file, ha='right', fontsize=8)
    cax = plt.subplot(gs[3, 1:3])
    plt.colorbar(cax=cax, ticks=ticks, orientation="horizontal")

    plt.tight_layout()

    fig.savefig(fname=("{}" + os.sep +
                       "percentiles_{}_{}.{}").format(plotdir,
                                                      percentile.replace(".", ""),
                                                      ref_file.split(".")[0] + "_" +
                                                      nonref_file.split(".")[0],
                                                      fformat))

    return


def __perc_format__(x):
    """
    reformats tick labels to percent
    --------------------------------
    returns % string
    """

    return r"{0:.0%}".format(x)


def trend(data, **kwargs):
    """
    produces trend map, restricted to threshold or with additional p-values
    -----------------------------------------------------------------------
    deligates figures
    """

    if "threshold" in data.keys():
        __threshold_trend_plot__(data["threshold"], **kwargs)
        data.pop("threshold")

    __standard_trend_plot__(data, **kwargs)

    return


def anomalytrend(data, **kwargs):
    """
    produces trend map, restricted to threshold or with additional p-values
    -----------------------------------------------------------------------
    deligates figures
    """

    trend(data, **kwargs)
    return


def climatologytrend(data, **kwargs):
    """
    produces trend map, restricted to threshold or with additional p-values
    -----------------------------------------------------------------------
    deligates figures
    """

    trend(data, **kwargs)
    return


def __standard_trend_plot__(data, **kwargs):
    """
    produces trend map with additional p-values
    -------------------------------------------
    produces the figures
    """
    # predefined settings
    numcols = 11
    figsize = (14, 5)
    xpos_add_txt = [.02, .48, .52, .98]
    extend = 'neither'

    # adjust min and max values
    if "vminmax" in kwargs:
        vminmax = kwargs.pop("vminmax")
        if all(np.isnan(vminmax)):
            vmin = None
            vmax = None
            ticks = None
        else:
            vmin = np.nanmin(vminmax)
            vmax = np.nanmax(vminmax)
            ticks = np.linspace(vmin, vmax, num=numcols)
            extend = "both"
    else:
        vmin = None
        vmax = None
        ticks = None

    # check for format
    fformat = kwargs.pop("fformat", "pdf")

    # check for plotdir
    plotdir = kwargs.pop("plotdir", ".")

    file = data["trend"].metadata.attributes["source_file"].split(os.sep)[-1]

    # adjust colorbar
    if "cmap" in kwargs:
        kwargs["cmap"] = cmap_switch(kwargs["cmap"], numcols)
    else:
        kwargs["cmap"] = None

    mst = mean_std_txt(data["trend"])

    # plot data
    fig = plt.figure(figsize=figsize)
    fig.suptitle(textwrap.fill("{} [{}] with accompanying p-values".format(
        data["trend"].long_name,
        data["trend"].units),
        160), fontsize=12)
    gs = gridspec.GridSpec(4, 2,
                           width_ratios=[1, 1],
                           height_ratios=[1, 80, 1, 4])
    __plot_map__(data["trend"], gs[1, 0], vmin=vmin, vmax=vmax, **kwargs)
    cax1 = plt.subplot(gs[3, 0])
    plt.colorbar(cax=cax1, ticks=ticks,
                 orientation="horizontal", extend=extend)
    kwargs["cmap"] = cm.get_cmap('Greens_r', 10)
    __plot_map__(data["p-value"], gs[1, 1], vmin=0, vmax=1, **kwargs)
    cax2 = plt.subplot(gs[3, 1])
    plt.colorbar(cax=cax2, extend="max", orientation="horizontal")
    fig.text(xpos_add_txt[0], .15, mst, ha='left', fontsize=8)
    fig.text(xpos_add_txt[1], .15, "trend", ha='right', fontsize=10)
    fig.text(xpos_add_txt[2], .15, "p-value", ha='left', fontsize=10)
    fig.text(xpos_add_txt[3], .15, file, ha='right', fontsize=8)

    plt.tight_layout()

    fig.savefig(fname=("{}" + os.sep +
                       "glob_trend_{}.{}").format(plotdir,
                                                  file.split(".")[0], fformat))

    return


def __threshold_trend_plot__(data, **kwargs):
    """
    produces trend map restricted to threshold
    ------------------------------------------
    produces the figures
    """
    # predefined settings
    numcols = 11
    figsize = (10, 5)
    xpos_add_txt = [.05, .80]
    extend = 'neither'

    # adjust min and max values
    if "vminmax" in kwargs:
        vminmax = kwargs.pop("vminmax")
        if all(np.isnan(vminmax)):
            vmin = None
            vmax = None
            ticks = None
        else:
            vmin = np.nanmin(vminmax)
            vmax = np.nanmax(vminmax)
            ticks = np.linspace(vmin, vmax, num=numcols)
            extend = "both"
    else:
        vmin = None
        vmax = None
        ticks = None

    # check for format
    fformat = kwargs.pop("fformat", "pdf")

    # check for plotdir
    plotdir = kwargs.pop("plotdir", ".")

    file = data.metadata.attributes["source_file"].split(os.sep)[-1]

    # adjust colorbar
    if "cmap" in kwargs:
        kwargs["cmap"] = cmap_switch(kwargs["cmap"], numcols)
    else:
        kwargs["cmap"] = None

    mst = mean_std_txt(data)

    # plot data
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 1)
    __plot_map__(data, gs[0, 0], vmin=vmin, vmax=vmax, **kwargs)
    fig.text(xpos_add_txt[0], .05, mst, ha='left', fontsize=8)
    fig.text(xpos_add_txt[1], .05, file, ha='right', fontsize=8)
    plt.colorbar(ticks=ticks, extend=extend)
    plt.title(textwrap.fill("{} [{}] (p-values <= {:.2f})".format(
        data.long_name,
        data.units,
        data.attributes["pthreshold"]), 80),
        pad=10, fontsize=12)
    plt.tight_layout()
    fig.savefig(fname=("{}" + os.sep +
                       "glob_trend_thresh_{}.{}").format(plotdir,
                                                         file.split(".")[0], fformat))
    return


def correlation(data, **kwargs):
    """
    produces correlation map, restricted to threshold or with add. p-values
    -----------------------------------------------------------------------
    deligates figures
    """

    if "threshold" in data.keys():
        __threshold_correlation_plot__(data["threshold"], **kwargs)
        data.pop("threshold")

    __standard_correlation_plot__(data, **kwargs)

    return


def anomalycorrelation(data, **kwargs):
    """
    produces correlation map, restricted to threshold or with add. p-values
    -----------------------------------------------------------------------
    deligates figures
    """

    correlation(data, **kwargs)
    return


def climatologycorrelation(data, **kwargs):
    """
    produces correlation map, restricted to threshold or with add. p-values
    -----------------------------------------------------------------------
    deligates figures
    """

    correlation(data, **kwargs)
    return


def __standard_correlation_plot__(data, **kwargs):
    """
    produces correlation map with additional p-values
    -------------------------------------------
    produces the figures
    """
    # predefined settings
    numcols = 11
    figsize = (14, 5)
    xpos_add_txt = [.02, .48, .52, .98]

    # adjust min and max values
    if "vminmax" in kwargs:
        vminmax = kwargs.pop("vminmax")
        if all(np.isnan(vminmax)):
            vmin = None
            vmax = None
            ticks = None
        else:
            vmin = np.nanmin(vminmax)
            vmax = np.nanmax(vminmax)
            ticks = np.linspace(vmin, vmax, num=numcols)
    else:
        vmin = None
        vmax = None
        ticks = None

    # check for format
    fformat = kwargs.pop("fformat", "pdf")

    # check for plotdir
    plotdir = kwargs.pop("plotdir", ".")

    files = get_filenames(
        data["correlation"].metadata.attributes["source_file"])

    # adjust colorbar
    if "cmap" in kwargs:
        kwargs["cmap"] = cmap_switch(kwargs["cmap"], numcols)
    else:
        kwargs["cmap"] = None

    mst = mean_std_txt(data["correlation"])

    # plot data
    fig = plt.figure(figsize=figsize)
    fig.suptitle(textwrap.fill("{} [{}] with accompanying p-values".format(
        data["correlation"].long_name,
        data["correlation"].units),
        160), fontsize=12)
    gs = gridspec.GridSpec(4, 2,
                           width_ratios=[1, 1],
                           height_ratios=[1, 80, 1, 4])
    __plot_map__(data["correlation"], gs[1, 0], vmin=vmin, vmax=vmax, **kwargs)
    cax1 = plt.subplot(gs[3, 0])
    plt.colorbar(cax=cax1, ticks=ticks,
                 orientation="horizontal")
    kwargs["cmap"] = cm.get_cmap('Greens_r', 10)
    __plot_map__(data["p-value"], gs[1, 1], vmin=0, vmax=1, **kwargs)
    cax2 = plt.subplot(gs[3, 1])
    plt.colorbar(cax=cax2, extend="max", orientation="horizontal")
    fig.text(xpos_add_txt[0], .15, mst, ha='left', fontsize=8)
    fig.text(xpos_add_txt[1], .15, "correlation", ha='right', fontsize=10)
    fig.text(xpos_add_txt[2], .15, "p-value", ha='left', fontsize=10)
    fig.text(xpos_add_txt[3], .15, "\n".join(files), ha='right', fontsize=8)

    plt.tight_layout()

    fig.savefig(fname=("{}" + os.sep +
                       "glob_correlation_{}.{}").format(plotdir,
                                                        clean_filename("_".join(
                                                            f.split(".")[0] for
                                                            f in files)), fformat))

    return


def __threshold_correlation_plot__(data, **kwargs):
    """
    produces correlation map restricted to threshold
    ------------------------------------------
    produces the figures
    """
    # predefined settings
    numcols = 11
    figsize = (10, 5)
    xpos_add_txt = [.05, .80]

    # adjust min and max values
    if "vminmax" in kwargs:
        vminmax = kwargs.pop("vminmax")
        if all(np.isnan(vminmax)):
            vmin = None
            vmax = None
            ticks = None
        else:
            vmin = np.nanmin(vminmax)
            vmax = np.nanmax(vminmax)
            ticks = np.linspace(vmin, vmax, num=numcols)
    else:
        vmin = None
        vmax = None
        ticks = None

    # check for format
    fformat = kwargs.pop("fformat", "pdf")

    # check for plotdir
    plotdir = kwargs.pop("plotdir", ".")

    files = get_filenames(data.metadata.attributes["source_file"])

    # adjust colorbar
    if "cmap" in kwargs:
        kwargs["cmap"] = cmap_switch(kwargs["cmap"], numcols)
    else:
        kwargs["cmap"] = None

    mst = mean_std_txt(data)

    # plot data
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 1)
    __plot_map__(data, gs[0, 0], vmin=vmin, vmax=vmax, **kwargs)
    fig.text(xpos_add_txt[0], .05, mst, ha='left', fontsize=8)
    fig.text(xpos_add_txt[1], .05, "\n".join(files), ha='right', fontsize=8)
    plt.colorbar(ticks=ticks)
    plt.title(textwrap.fill("{} [{}] (p-values <= {:.2f})".format(
        data.long_name, data.units,
        data.attributes["pthreshold"]), 80),
        pad=10, fontsize=12)
    plt.tight_layout()
    fig.savefig(fname=("{}" + os.sep +
                       "glob_correlation_thresh_{}.{}").format(plotdir,
                                                               clean_filename("_".join(
                                                                   f.split(".")[0] for
                                                                   f in files)), fformat))
    return
