#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 08:06:54 2019

@author: bmueller
"""

import iris.plot as iplt
import iris
import cartopy.crs as ccrs
import logging
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import textwrap
from .utilities import adjust_minmax, cmap_switch, mean_std_txt
from .utilities import get_filenames, clean_filename
import random
import string

logger = logging.getLogger(os.path.basename(__file__))


def __plot_map__(data, gs, **kwargs):
    """ 
    basic 2D plotting routine
    -------------------------
    unified 2D map plotting routine
    """
    
    ax=plt.subplot(gs, projection = ccrs.Robinson()) # does not work
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
    def randomString(stringLength=10):
        """Generate a random string of fixed length """
        letters = string.ascii_lowercase
        return ''.join(random.choice(letters) for i in range(stringLength))
    
    figsize = (10,5)
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 1)
    __plot_map__(data, gs[0,0])
    plt.colorbar()
    fig.savefig(fname = ("{}" + os.sep +
                         "simple_plot_{}.{}").format(kwargs["plotdir"],
                                         randomString(), "png"))
    return


def glob_temp_mean(data, **kwargs):
    """
    produces global_temporal_mean plots
    -----------------------------------
    produces the figure
    """
    # predefined settings
    numcols = 11
    figsize = (10,5)
    xpos_add_txt = [.05,.80]
    
    # adjust min and max values
    if "vminmax" in kwargs:
        vminmax = kwargs.pop("vminmax")
        vmin = np.nanmin(vminmax)
        vmax = np.nanmax(vminmax)
        ticks = np.linspace(vmin, vmax, num=numcols)
    else: 
        vmin = None
        vmax = None
        ticks = None
        
    # check for format
    fformat = kwargs.pop("fformat","pdf")
        
    # check for plotdir
    plotdir = kwargs.pop("plotdir",".")
        
        
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
    __plot_map__(data, gs[0,0], vmin=vmin, vmax=vmax, **kwargs)
    fig.text(xpos_add_txt[0], .05, mst, ha='left', fontsize=8)
    fig.text(xpos_add_txt[1], .05, file, ha='right', fontsize=8)
    plt.colorbar(ticks=ticks)
    plt.title("Temporal Mean of {} [{}]".format(data.long_name,
              data.units),
              pad=10, fontsize=12)
    plt.tight_layout()
    fig.savefig(fname = ("{}" + os.sep +
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
    figsize = (10,5)
    xpos_add_txt = [.05,.80]
    
    # adjust min and max values
    if "vminmax" in kwargs:
        vmax = np.abs(np.diff([np.nanmin(kwargs["vminmax"]), 
                                np.nanmax(kwargs["vminmax"])]))
        vmin = -vmax
        kwargs.pop("vminmax")
    
    else:
        vmin = data.collapsed(["longitude", "latitude"], iris.analysis.MIN)
        vmax = data.collapsed(["longitude", "latitude"], iris.analysis.MAX)
        
        vmin, vmax = adjust_minmax([vmin.data, vmax.data], symmetric = True)
    
    ticks = np.linspace(vmin, vmax, num=numcols+1)
        
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
    __plot_map__(data, gs[0,0], vmin=vmin, vmax=vmax, **kwargs)
    fig.text(xpos_add_txt[0], .05, mst, ha='left', fontsize=8)
    fig.text(xpos_add_txt[1], .05, "\n".join(files), ha='right', fontsize=8)
    plt.colorbar(ticks=ticks)
    plt.title(textwrap.fill(
            "Temporal Mean absolute difference of {} [{}]".format(
                            data.long_name,
                            data.units), 60),
                pad=10, fontsize=12)
    plt.tight_layout()
    fig.savefig(fname = ("{}" + os.sep +
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
    figsize = (10,5)
    xpos_add_txt = [.05,.80]
    extend = None
    
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
        
        vmin, vmax = adjust_minmax([vmin.data, vmax.data], symmetric = True)
    
    ticks = np.linspace(vmin, vmax, num=numcols+1)
        
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
    __plot_map__(data, gs[0,0], vmin=vmin, vmax=vmax, **kwargs)
    fig.text(xpos_add_txt[0], .05, mst, ha='left', fontsize=8)
    fig.text(xpos_add_txt[1], .05, "\n".join(files), ha='right', fontsize=8)
    plt.colorbar(ticks=ticks, extend=extend)
    plt.title(textwrap.fill(
            "Temporal Mean relative difference of {} [{}]".format(
                            data.long_name,
                            data.units), 60),
                pad=10, fontsize=12)
    plt.tight_layout()
    fig.savefig(fname = ("{}" + os.sep +
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
    figsize = (10,5)
    
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

    files = [pcol.split(".") [0] for pcol in pdf.columns[1:]]

    fig.savefig(fname = ("{}" + os.sep +
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
    figsize = (14,5)
    xpos_add_txt = [.02,.48,.52,.98]
    
    # adjust min and max values
    if "vminmax" in kwargs:
        vminmax = kwargs.pop("vminmax")
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
    __plot_map__(ref, gs[1,0:2], vmin=vmin, vmax=vmax, **kwargs)
    __plot_map__(nonref, gs[1,2:4], vmin=vmin, vmax=vmax, **kwargs)
    fig.text(xpos_add_txt[0], .15, mst_ref, ha='left', fontsize=8)
    fig.text(xpos_add_txt[1], .15, ref_file, ha='right', fontsize=8)
    fig.text(xpos_add_txt[2], .15, mst_nonref, ha='left', fontsize=8)
    fig.text(xpos_add_txt[3], .15, nonref_file, ha='right', fontsize=8)
    cax = plt.subplot(gs[3,1:3])
    plt.colorbar(cax = cax, ticks=ticks, orientation="horizontal")
    
    plt.tight_layout()
    
    fig.savefig(fname = ("{}" + os.sep +
                         "percentiles_{}_{}.{}").format(plotdir,
                                         percentile.replace(".", ""),
                                         ref_file.split(".")[0] + "_" +
                                         nonref_file.split(".")[0],
                                         fformat))
    
    return 

def __perc_format__(x, ticknum):
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
    produces the figures
    """
    #TODO: implement
    simple_plot(data["trend"], **kwargs)
    simple_plot(data["p-value"], **kwargs)
    try:
        simple_plot(data["threshold"], **kwargs)
    except:
        pass
    return 

def anomalytrend(data, **kwargs):
    """
    produces trend map, restricted to threshold or with additional p-values
    -----------------------------------------------------------------------
    produces the figures
    """
    #TODO: implement
    simple_plot(data["trend"], **kwargs)
    simple_plot(data["p-value"], **kwargs)
    try:
        simple_plot(data["threshold"], **kwargs)
    except:
        pass
    return