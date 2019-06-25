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
import numpy as np

logger = logging.getLogger(os.path.basename(__file__))

def glob_temp_mean(data, **kwargs):
    """
    produces global_temporal_mean plots
    -----------------------------------
    returns a list of mean cubes
    """
    # predefined settings
    numcols = 11
    figsize = (10,5)
    xpos_add_txt = [.02,.80]
    
    # adjust min and max values
    if "vminmax" in kwargs:
        vminmax = kwargs["vminmax"]
        kwargs.pop("vminmax")
        vmin = np.nanmin(vminmax)
        vmax = np.nanmax(vminmax)
        ticks = np.linspace(vmin, vmax, num = numcols)
    else: 
        vmin = None
        vmax = None
        ticks = None
        
    # check for format
    if "fformat" in kwargs:
        fformat = kwargs["fformat"]
        kwargs.pop("fformat")
    else:
        fformat = "pdf"
        
    # check for plotdir
    if "plotdir" in kwargs:
        plotdir = kwargs["plotdir"]
        kwargs.pop("plotdir")
    else:
        plotdir = "."
        
        
    file = data.metadata.attributes["source_file"].split(os.sep)[-1]
        
    # adjust colorbar
    if "cmap" in kwargs:
        kwargs["cmap"] = plt.cm.get_cmap(kwargs["cmap"], numcols-1)
    else: 
        kwargs["cmap"] = None
        
    weights = iris.analysis.cartography.area_weights(data)
    mean = data.collapsed(["longitude", "latitude"],
                          iris.analysis.MEAN, weights = weights)
    std = data.collapsed(["longitude", "latitude"],
                          iris.analysis.STD_DEV)
    
    mean_std_txt = ("spatial statistics: \n" +
                    " \u03BC: {:.3g} \n" + 
                    " \u03C3: {:.3g}").format(mean.data, std.data)
    
    # plot data
    fig = plt.figure(figsize=figsize)
    __plot_map__(data, vmin=vmin, vmax=vmax, **kwargs)
    fig.text(xpos_add_txt[0], .05, mean_std_txt, ha='left')
    fig.text(xpos_add_txt[1], .05, file, ha='right')
    plt.colorbar(ticks=ticks)
    plt.title("Temporal Mean of {} [{}]".format(data.long_name,
              data.units),
              pad=10, fontsize=12)
    plt.tight_layout()
    fig.savefig(fname = ("{}" + os.sep +
                         "glob_temp_mean_{}.{}").format(plotdir,
                                         file.split(".")[0], fformat))
    
    return 

def __plot_map__(data, **kwargs):
    """ 
    basic 2D plotting routine
    -------------------------
    unified 2D map plotting routine
    """
    
    ax = plt.axes(projection=ccrs.Robinson())
    ax.coastlines(linewidth=0.75, color='navy')
    ax.gridlines(crs=ccrs.PlateCarree(), linestyle='-')
    iplt.pcolormesh(data, **kwargs)
    
    return ax