#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to create table of scalar data.

Description
-----------
Calculate and plot ZEC (Zero Emission Commitment) Temperature.
Requires input data from 1pctCO2 and a dedicated simulation from ZECMIP.

Configuration options in recipe
-------------------------------
option : bool, optional (default: True)
    What it does!


"""

import logging
import os
from copy import deepcopy
from pprint import pformat

import iris
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    io,
    run_diagnostic,
    select_metadata,
    sorted_metadata,
)

logger = logging.getLogger(os.path.basename(__file__))

def calculate_zec(cfg):
    """Calculate ZEC for each model."""
    zec = {}
    input_data = cfg['input_data'].values()
    input_data = sorted_metadata(input_data, ['short_name', 'exp', 'dataset'])
    # I will need to foolproof this more. ZECMIP runs can have different exps
    # base temp should always come from 1pctCO2
    # probably include dimensions checks here since base temp should be 1D
    #onepct_data = select_metadata(input_data, variable_group='tas_base', exp='1pctCO2')
    #logger.info(onepct_data)
    zecmip_data = select_metadata(input_data, variable_group='tas')
    for data in zecmip_data:
        name = data['dataset']
        logger.info('Processing %s', name)
        anom = select_metadata(input_data, short_name='tas',
                               exp = '1pctCO2', dataset=name)
        if not anom:
            raise ValueError("No anomaly base available for dataset %s",
                             name)
        
        tas = iris.load_cube(data['filename'])
        tas_anom = iris.load_cube(anom[0]['filename'])
        logger.info(tas_anom.data)
        zec_model = deepcopy(tas)
        # Fix time to start at 0
        fix_time(zec_model)
        zec_model = zec_model - tas_anom.data
        zec[name] = zec_model
    return zec

def fix_time(cube):
    """Fix time to start at year 0."""
    time_coord = iris.coords.DimCoord(
    np.arange(cube.coord('time').shape[0]),
    var_name='time',
    standard_name='time',
    long_name='time',
    units='years',
    )
    cube.remove_coord('time')
    cube.add_dim_coord(time_coord, 0)
    #return cube

def plot_zec_timeseries_individual(zec, cfg):
    """Plot individual ZEC timeseries per model."""
    for model in zec:
        fig, axes = plt.subplots(figsize=(9,7))
        cube = zec[model]
        iris.plot.plot(cube, axes=axes)
        axes.axhline(color="lightgrey", linestyle="--")     
        axes.set_title('%s ZEC'%model)
        axes.set_xlabel('Time [yr]')
        axes.set_ylabel('ZEC [K]')
        axes.set_xlim([0, np.max(cube.coord('time').points[-1])])

        # Save plot
        plot_path = get_plot_filename('zec_timeseries_%s'%model, cfg)
        fig.savefig(plot_path)
        plt.close()            
    

def plot_zec_timeseries_mmm(zec, cfg):
    """Plot all ZEC timeseries with multi-model mean."""
    fig, axes = plt.subplots(figsize=(10,6))
    for model in zec:
        iris.plot.plot(zec[model], axes=axes, label=model)
        axes.axhline(color="lightgrey", linestyle="--")     
    axes.set_title('ZEC')
    axes.set_xlabel('Time [yr]')
    axes.set_ylabel('ZEC [K]')
    axes.set_xlim([0, 100])
    axes.legend(bbox_to_anchor=(0.10, -0.09), loc="upper left",
            ncol=3, handlelength=3.5, fontsize=12)

    # Save plot
    plot_path = get_plot_filename('zec_timeseries_all_models', cfg)
    plt.tight_layout()
    fig.savefig(plot_path)
    plt.close()      

def calc_zec_x(zec, x):
    zec_x = {}
    for model in zec:
        # ZEC_X is the 20-year anomaly centered at year x
        zec_x[model] = np.mean(zec[model].data[x-10:x+9])
    return zec_x

def plot_zec_x_bar(zec_x, x_i, cfg):
    """Plot a barplot of all ZEC_x in ascending order"""
    logger.info('WIP')
    # Sort by value
    data = dict(sorted(zec_x.items(), key=lambda x: x[1]))
    # Make plot
    fig, axes = plt.subplots(figsize=(10,6))
    logger.info(list(data.keys()))
    logger.info(list(data.values()))
    axes.bar(list(data.keys()), list(data.values()))


    axes.axhline(color="lightgrey", linestyle="--")

    axes.set_ylabel(r'ZEC$_{%s}$'%str(x_i))
    # Save plot
    plot_path = get_plot_filename('zec_%s_barplot'%str(x_i), cfg)
    plt.tight_layout()
    fig.savefig(plot_path)
    plt.close()   

def main(cfg):
    """Execute diagnostic."""
    # Calculate ZEC
    zec = calculate_zec(cfg)

    # Plot ZEC Timeseries
    plot_zec_timeseries_individual(zec, cfg)
    plot_zec_timeseries_mmm(zec, cfg)
    
    # Calculate ZEC_X
    x = [50]
    for x_i in x:
        zec_x = calc_zec_x(zec, x_i)
        plot_zec_x_bar(zec_x, x_i, cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)