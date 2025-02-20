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

import iris
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    io,
    run_diagnostic,
    select_metadata,
)

logger = logging.getLogger(os.path.basename(__file__))


def _get_default_cfg(cfg):
    """Get default options for configuration dictionary."""
    cfg = deepcopy(cfg)

    cfg.setdefault('zec_x', [50])
    cfg.setdefault(
        'experiments', {
            'reference': ['esm-flat10', '1pctCO2'],
            'simulation': ['esm-flat10-zec', 'esm-1pct-brch-1000PgC']
        })
    return cfg


def group_data(cfg):
    """Groups input data in preparation of ZEC calculation."""
    input_data = cfg['input_data'].values()
    group_data = group_metadata(input_data, 'exp')
    for exp in group_data:
        logger.info(exp)
        if exp in cfg['experiments']['simulation']:
            zecmip_data = select_metadata(input_data,
                                          short_name='tas',
                                          exp=exp)
        elif exp in cfg['experiments']['reference']:
            anom = select_metadata(input_data, short_name='tas', exp=exp)
        else:
            logger.info('in raise error')
            raise ValueError(
                f"{exp} is not a valid experiment for calculating ZEC, "
                f"please check the configuration value of 'experiments'. "
                f"Current accepted experiments are {cfg['experiments']}")
    return zecmip_data, anom


def calculate_zec(cfg):
    """Calculate ZEC for each model."""
    zec = {}
    zecmip_data, anom = group_data(cfg)
    for data in zecmip_data:
        # Account for ensembles by using alias, remove exp name
        name = data['alias'].replace('_' + data['exp'], '')
        logger.info('Processing %s' % name)
        tas = iris.load_cube(data['filename'])
        # Match the correct anomaly data
        match_anom = select_metadata(anom,
                                     dataset=data['dataset'],
                                     ensemble=data['ensemble'])
        tas_anom = iris.load_cube(match_anom[0]['filename'])
        zec_model = deepcopy(tas)
        # Fix time to start at 0
        fix_time(zec_model)
        zec[name] = zec_model - tas_anom.data
    save_zec(zec, cfg)
    return zec


def save_zec(zec, cfg):
    """Write ZEC data to file."""
    netcdf_path = get_diagnostic_filename("zec", cfg)
    var_attrs = {
        'short_name': 'zec',
        'long_name': 'Zero Emissions Commitment (ZEC)',
        'units': 'K'
    }
    io.save_1d_data(zec, netcdf_path, 'time', var_attrs)
    write_provenance(cfg, netcdf_path)


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


def plot_zec_timeseries(zec, cfg):
    """Plot all ZEC timeseries."""
    fig, axes = plt.subplots(figsize=(10, 6))
    for model in zec:
        iris.plot.plot(zec[model], axes=axes, label=model)
        axes.axhline(color="lightgrey", linestyle="--")
    axes.set_title('ZEC')
    axes.set_xlabel('Time [yr]')
    axes.set_ylabel('ZEC [K]')
    axes.set_xlim([0, 100])
    axes.legend(bbox_to_anchor=(0.10, -0.09),
                loc="upper left",
                ncol=3,
                handlelength=3.5,
                fontsize=12)

    # Save plot
    plot_path = get_plot_filename('zec_timeseries_all_models', cfg)
    plt.tight_layout()
    fig.savefig(plot_path)
    plt.close()
    logger.info("Wrote %s", plot_path)
    prov_dict = {'caption': 'ZEC timeseries', 'plot_type': 'times'}
    write_provenance(cfg, plot_path, prov_dict)


def calc_zec_x(zec, x):
    zec_x = {}
    for model in zec:
        # ZEC_X is the 20-year anomaly centered at year x
        zec_x[model] = np.mean(zec[model].data[x - 10:x + 9])
    return zec_x


def write_zec_x(zec_x, x, cfg):
    netcdf_path = get_diagnostic_filename("zec_%s" % str(x), cfg)
    var_attrs = {
        'short_name': 'zec',
        'long_name': 'Zero Emissions Commitment (ZEC)',
        'units': 'K'
    }
    io.save_scalar_data(zec_x, netcdf_path, var_attrs)
    write_provenance(cfg, netcdf_path)


def plot_zec_x_bar(zec_x, x_i, cfg):
    """Plot a barplot of all ZEC_x in ascending order."""
    # Sort by value
    data = dict(sorted(zec_x.items(), key=lambda x: x[1]))
    labels = [x.replace('_', ' \n ') for x in data]
    # Make plot
    fig, axes = plt.subplots(figsize=(10, 6))
    axes.bar(labels, list(data.values()))
    axes.axhline(color="lightgrey", linestyle="--")
    axes.set_ylabel(r'ZEC$_{%s}$' % str(x_i))
    # Save plot
    plot_path = get_plot_filename('zec_%s_barplot' % str(x_i), cfg)
    plt.tight_layout()
    fig.savefig(plot_path)
    plt.close()
    logger.info("Wrote %s", plot_path)
    prov_dict = {'caption': 'Barplot of ZEC_%s' % str(x_i), 'plot_type': 'bar'}
    write_provenance(cfg, plot_path, prov_dict)


def write_provenance(cfg, path, add_dict={}):
    """Helper function to write provenance record."""
    input_data = cfg['input_data'].values()
    provenance_record = {
        'authors': ['gier_bettina'],
        'references': ['macdougall20'],  # add sanderson24gmd later
        'ancestors': [d["filename"] for d in input_data],
    }
    for key in add_dict:
        provenance_record[key] = add_dict[key]
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(path, provenance_record)


def main(cfg):
    """Execute diagnostic."""
    # Calculate ZEC
    zec = calculate_zec(cfg)

    # Plot ZEC Timeseries
    plot_zec_timeseries(zec, cfg)

    # Calculate ZEC_X
    x = cfg['zec_x'] if isinstance(cfg['zec_x'], list) else [cfg['zec_x']]
    for x_i in x:
        zec_x = calc_zec_x(zec, x_i)
        write_zec_x(zec_x, x_i, cfg)
        plot_zec_x_bar(zec_x, x_i, cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        config = _get_default_cfg(config)
        main(config)
