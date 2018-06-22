#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Plot figure 9.42a of IPCC AR5 chapter 9.

###############################################################################
ipcc_ar5/ch09_fig09-42a.py
Author: Manuel Schlund (DLR, Germany)
CRESCENDO project
###############################################################################

Description
-----------
    Calculate and plot the equilibrium climate sensitivity (ECS) vs. the global
    mean surface temperature (GMSAT) for several CMIP5 models (see IPCC AR5 WG1
    ch. 9, fig. 9.42a).

Configuration options
---------------------
    plot_ecs_regression : Switch to plot the linear regressions needed for the
                          ECS calculations

###############################################################################

"""


import datetime
import logging
import os

from scipy import stats
import iris
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n

logger = logging.getLogger(os.path.basename(__file__))


def calculate_ecs(cfg, datasets, variables):
    """Calculate ECS (Andrews et al. 2012)."""
    style_file = e.plot.get_path_to_mpl_style('default.mplstyle')
    plt.style.use(style_file)
    fig, axes = plt.subplots()

    # Iterate through datasets
    for dataset_path in \
            datasets.get_path_list(short_name=variables.tas, exp=DIFF):
        dataset = datasets.get_dataset(dataset_path)
        data_tas = datasets.get_data(dataset_path=dataset_path)
        data_rtmt = datasets.get_data(short_name=variables.rtmt, exp=DIFF,
                                      dataset=dataset)
        reg_stats = stats.linregress(data_tas, data_rtmt)
        data_ecs = -reg_stats.intercept / (2 * reg_stats.slope)
        datasets.add_dataset(variables.ecs + DIFF + dataset,
                             data_ecs, short_name=variables.ecs, exp=DIFF,
                             dataset=dataset)

        # Plot ECS regression if desired
        if not (cfg[n.WRITE_PLOTS] and cfg.get('plot_ecs_regression')):
            continue

        # Plot data
        axes.plot(data_tas, data_rtmt, linestyle='none',
                  markeredgecolor='b', markerfacecolor='none',
                  marker='s')

        # Plot regerssion line
        x_reg = np.linspace(-1.0, 8.0, 2)
        y_reg = reg_stats.slope * x_reg + reg_stats.intercept
        axes.plot(x_reg, y_reg, color='k', linestyle='-')

        # Options
        axes.set_title(dataset)
        axes.set_xlabel(variables.TAS.standard_name + " / " +
                        variables.TAS.units)
        axes.set_ylabel(variables.RTMT.standard_name + " / " +
                        variables.RTMT.units)
        axes.set_xlim(0.0, 7.0)
        axes.set_ylim(-2.0, 10.0)
        axes.axhline(linestyle='dotted', c='black')
        axes.text(0.05, 0.05,
                  "r = {:.2f}".format(reg_stats.rvalue),
                  transform=axes.transAxes)
        axes.text(0.05, 0.9,
                  r"$\alpha$ = {:.2f},  ".format(-reg_stats.slope) +
                  "F = {:.2f},  ".format(reg_stats.intercept) +
                  "ECS = {:.2f}".format(data_ecs),
                  transform=axes.transAxes)

        # Save plot
        filepath = os.path.join(cfg[n.PLOT_DIR],
                                dataset + '.' + cfg[n.OUTPUT_FILE_TYPE])
        fig.savefig(filepath, bbox_inches='tight', orientation='landscape')
        logger.info("Writing %s", filepath)
        axes.cla()
    plt.close()


def plot_data(cfg, datasets, variables):
    """Plot data."""
    if not cfg[n.WRITE_PLOTS]:
        return None

    # Setup matplotlib
    style_file = e.plot.get_path_to_mpl_style('default.mplstyle')
    plt.style.use(style_file)
    fig, axes = plt.subplots()
    for dataset_path in \
            datasets.get_path_list(short_name=variables.ecs, exp=DIFF):
        dataset = datasets.get_dataset(dataset_path)
        data_ecs = datasets.get_data(dataset_path=dataset_path)
        data_tas_hist = datasets.get_data(short_name=variables.tas,
                                          exp=HISTORICAL, dataset=dataset)
        data_tas_pic = datasets.get_data(short_name=variables.tas,
                                         exp=PICONTROL_TEMP_MEAN,
                                         dataset=dataset)
        style = e.plot.get_dataset_style(dataset)

        # Plot
        axes.plot(data_ecs, data_tas_hist, linestyle='none',
                  markeredgecolor=style['color'],
                  markerfacecolor=style['facecolor'],
                  marker=style['mark'], markersize=10, label=dataset)
        axes.plot(data_ecs, data_tas_pic, linestyle='none',
                  markeredgecolor=style['color'],
                  markerfacecolor=style['facecolor'],
                  marker=style['mark'], markersize=6, label='_' + dataset)

    # Options
    axes.set_title("GMSAT vs. ECS for CMIP5 datasets")
    axes.set_xlabel(variables.ECS.standard_name + " / " +
                    variables.ECS.units)
    axes.set_ylabel(variables.TAS.standard_name + " / " +
                    variables.TAS.units)
    axes.set_xlim(1.5, 5.0)
    legend = axes.legend(loc='center left',
                         bbox_to_anchor=(1.05, 0.5), borderaxespad=0.0,
                         ncol=2)

    # Save plot
    filename = 'ch09_fig09-42a.' + cfg[n.OUTPUT_FILE_TYPE]
    filepath = os.path.join(cfg[n.PLOT_DIR], filename)
    fig.savefig(filepath, additional_artists=[legend],
                bbox_inches='tight', orientation='landscape')
    logger.info("Writing %s", filepath)
    axes.cla()
    plt.close()
    return None


def write_data(cfg, datasets, variables):
    """Write netcdf file."""
    if cfg[n.WRITE_PLOTS]:
        data_ecs = datasets.get_data_list(short_name=variables.ecs, exp=DIFF)
        models = [dataset_info[n.DATASET] for dataset_info in
                  datasets.get_dataset_info_list(short_name=variables.ecs,
                                                 exp=DIFF)]
        dataset_coord = iris.coords.AuxCoord(models, long_name='models')
        time_now = datetime.datetime.utcnow().isoformat(' ') + 'UTC'
        attr = {'created_by': 'ESMValTool version {}'.format(cfg[n.VERSION]) +
                              ', diagnostic {}'.format(cfg[n.SCRIPT]),
                'creation_date': time_now}
        cube = iris.cube.Cube(data_ecs, long_name=variables.ECS.long_name,
                              var_name=variables.ecs,
                              units=variables.ECS.units, attributes=attr)
        cube.add_aux_coord(dataset_coord, 0)

        # Save file
        filepath = os.path.join(cfg[n.WORK_DIR], variables.ecs + '.nc')
        iris.save(cube, filepath)
        logger.info("Writing %s", filepath)


###########################################################################
# Setup diagnostic
###########################################################################

# Variables
ECS = e.Variable('ecs', 'ecs', 'equilibrium climate sensitivity', 'K')

# Experiments
PICONTROL = 'piControl'
HISTORICAL = 'historical'
ABRUPT4XCO2 = 'abrupt4xCO2'
DIFF = 'difference of abrupt4xCO2 and piControl'
PICONTROL_TEMP_MEAN = 'total temporal mean of piControl'


def main(cfg):
    """Run the diagnostic.

    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the recipe.

    """
    ###########################################################################
    # Read recipe data
    ###########################################################################

    # Dataset data containers
    data = e.Datasets(cfg)
    logging.debug("Found datasets in recipe:\n%s", data)

    # Variables
    var = e.Variables(cfg)
    logging.debug("Found variables in recipe:\n%s", var)
    var.add_var(ecs=ECS)

    ###########################################################################
    # Read data
    ###########################################################################

    # Create iris cube for each dataset
    for dataset_path in data:
        cube = iris.load(dataset_path, var.standard_names())[0]

        # Global mean
        for coord in [cube.coord(n.LAT), cube.coord(n.LON)]:
            if not coord.has_bounds():
                coord.guess_bounds()
        area_weights = iris.analysis.cartography.area_weights(cube)
        cube = cube.collapsed([n.LAT, n.LON], iris.analysis.MEAN,
                              weights=area_weights)

        # Historical: total temporal mean; else: annual mean
        if data.get_exp(dataset_path) == HISTORICAL:
            cube = cube.collapsed([n.TIME], iris.analysis.MEAN)
        else:
            cube = cube.aggregated_by(n.YEAR, iris.analysis.MEAN)

        data.set_data(cube.data, dataset_path=dataset_path)

    ###########################################################################
    # Process data
    ###########################################################################

    # Substract piControl experiment from abrupt4xCO2 experiment and add total
    # temporal mean of piControl
    for dataset_path in data.get_path_list(exp=PICONTROL):
        dataset = data.get_dataset(dataset_path)
        short_name = data.get_short_name(dataset_path)
        new_data = data.get_data(dataset_path=dataset_path)
        data_diff = data.get_data(short_name=short_name, exp=ABRUPT4XCO2,
                                  dataset=dataset) - new_data
        data.add_dataset(short_name + DIFF + dataset, data_diff,
                         short_name=short_name, exp=DIFF, dataset=dataset)
        data.add_dataset(short_name + PICONTROL_TEMP_MEAN + dataset,
                         np.mean(new_data), short_name=short_name,
                         exp=PICONTROL_TEMP_MEAN, dataset=dataset)

    # Calculate ECS (cf. Andrews et al. 2015)
    calculate_ecs(cfg, data, var)

    ###########################################################################
    # Plot data
    ###########################################################################

    plot_data(cfg, data, var)

    ###########################################################################
    # Write nc file
    ###########################################################################

    write_data(cfg, data, var)


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
