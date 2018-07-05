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
                          ECS calculations.
    plot_name           : Name of plot file.
    save                : Keyword arguments for the fig.saveplot() function.
    axes_functions      : Plot appearance functions.

###############################################################################

"""


import datetime
import logging
import os

import iris
import numpy as np
from scipy import stats

import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n

logger = logging.getLogger(os.path.basename(__file__))


def calculate_ecs(cfg, datasets, variables):
    """Calculate ECS (Andrews et al. 2012)."""
    for dataset_path in \
            datasets.get_path_list(short_name=variables.tas, exp=DIFF):
        dataset = datasets.get_info(n.DATASET, dataset_path)
        data_tas = datasets.get_data(dataset_path)
        data_rtmt = datasets.get_data(short_name=variables.rtmt, exp=DIFF,
                                      dataset=dataset)
        reg = stats.linregress(data_tas, data_rtmt)
        datasets.add_dataset(variables.ecs + DIFF + dataset,
                             -reg.intercept / (2 * reg.slope),
                             short_name=variables.ecs, exp=DIFF,
                             dataset=dataset)

        # Plot ECS regression if desired
        if not (cfg[n.WRITE_PLOTS] and cfg.get('plot_ecs_regression')):
            continue
        filepath = os.path.join(cfg[n.PLOT_DIR],
                                dataset + '.' + cfg[n.OUTPUT_FILE_TYPE])
        plot_kwargs = {'linestyle': 'none',
                       'markeredgecolor': 'b',
                       'markerfacecolor': 'none',
                       'marker': 's'}

        # Regression line
        x_reg = np.linspace(-1.0, 8.0, 2)
        y_reg = reg.slope * x_reg + reg.intercept
        reg_plot_kwargs = {'color': 'k', 'linestyle': '-'}

        # Plot data
        text = 'r = {:.2f}, '.format(reg.rvalue) + \
               r'$\alpha$ = {:.2f}, '.format(-reg.slope) + \
               'F = {:.2f}, '.format(reg.intercept) + \
               'ECS = {:.2f}'.format(-reg.intercept / (2 * reg.slope))
        e.plot.scatterplot(
            [data_tas, x_reg],
            [data_rtmt, y_reg],
            filepath,
            plot_kwargs=[plot_kwargs, reg_plot_kwargs],
            save_kwargs=cfg.get('save'),
            axes_functions={'set_title': dataset,
                            'set_xlabel': (variables.tas + " / " +
                                           variables.TAS.units),
                            'set_ylabel': (variables.rtmt + " / " +
                                           variables.RTMT.units),
                            'set_xlim': [0.0, 7.0],
                            'set_ylim': [-2.0, 10.0],
                            'text': {'args': [0.05, 0.9, text],
                                     'kwargs': {'transform': 'transAxes'}}})


def plot_data(cfg, datasets, variables):
    """Plot data."""
    if not cfg[n.WRITE_PLOTS]:
        return
    filepath = os.path.join(cfg[n.PLOT_DIR],
                            cfg.get('plot_name', 'ch09_fig09-42a') + '.' +
                            cfg[n.OUTPUT_FILE_TYPE])
    x_data = []
    y_data = []
    dataset_names = []
    plot_kwargs = []
    names = datasets.get_info_list(n.DATASET, short_name=variables.ecs,
                                   exp=DIFF)
    ecs_data = datasets.get_data_list(short_name=variables.ecs, exp=DIFF)

    # Historical
    x_data.extend(ecs_data)
    y_data.extend(datasets.get_data_list(short_name=variables.tas,
                                         exp=HISTORICAL))
    dataset_names.extend(names)
    for name in names:
        plot_kwargs.append({'label': name, 'linestyle': 'none',
                            'markersize': 10})

    # piControl
    x_data.extend(ecs_data)
    y_data.extend(datasets.get_data_list(short_name=variables.tas,
                                         exp=PICONTROL_TEMP_MEAN))
    dataset_names.extend(names)
    for name in names:
        plot_kwargs.append({'label': '_' + name, 'linestyle': 'none',
                            'markersize': 6})

    # Plot data
    e.plot.multi_dataset_scatterplot(
        x_data,
        y_data,
        dataset_names,
        filepath,
        plot_kwargs=plot_kwargs,
        save_kwargs=cfg.get('save'),
        axes_functions=cfg.get('axes_functions'))
    return


def write_data(cfg, datasets, variables):
    """Write netcdf file."""
    if cfg[n.WRITE_PLOTS]:
        data_ecs = datasets.get_data_list(short_name=variables.ecs, exp=DIFF)
        models = datasets.get_info_list(n.DATASET, short_name=variables.ecs,
                                        exp=DIFF)
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
        if data.get_info(n.EXP, dataset_path) == HISTORICAL:
            cube = cube.collapsed([n.TIME], iris.analysis.MEAN)
        else:
            cube = cube.aggregated_by(n.YEAR, iris.analysis.MEAN)

        data.set_data(cube.data, dataset_path)

    ###########################################################################
    # Process data
    ###########################################################################

    # Substract piControl experiment from abrupt4xCO2 experiment and add total
    # temporal mean of piControl
    for dataset_path in data.get_path_list(exp=PICONTROL):
        dataset = data.get_info(n.DATASET, dataset_path)
        short_name = data.get_info(n.SHORT_NAME, dataset_path)
        new_data = data.get_data(dataset_path)
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
