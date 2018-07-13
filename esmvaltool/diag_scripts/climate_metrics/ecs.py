#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Calculate ECS following Andrews et al. (2012).

###############################################################################
climate_metrics/ecs.py
Author: Manuel Schlund (DLR, Germany)
CRESCENDO project
###############################################################################

Description
-----------
    Calculate the equilibrium climate sensitivity (ECS) using the regression
    method proposed by Andrews et al. (2012).

Configuration options
---------------------
    plot_regression : Switch to plot the linear regression.
    output_name     : Name of the output files.

###############################################################################

"""


import logging
import os
from collections import OrderedDict
from datetime import datetime

import iris
import numpy as np
from scipy import stats

import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n

logger = logging.getLogger(os.path.basename(__file__))


def plot_ecs_regression(cfg, dataset_name, data, variables, regression_stats):
    """Plot linear regression used to calculate ECS."""
    if not (cfg[n.WRITE_PLOTS] and cfg.get('plot_ecs_regression')):
        return
    ecs = -regression_stats.intercept / (2 * regression_stats.slope)
    filepath = os.path.join(cfg[n.PLOT_DIR],
                            dataset_name + '.' + cfg[n.OUTPUT_FILE_TYPE])

    # Regression line
    x_reg = np.linspace(-1.0, 8.0, 2)
    y_reg = regression_stats.slope * x_reg + regression_stats.intercept

    # Plot data
    text = 'r = {:.2f}, '.format(regression_stats.rvalue) + \
           r'$\alpha$ = {:.2f}, '.format(-regression_stats.slope) + \
           'F = {:.2f}, '.format(regression_stats.intercept) + \
           'ECS = {:.2f}'.format(ecs)
    e.plot.scatterplot(
        [data[0], x_reg],
        [data[1], y_reg],
        filepath,
        plot_kwargs=[{'linestyle': 'none',
                      'markeredgecolor': 'b',
                      'markerfacecolor': 'none',
                      'marker': 's'},
                     {'color': 'k',
                      'linestyle': '-'}],
        save_kwargs={
            'bbox_inches': 'tight',
            'orientation': 'landscape'},
        axes_functions={
            'set_title': dataset_name,
            'set_xlabel': 'tas / ' + variables.units('tas'),
            'set_ylabel': 'rtmt / ' + variables.units('rtmt'),
            'set_xlim': [0.0, 7.0],
            'set_ylim': [-2.0, 10.0],
            'text': {'args': [0.05, 0.9, text],
                     'kwargs': {'transform': 'transAxes'}}})

    # Write netcdf file for every plot
    if not cfg[n.WRITE_NETCDF]:
        return
    tas_coord = iris.coords.AuxCoord(data[0], **variables.iris_dict('tas'))
    attr = {'model': dataset_name,
            'regression_r_value': regression_stats.rvalue,
            'regression_slope': regression_stats.slope,
            'regression_interception': regression_stats.intercept,
            'climate_sensitivity': -regression_stats.slope,
            'ECS': ecs,
            'created_by': 'ESMValTool version {}'.format(cfg[n.VERSION]) +
                          ', diagnostic {}'.format(cfg[n.SCRIPT]),
            'creation_date': datetime.utcnow().isoformat(' ') + 'UTC'}
    cube = iris.cube.Cube(data[1],
                          attributes=attr,
                          aux_coords_and_dims=[(tas_coord, 0)],
                          **variables.iris_dict('rtmt'))
    filepath = os.path.join(cfg[n.WORK_DIR],
                            'ecs_regression_' + dataset_name + '.nc')
    iris.save(cube, filepath)
    logger.info("Writing %s", filepath)
    return


###############################################################################
# Setup diagnostic
###############################################################################

# Variables
ECS = e.Variable('ecs',
                 'equilibrium_climate_sensitivity',
                 'Change in global mean surface temperature at equilibrium '
                 'caused by a doubling of the atmospheric CO2 concentration',
                 'K')

# Experiments
PICONTROL = 'piControl'
ABRUPT4XCO2 = 'abrupt4xCO2'


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
    var.add_vars(ecs=ECS)

    # Check for tas and rtmt
    if not var.vars_available('tas', 'rtmt'):
        raise ValueError("This diagnostic needs 'tas' and 'rtmt' variables")

    ###########################################################################
    # Read data
    ###########################################################################

    # Create iris cube for each dataset and save annual means
    for dataset_path in data:
        cube = iris.load(dataset_path, var.standard_names())[0]
        cube = cube.aggregated_by(n.YEAR, iris.analysis.MEAN)
        data.set_data(cube.data, dataset_path)

    ###########################################################################
    # Process data
    ###########################################################################
    data_ecs = OrderedDict()

    for dataset_path in \
            data.get_path_list(short_name='tas', exp=PICONTROL):

        # Substract piControl experiment from abrupt4xCO2 experiment
        dataset = data.get_info(n.DATASET, dataset_path)
        data_rtmt_pic = data.get_data(short_name='rtmt', exp=PICONTROL,
                                      dataset=dataset)
        data_tas = data.get_data(short_name='tas', exp=ABRUPT4XCO2,
                                 dataset=dataset) - data.get_data(dataset_path)
        data_rtmt = data.get_data(short_name='rtmt', exp=ABRUPT4XCO2,
                                  dataset=dataset) - data_rtmt_pic

        # Perform linear regression
        reg = stats.linregress(data_tas, data_rtmt)

        # Plot ECS regression if desired
        plot_ecs_regression(cfg, dataset, [data_tas, data_rtmt], var, reg)

        # Save data
        data_ecs[dataset] = -reg.intercept / (2 * reg.slope)

    ###########################################################################
    # Write data
    ###########################################################################
    if cfg[n.WRITE_NETCDF]:
        dataset_coord = iris.coords.AuxCoord(list(data_ecs),
                                             long_name='datasets')
        attr = {'created_by': 'ESMValTool version {}'.format(cfg[n.VERSION]) +
                              ', diagnostic {}'.format(cfg[n.SCRIPT]),
                'creation_date': datetime.utcnow().isoformat(' ') + 'UTC'}
        cube = iris.cube.Cube(list(data_ecs.values()),
                              long_name=var.long_name('ecs'),
                              var_name='ecs',
                              units=var.units('ecs'),
                              aux_coords_and_dims=[(dataset_coord, 0)],
                              attributes=attr)

        # Save file
        filepath = os.path.join(cfg[n.WORK_DIR],
                                cfg.get('output_name', 'ecs') + '.nc')
        iris.save(cube, filepath)
        logger.info("Writing %s", filepath)


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
