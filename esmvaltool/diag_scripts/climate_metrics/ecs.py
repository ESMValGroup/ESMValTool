#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to calculate ECS following Andrews et al. (2012).

Description
-----------
Calculate the equilibrium climate sensitivity (ECS) using the regression method
proposed by Andrews et al. (2012).

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recips
-------------------------------
plot_ecs_regression : bool, optional (default: False)
    Plot the linear regression graph.
read_external_file : str, optional
    Read ECS from external file.

"""

import logging
import os
from pprint import pformat

import iris
import numpy as np
import yaml
from scipy import stats

from esmvaltool.diag_scripts.shared import (
    extract_variables, group_metadata, plot, run_diagnostic, save_iris_cube,
    save_scalar_data, select_metadata, variables_available)

logger = logging.getLogger(os.path.basename(__file__))


def read_external_file(cfg):
    """Read external file to get ECS."""
    ecs = {}
    clim_sens = {}
    if not cfg.get('read_external_file'):
        return (ecs, clim_sens)
    base_dir = os.path.dirname(__file__)
    filepath = os.path.join(base_dir, cfg['read_external_file'])
    if os.path.isfile(filepath):
        with open(filepath, 'r') as infile:
            external_data = yaml.safe_load(infile)
    else:
        logger.error("Desired external file %s does not exist", filepath)
        return (ecs, clim_sens)
    ecs = external_data.get('ecs', {})
    clim_sens = external_data('climate_sensitivity', {})
    logger.info("External file %s", filepath)
    logger.info("Found ECS:")
    logger.info("%s", pformat(ecs))
    logger.info("Found climate sensitivities:")
    logger.info("%s", pformat(clim_sens))
    return (ecs, clim_sens)


def plot_ecs_regression(cfg, dataset_name, tas_cube, rtmt_cube,
                        regression_stats):
    """Plot linear regression used to calculate ECS."""
    if not (cfg['write_plots'] and cfg.get('plot_ecs_regression')):
        return
    ecs = -regression_stats.intercept / (2 * regression_stats.slope)
    filepath = os.path.join(cfg['plot_dir'],
                            dataset_name + '.' + cfg['output_file_type'])

    # Regression line
    x_reg = np.linspace(-1.0, 8.0, 2)
    y_reg = regression_stats.slope * x_reg + regression_stats.intercept

    # Plot data
    text = 'r = {:.2f}, '.format(regression_stats.rvalue) + \
           r'$\lambda$ = {:.2f}, '.format(-regression_stats.slope) + \
           'F = {:.2f}, '.format(regression_stats.intercept) + \
           'ECS = {:.2f}'.format(ecs)
    plot.scatterplot(
        [tas_cube.data, x_reg],
        [rtmt_cube.data, y_reg],
        filepath,
        plot_kwargs=[{
            'linestyle': 'none',
            'markeredgecolor': 'b',
            'markerfacecolor': 'none',
            'marker': 's',
        }, {
            'color': 'k',
            'linestyle': '-',
        }],
        save_kwargs={
            'bbox_inches': 'tight',
            'orientation': 'landscape',
        },
        axes_functions={
            'set_title': dataset_name,
            'set_xlabel': 'tas / ' + tas_cube.units.name,
            'set_ylabel': 'rtmt / ' + rtmt_cube.units.name,
            'set_xlim': [0.0, 7.0],
            'set_ylim': [-2.0, 10.0],
            'text': {
                'args': [0.05, 0.9, text],
                'kwargs': {
                    'transform': 'transAxes'
                },
            },
        },
    )

    # Write netcdf file for every plot
    if not cfg['write_netcdf']:
        return
    tas_coord = iris.coords.AuxCoord(
        tas_cube.data,
        **extract_variables(cfg, as_iris=True)['tas'])
    attr = {
        'model': dataset_name,
        'regression_r_value': regression_stats.rvalue,
        'regression_slope': regression_stats.slope,
        'regression_interception': regression_stats.intercept,
        'climate_sensitivity': -regression_stats.slope,
        'ECS': ecs,
    }
    cube = iris.cube.Cube(
        rtmt_cube.data,
        attributes=attr,
        aux_coords_and_dims=[(tas_coord, 0)],
        **extract_variables(cfg, as_iris=True)['rtmt'])
    filepath = os.path.join(cfg['work_dir'],
                            'ecs_regression_' + dataset_name + '.nc')
    save_iris_cube(cube, filepath, cfg)
    return


def main(cfg):
    """Run the diagnostic."""
    input_data = cfg['input_data'].values()

    # Read external file if desired
    if cfg.get('read_external_file'):
        (ecs, clim_sens) = read_external_file(cfg)
    else:
        if not variables_available(cfg, ['tas', 'rtmt']):
            raise ValueError("This diagnostic needs 'tas' and 'rtmt' "
                             "variables if 'read_external_file' is not given")
        ecs = {}
        clim_sens = {}

    # Read data
    tas_data = select_metadata(input_data, short_name='tas')
    rtmt_data = select_metadata(input_data, short_name='rtmt')

    # Iterate over all datasets and save ECS
    for (dataset, data) in group_metadata(tas_data, 'dataset').items():
        logger.info("Processing %s", dataset)
        paths = {
            'tas_4x':
            select_metadata(data, exp='abrupt4xCO2'),
            'tas_pi':
            select_metadata(data, exp='piControl'),
            'rtmt_4x':
            select_metadata(rtmt_data, dataset=dataset, exp='abrupt4xCO2'),
            'rtmt_pi':
            select_metadata(rtmt_data, dataset=dataset, exp='piControl'),
        }
        cubes = {}
        for (key, path) in paths.items():
            cube = iris.load_cube(path[0]['filename'])
            cube = cube.aggregated_by('year', iris.analysis.MEAN)
            cubes[key] = cube

        # Substract piControl run from abrupt4xCO2 experiment
        cubes['tas_4x'].data -= cubes['tas_pi'].data
        cubes['rtmt_4x'].data -= cubes['rtmt_pi'].data

        # Perform linear regression
        reg = stats.linregress(cubes['tas_4x'].data, cubes['rtmt_4x'].data)

        # Plot ECS regression if desired
        plot_ecs_regression(cfg, dataset, cubes['tas_4x'], cubes['rtmt_4x'],
                            reg)

        # Save data
        if cfg.get('read_external_file') and dataset in ecs:
            logger.info(
                "Overwriting external given ECS and climate "
                "sensitivity for %s", dataset)
        ecs[dataset] = -reg.intercept / (2 * reg.slope)
        clim_sens[dataset] = -reg.slope

    # Write data
    path = os.path.join(cfg['work_dir'], 'ecs.nc')
    var_attrs = {
        'short_name': 'ecs',
        'standard_name': 'equilibrium_climate_sensitivity',
        'long_name': 'Equilibrium Climate Sensitivity (ECS)',
        'units': cubes['tas_4x'].units,
    }
    save_scalar_data(ecs, path, cfg, var_attrs)
    path = os.path.join(cfg['work_dir'], 'lambda.nc')
    var_attrs = {
        'short_name': 'lambda',
        'standard_name': 'climate_sensitivity',
        'long_name': 'Climate Sensitivity',
        'units': cubes['rtmt_4x'].units / cubes['tas_4x'].units,
    }
    save_scalar_data(clim_sens, path, cfg, var_attrs)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
