#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to reproduce Cox et al. (2018).

Description
-----------
Plot equilibrium climate sensitivity ECS vs. temperature variability metric psi
to establish an emergent relationship for ECS.

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recips
-------------------------------
confidence_level : float, optional (default: 0.66)
    Confidence level for ECS error estimation.

"""

import logging
import os

import iris
import matplotlib  # noqa
import matplotlib.pyplot as plt  # noqa
import numpy as np
from iris import Constraint

import esmvaltool.diag_scripts.emergent_constraints as ec
from esmvaltool.diag_scripts.shared import (
    extract_variables, get_all_ancestor_files, get_file_from_ancestors,
    group_metadata, plot, run_diagnostic, save_iris_cube, variables_available)

matplotlib.use('Agg')  # noqa

logger = logging.getLogger(os.path.basename(__file__))
plt.style.use(plot.get_path_to_mpl_style())

COLOR_SMALL_LAMBDA = '#800060'
COLOR_LARGE_LAMBDA = '#009900'


def _get_model_color(model, lambda_cube):
    """Get color of model dependent on climate sensitivity."""
    clim_sens = lambda_cube.extract(Constraint(dataset=model)).data
    if clim_sens < 1.0:
        col = '#800060'
    else:
        col = '#009900'
    return col


def _plot_model_point(axes, model, psi_cube, ecs_cube, lambda_cube):
    """Plot a single model point for emergent relationship."""
    col = _get_model_color(model, lambda_cube)
    style = plot.get_dataset_style(model, 'cox18nature.yml')
    axes.plot(
        psi_cube.extract(Constraint(dataset=model)).data,
        ecs_cube.extract(Constraint(dataset=model)).data,
        linestyle='none',
        marker=style['mark'],
        markeredgecolor=col,
        markerfacecolor=col,
        markersize=style['size'])


def plot_temperature_anomaly(cfg, tas_cubes, lambda_cube):
    """Plot temperature anomaly versus time."""
    if not cfg['write_plots']:
        return
    (fig, axes) = plt.subplots()
    models = list(tas_cubes.keys())


def plot_emergent_relationship(cfg, psi_cube, ecs_cube, lambda_cube, obs_cube):
    """Plot emergent relationship."""
    if not cfg['write_plots']:
        return
    (fig, axes) = plt.subplots()
    models = psi_cube.coord('dataset').points
    obs_mean = np.mean(obs_cube.data)
    obs_std = np.std(obs_cube.data)

    # Calculate regression line
    lines = ec.regression_line(psi_cube.data, ecs_cube.data)

    # Plot points
    for model in models:
        _plot_model_point(axes, model, psi_cube, ecs_cube, lambda_cube)

    # Plot lines
    axes.set_xlim(auto=False)
    axes.set_ylim(auto=False)
    axes.plot(
        lines['x'],
        lines['y_best_estim'],
        color='black',
        linestyle='dashdot',
        label='Linear regression')
    axes.plot(
        lines['x'], lines['y_minus_err'], color='black', linestyle='dashed')
    axes.plot(
        lines['x'], lines['y_plus_err'], color='black', linestyle='dashed')
    axes.axvline(
        obs_mean,
        color='blue',
        linestyle='dashdot',
        label='Observational constraint')
    axes.axvline(obs_mean - obs_std, color='blue', linestyle='dashed')
    axes.axvline(obs_mean + obs_std, color='blue', linestyle='dashed')

    # Plot appearance
    axes.set_title('Emergent relationship fit')
    axes.set_xlabel(r'$\Psi$ / K')
    axes.set_ylabel('ECS / K')
    legend = axes.legend(loc='upper left')

    # Save plot
    path = os.path.join(
        cfg['plot_dir'],
        'emergent_relationship.{}'.format(cfg['output_file_type']))
    fig.savefig(
        path,
        additional_artists=[legend],
        bbox_inches='tight',
        orientation='landscape')
    logger.info("Wrote %s", path)
    plt.close()


def main(cfg):
    """Run the diagnostic."""
    input_data = cfg['input_data'].values()

    # Check if tas is available
    if not variables_available(cfg, ['tas']):
        raise ValueError("This diagnostic needs 'tas' variable")

    # Get tas data
    tas_cubes = {}
    for (dataset, [data]) in group_metadata(input_data, 'dataset').items():
        cube = iris.load_cube(data['filename'])
        cube = cube.aggregated_by
        tas_cubes[dataset] = cube

    # Get paths to data (ECS, lambda and psi)
    ecs_filepath = get_file_from_ancestors(cfg, 'ecs.nc')
    lambda_filepath = get_file_from_ancestors(cfg, 'lambda.nc')
    psi_filepath = get_file_from_ancestors(cfg, 'psi.nc')
    psi_files = get_all_ancestor_files(cfg, pattern='psi_*.nc')

    # Load cubes
    psi_cube_all = iris.load_cube(psi_filepath)
    psi_cube = psi_cube_all.extract(ec.iris_constraint_no_obs(cfg))
    ecs_cube = iris.load_cube(ecs_filepath)
    lambda_cube = iris.load_cube(lambda_filepath)
    ec.check_dataset_dimensions(psi_cube, ecs_cube, lambda_cube)
    psi_cubes = {}
    obs_cubes = {}
    for path in psi_files:
        cube = iris.load_cube(path)
        if cube.attributes['project'] == 'OBS':
            obs_cubes[cube.attributes['dataset']] = cube
        else:
            psi_cubes[cube.attributes['dataset']] = cube

    # Plots
    for obs_cube in obs_cubes.values():
        plot_temperature_anomaly(cfg, tas_cubes, lambda_cube)
        plot_emergent_relationship(cfg, psi_cube, ecs_cube, lambda_cube,
                                   obs_cube)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
