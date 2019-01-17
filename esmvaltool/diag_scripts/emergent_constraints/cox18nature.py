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

Configuration options in recipe
-------------------------------
confidence_level : float, optional (default: 0.66)
    Confidence level for ECS error estimation.

"""

import logging
import os

import iris
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

import esmvaltool.diag_scripts.emergent_constraints as ec
from esmvaltool.diag_scripts.shared import (
    get_ancestor_file, get_plot_filename, group_metadata,
    iris_project_constraint, match_dataset_coordinates, netcdf_to_metadata,
    plot, run_diagnostic, variables_available)

logger = logging.getLogger(os.path.basename(__file__))
plt.style.use(plot.get_path_to_mpl_style())

COLOR_SMALL_LAMBDA = '#800060'
COLOR_LARGE_LAMBDA = '#009900'
(FIG, AXES) = plt.subplots()


def _get_model_color(model, lambda_cube):
    """Get color of model dependent on climate sensitivity."""
    clim_sens = lambda_cube.extract(iris.Constraint(dataset=model)).data
    if clim_sens < 1.0:
        col = COLOR_SMALL_LAMBDA
    else:
        col = COLOR_LARGE_LAMBDA
    return col


def _plot_model_point(model, psi_cube, ecs_cube, lambda_cube):
    """Plot a single model point for emergent relationship."""
    col = _get_model_color(model, lambda_cube)
    style = plot.get_dataset_style(model, 'cox18nature.yml')
    AXES.plot(
        psi_cube.extract(iris.Constraint(dataset=model)).data,
        ecs_cube.extract(iris.Constraint(dataset=model)).data,
        linestyle='none',
        marker=style['mark'],
        markeredgecolor=col,
        markerfacecolor=col,
        markersize=style['size'])


def _get_line_plot_legend():
    """Add legend for line plots."""
    color_obs = plot.get_dataset_style('OBS', 'cox18nature.yml')['color']
    handles = [
        mlines.Line2D([], [],
                      color=COLOR_SMALL_LAMBDA,
                      label=r'$\lambda < 1.0$ Wm$^{-2}$K$^{-1}$'),
        mlines.Line2D([], [],
                      color=COLOR_LARGE_LAMBDA,
                      label=r'$\lambda > 1.0$ Wm$^{-2}$K$^{-1}$'),
        mlines.Line2D([], [],
                      linestyle='none',
                      marker='o',
                      markeredgecolor=color_obs,
                      markerfacecolor=color_obs,
                      label='Observations'),
    ]
    return AXES.legend(handles=handles, loc='upper left')


def _save_fig(cfg, basename, legend=None):
    """Save matplotlib figure."""
    path = get_plot_filename(basename, cfg)
    if legend is None:
        legend = []
    else:
        legend = [legend]
    FIG.savefig(
        path,
        additional_artists=legend,
        bbox_inches='tight',
        orientation='landscape')
    logger.info("Wrote %s", path)
    AXES.cla()


def get_external_cubes(cfg):
    """Get external cubes for psi, ECS and lambda."""
    cubes = []
    for filename in ('psi.nc', 'ecs.nc', 'lambda.nc'):
        filepath = get_ancestor_file(cfg, filename)
        cube = iris.load_cube(filepath)
        cube = cube.extract(iris_project_constraint(['OBS'], cfg, negate=True))
        cubes.append(cube)
    cubes = match_dataset_coordinates(cubes)
    return (cubes[0], cubes[1], cubes[2])


def plot_temperature_anomaly(cfg, tas_cubes, lambda_cube, obs_name):
    """Plot temperature anomaly versus time."""
    if not cfg['write_plots']:
        return
    models = lambda_cube.coord('dataset').points

    # Plot lines
    base_constraint = iris.Constraint(year=lambda cell: 1961 <= cell <= 1990)
    for model in models:
        col = _get_model_color(model, lambda_cube)
        cube = tas_cubes[model]
        cube.data -= np.mean(cube.extract(base_constraint).data)
        AXES.plot(cube.coord('year').points, cube.data, color=col)
    obs_style = plot.get_dataset_style('OBS', 'cox18nature.yml')
    obs_cube = tas_cubes[obs_name]
    obs_cube.data -= np.mean(obs_cube.extract(base_constraint).data)
    AXES.plot(
        obs_cube.coord('year').points,
        obs_cube.data,
        linestyle='none',
        marker='o',
        markeredgecolor=obs_style['color'],
        markerfacecolor=obs_style['color'])

    # Plot appearance
    AXES.set_title('Simulation of global warming record')
    AXES.set_xlabel('Year')
    AXES.set_ylabel('Temperature anomaly / K')
    legend = _get_line_plot_legend()
    _save_fig(cfg, 'temperature_anomaly_{}'.format(obs_name), legend)


def plot_psi(cfg, psi_cubes, lambda_cube, obs_name):
    """Plot temperature variability metric psi versus time."""
    if not cfg['write_plots']:
        return
    models = lambda_cube.coord('dataset').points

    # Plot lines
    for model in models:
        col = _get_model_color(model, lambda_cube)
        cube = psi_cubes[model]
        AXES.plot(cube.coord('year').points, cube.data, color=col)
    obs_style = plot.get_dataset_style('OBS', 'cox18nature.yml')
    obs_cube = psi_cubes[obs_name]
    AXES.plot(
        obs_cube.coord('year').points,
        obs_cube.data,
        linestyle='none',
        marker='o',
        markeredgecolor=obs_style['color'],
        markerfacecolor=obs_style['color'])

    # Plot appearance
    AXES.set_title('Metric of variability versus time')
    AXES.set_xlabel('Year')
    AXES.set_ylabel(r'$\Psi$ / K')
    legend = _get_line_plot_legend()
    _save_fig(cfg, 'temperature_variability_metric_{}'.format(obs_name),
              legend)


def plot_emergent_relationship(cfg, psi_cube, ecs_cube, lambda_cube, obs_cube):
    """Plot emergent relationship."""
    if not cfg['write_plots']:
        return
    models = psi_cube.coord('dataset').points
    obs_mean = np.mean(obs_cube.data)
    obs_std = np.std(obs_cube.data)

    # Calculate regression line
    lines = ec.regression_line(psi_cube.data, ecs_cube.data)
    logger.info("Found emergent relationship with slope %.2f (r = %.2f)",
                lines['slope'], lines['rvalue'])

    # Plot points
    for model in models:
        _plot_model_point(model, psi_cube, ecs_cube, lambda_cube)

    # Plot lines
    AXES.set_xlim(auto=False)
    AXES.set_ylim(auto=False)
    AXES.plot(
        lines['x'],
        lines['y_best_estim'],
        color='black',
        linestyle='dashdot',
        label='Linear regression')
    AXES.plot(
        lines['x'], lines['y_minus_err'], color='black', linestyle='dashed')
    AXES.plot(
        lines['x'], lines['y_plus_err'], color='black', linestyle='dashed')
    AXES.axvline(
        obs_mean,
        color='blue',
        linestyle='dashdot',
        label='Observational constraint')
    AXES.axvline(obs_mean - obs_std, color='blue', linestyle='dashed')
    AXES.axvline(obs_mean + obs_std, color='blue', linestyle='dashed')

    # Plot appearance
    AXES.set_title('Emergent relationship fit')
    AXES.set_xlabel(r'$\Psi$ / K')
    AXES.set_ylabel('ECS / K')
    legend = AXES.legend(loc='upper left')

    # Save plot
    _save_fig(
        cfg, 'emergent_relationship_{}'.format(obs_cube.attributes['dataset']),
        legend)


def plot_pdf(cfg, psi_cube, ecs_cube, obs_cube):
    """Plot probability density function of ECS."""
    if not cfg['write_plots']:
        return
    obs_mean = np.mean(obs_cube.data)
    obs_std = np.std(obs_cube.data)

    # Calculate PDF
    (ecs_lin, ecs_pdf) = ec.gaussian_pdf(psi_cube.data, ecs_cube.data,
                                         obs_mean, obs_std)

    # Plot
    AXES.plot(
        ecs_lin,
        ecs_pdf,
        color='black',
        linewidth=2.0,
        label='Emergent constraint')
    AXES.hist(
        ecs_cube.data,
        bins=6,
        range=(2.0, 5.0),
        density=True,
        color='orange',
        label='CMIP5 models')

    # Plot appearance
    AXES.set_title('PDF of emergent constraint')
    AXES.set_xlabel('ECS / K')
    AXES.set_ylabel('Probability density')
    legend = AXES.legend(loc='upper left')

    # Save plot
    _save_fig(cfg, 'pdf_{}'.format(obs_cube.attributes['dataset']), legend)


def plot_cdf(cfg, psi_cube, ecs_cube, obs_cube):
    """Plot cumulative distribution function of ECS."""
    if not cfg['write_plots']:
        return
    obs_mean = np.mean(obs_cube.data)
    obs_std = np.std(obs_cube.data)
    confidence_level = cfg.get('confidence_level', 0.66)
    conf_low = (1.0 - confidence_level) / 2.0
    conf_high = (1.0 + confidence_level) / 2.0

    # Calculate PDF and CDF
    (ecs_lin, ecs_pdf) = ec.gaussian_pdf(psi_cube.data, ecs_cube.data,
                                         obs_mean, obs_std)
    ecs_cdf = ec.cdf(ecs_lin, ecs_pdf)

    # Plot
    AXES.plot(
        ecs_lin,
        ecs_cdf,
        color='black',
        linewidth=2.0,
        label='Emergent constraint')
    AXES.hist(
        ecs_cube.data,
        bins=6,
        range=(2.0, 5.0),
        cumulative=True,
        density=True,
        color='orange',
        label='CMIP5 models')
    AXES.axhline(conf_low, color='black', linestyle='dashdot')
    AXES.axhline(conf_high, color='black', linestyle='dashdot')

    # Plot appearance
    AXES.set_title('CDF of emergent constraint')
    AXES.set_xlabel('ECS / K')
    AXES.set_ylabel('CDF')
    legend = AXES.legend(loc='upper left')

    # Save plot
    _save_fig(cfg, 'cdf_{}'.format(obs_cube.attributes['dataset']), legend)


def get_ecs_range(cfg, psi_cube, ecs_cube, obs_cube):
    """Get constrained ecs range."""
    confidence_level = cfg.get('confidence_level', 0.66)
    conf_low = (1.0 - confidence_level) / 2.0
    conf_high = (1.0 + confidence_level) / 2.0

    # Calculate PDF and CDF
    (ecs_lin, ecs_pdf) = ec.gaussian_pdf(psi_cube.data, ecs_cube.data,
                                         np.mean(obs_cube.data),
                                         np.std(obs_cube.data))
    ecs_cdf = ec.cdf(ecs_lin, ecs_pdf)

    # Calculate constrained ECS range
    ecs_mean = ecs_lin[np.argmax(ecs_pdf)]
    ecs_index_range = np.where((ecs_cdf >= conf_low) &
                               (ecs_cdf <= conf_high))[0]
    ecs_range = ecs_lin[ecs_index_range]
    ecs_low = min(ecs_range)
    ecs_high = max(ecs_range)
    return (ecs_mean, ecs_low, ecs_high)


def main(cfg):
    """Run the diagnostic."""
    if not variables_available(cfg, ['tas']):
        raise ValueError("This diagnostic needs 'tas' variable")

    # Get tas data
    tas_cubes = {}
    tas_obs = []
    for (dataset, [data]) in group_metadata(cfg['input_data'].values(),
                                            'dataset').items():
        cube = iris.load_cube(data['filename'])
        cube = cube.aggregated_by('year', iris.analysis.MEAN)
        tas_cubes[dataset] = cube
        if data['project'] == 'OBS':
            tas_obs.append(dataset)

    # Get time-dependent psi data
    psi_cubes = {}
    psi_obs = []
    psi_data = netcdf_to_metadata(cfg, pattern='psi_*.nc')
    for (dataset, [data]) in group_metadata(psi_data, 'dataset').items():
        cube = iris.load_cube(data['filename'])
        cube = cube.aggregated_by('year', iris.analysis.MEAN)
        psi_cubes[dataset] = cube
        if data['project'] == 'OBS':
            psi_obs.append(dataset)

    # Get psi, ECS and psi for models
    (psi_cube, ecs_cube, lambda_cube) = get_external_cubes(cfg)

    # Plots
    for obs_name in tas_obs:
        logger.info("Observation for tas: %s", obs_name)
        plot_temperature_anomaly(cfg, tas_cubes, lambda_cube, obs_name)
    for obs_name in psi_obs:
        logger.info("Observation for psi: %s", obs_name)
        plot_psi(cfg, psi_cubes, lambda_cube, obs_name)
        obs_cube = psi_cubes[obs_name]
        plot_emergent_relationship(cfg, psi_cube, ecs_cube, lambda_cube,
                                   obs_cube)
        plot_pdf(cfg, psi_cube, ecs_cube, obs_cube)
        plot_cdf(cfg, psi_cube, ecs_cube, obs_cube)

        # Print ECS range
        ecs_range = get_ecs_range(cfg, psi_cube, ecs_cube, obs_cube)
        logger.info("Observational constraint: Ψ = (%.2f ± %.2f) K",
                    np.mean(obs_cube.data), np.std(obs_cube.data))
        logger.info(
            "Constrained ECS range: (%.2f - %.2f) K with best "
            "estimate %.2f K", ecs_range[1], ecs_range[2], ecs_range[0])


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
        plt.close()
