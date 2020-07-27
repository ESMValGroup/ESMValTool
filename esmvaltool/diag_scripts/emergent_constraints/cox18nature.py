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
import iris.coord_categorisation
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

import esmvaltool.diag_scripts.emergent_constraints as ec
import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            get_plot_filename, group_metadata,
                                            io, plot, run_diagnostic,
                                            select_metadata)

logger = logging.getLogger(os.path.basename(__file__))
plt.style.use(plot.get_path_to_mpl_style())

COLOR_SMALL_LAMBDA = '#800060'
COLOR_LARGE_LAMBDA = '#009900'
(FIG, AXES) = plt.subplots()

ECS_ATTRS = {
    'short_name': 'ecs',
    'long_name': 'Equilibrium Climate Sensitivity (ECS)',
    'units': 'K',
}
TASA_ATTRS = {
    'short_name': 'tasa',
    'long_name': 'Near-Surface Air Temperature Anomaly',
    'units': 'K',
}
PSI_ATTRS = {
    'short_name': 'psi',
    'long_name': 'Temperature variability metric',
    'units': 'K',
}


def _get_ancestor_files(cfg, obs_name, projects=None):
    """Get ancestor files for provenance."""
    if projects is None:
        projects = _get_project(cfg)
    if isinstance(projects, str):
        projects = [projects]
    datasets = []
    for project in projects:
        datasets.extend(
            select_metadata(cfg['input_data'].values(), project=project))
    datasets.extend(
        select_metadata(cfg['input_data'].values(), dataset=obs_name))
    return [d['filename'] for d in datasets]


def _get_model_color(model, lambda_cube):
    """Get color of model dependent on climate feedback parameter."""
    clim_sens = lambda_cube.extract(iris.Constraint(dataset=model)).data
    if clim_sens < 1.0:
        col = COLOR_SMALL_LAMBDA
    else:
        col = COLOR_LARGE_LAMBDA
    return col


def _plot_model_point(model, psi_cube, ecs_cube, lambda_cube):
    """Plot a single model point for emergent relationship."""
    col = _get_model_color(model, lambda_cube)
    style = plot.get_dataset_style(model, 'cox18nature')
    AXES.plot(psi_cube.extract(iris.Constraint(dataset=model)).data,
              ecs_cube.extract(iris.Constraint(dataset=model)).data,
              linestyle='none',
              marker=style['mark'],
              markeredgecolor=col,
              markerfacecolor=col,
              markersize=style['size'])


def _get_line_plot_legend():
    """Add legend for line plots."""
    color_obs = plot.get_dataset_style('OBS', 'cox18nature')['color']
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


def _get_project(cfg):
    """Extract project from cfg."""
    input_data = cfg['input_data'].values()
    projects = list(group_metadata(input_data, 'project').keys())
    projects = [p for p in projects if 'obs' not in p.lower()]
    if len(projects) == 1:
        return projects[0]
    return projects


def _save_fig(cfg, basename, legend=None):
    """Save matplotlib figure."""
    path = get_plot_filename(basename, cfg)
    if legend is None:
        legend = []
    else:
        legend = [legend]
    FIG.savefig(path,
                additional_artists=legend,
                bbox_inches='tight',
                orientation='landscape')
    logger.info("Wrote %s", path)
    AXES.cla()
    return path


def get_external_cubes(cfg):
    """Get external cubes for psi, ECS and lambda."""
    cubes = iris.cube.CubeList()
    for filename in ('psi.nc', 'ecs.nc', 'lambda.nc'):
        filepath = io.get_ancestor_file(cfg, filename)
        cube = iris.load_cube(filepath)
        cube = cube.extract(
            ih.iris_project_constraint(['OBS'], cfg, negate=True))
        cubes.append(cube)
    cubes = ih.intersect_dataset_coordinates(cubes)
    return (cubes[0], cubes[1], cubes[2])


def get_provenance_record(caption, statistics, plot_type, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'ancestors': ancestor_files,
        'authors': ['schlund_manuel'],
        'caption': caption,
        'domains': ['global'],
        'plot_type': plot_type,
        'realms': ['atmos'],
        'references': ['cox18nature'],
        'statistics': statistics,
        'themes': ['EC'],
    }
    return record


def get_psi(cfg):
    """Get time-dependent ``psi`` data."""
    psi_cubes = {}
    psi_obs = []
    for (dataset, [data]) in group_metadata(
            io.netcdf_to_metadata(cfg, pattern='psi_*.nc'), 'dataset').items():
        cube = iris.load_cube(data['filename'])
        cube = cube.aggregated_by('year', iris.analysis.MEAN)
        psi_cubes[dataset] = cube
        if data['project'] == 'OBS':
            psi_obs.append(dataset)
    return (psi_cubes, psi_obs)


def get_tas(input_data):
    """Get time-dependent ``tas`` data."""
    tas_cubes = {}
    tas_obs = []
    for (dataset, [data]) in group_metadata(input_data, 'dataset').items():
        cube = iris.load_cube(data['filename'])
        iris.coord_categorisation.add_year(cube, 'time')
        cube = cube.aggregated_by('year', iris.analysis.MEAN)
        tas_cubes[dataset] = cube
        if data['project'] == 'OBS':
            tas_obs.append(dataset)
    return (tas_cubes, tas_obs)


def plot_temperature_anomaly(cfg, tas_cubes, lambda_cube, obs_name):
    """Plot temperature anomaly versus time."""
    for cube in tas_cubes.values():
        cube.data -= np.mean(
            cube.extract(
                iris.Constraint(year=lambda cell: 1961 <= cell <= 1990)).data)

    # Save netcdf file and provenance
    filename = 'temperature_anomaly_{}'.format(obs_name)
    netcdf_path = get_diagnostic_filename(filename, cfg)
    io.save_1d_data(tas_cubes, netcdf_path, 'year', TASA_ATTRS)
    project = _get_project(cfg)
    provenance_record = get_provenance_record(
        "Simulated change in global temperature from {} models (coloured "
        "lines), compared to the global temperature anomaly from the {} "
        "dataset (black dots). The anomalies are relative to a baseline "
        "period of 1961-1990.".format(project, obs_name), ['anomaly'],
        ['times'], _get_ancestor_files(cfg, obs_name))

    # Plot
    if cfg['write_plots']:
        models = lambda_cube.coord('dataset').points

        # Plot lines
        for model in models:
            cube = tas_cubes[model]
            AXES.plot(cube.coord('year').points,
                      cube.data,
                      color=_get_model_color(model, lambda_cube))
        obs_style = plot.get_dataset_style('OBS', 'cox18nature')
        obs_cube = tas_cubes[obs_name]
        AXES.plot(obs_cube.coord('year').points,
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

        # Save plot
        provenance_record['plot_file'] = _save_fig(cfg, filename, legend)

    # Write provenance
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)


def plot_psi(cfg, psi_cubes, lambda_cube, obs_name):
    """Plot temperature variability metric psi versus time."""
    filename = 'temperature_variability_metric_{}'.format(obs_name)
    netcdf_path = get_diagnostic_filename(filename, cfg)
    io.save_1d_data(psi_cubes, netcdf_path, 'year', PSI_ATTRS)
    project = _get_project(cfg)
    provenance_record = get_provenance_record(
        "Psi metric of variability versus time, from the {0} models "
        "(coloured lines), and the {1} observational data (black circles). "
        "The psi values are calculated for windows of width {2} yr, after "
        "linear de-trending in each window. These {2}-yr windows are shown "
        "for different end times.".format(project, obs_name,
                                          cfg.get('window_length', 55)),
        ['corr', 'var'], ['times'], _get_ancestor_files(cfg, obs_name))

    # Plot
    if cfg['write_plots']:
        models = lambda_cube.coord('dataset').points

        # Plot lines
        for model in models:
            cube = psi_cubes[model]
            AXES.plot(cube.coord('year').points,
                      cube.data,
                      color=_get_model_color(model, lambda_cube))
        obs_style = plot.get_dataset_style('OBS', 'cox18nature')
        obs_cube = psi_cubes[obs_name]
        AXES.plot(obs_cube.coord('year').points,
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

        # Save plot
        provenance_record['plot_file'] = _save_fig(cfg, filename, legend)

    # Write provenance
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)


def plot_emergent_relationship(cfg, psi_cube, ecs_cube, lambda_cube, obs_cube):
    """Plot emergent relationship."""
    filename = 'emergent_relationship_{}'.format(
        obs_cube.attributes['dataset'])
    cube = ecs_cube.copy()
    cube.add_aux_coord(
        iris.coords.AuxCoord(psi_cube.data, **ih.convert_to_iris(PSI_ATTRS)),
        0)
    netcdf_path = get_diagnostic_filename(filename, cfg)
    io.iris_save(cube, netcdf_path)
    provenance_record = get_provenance_record(
        "Emergent relationship between ECS and the psi metric. The black dot-"
        "dashed line shows the best-fit linear regression across the model "
        "ensemble, with the prediction error for the fit given by the black "
        "dashed lines. The vertical blue lines show the observational "
        "constraint from the {} observations: the mean (dot-dashed line) and "
        "the mean plus and minus one standard deviation (dashed lines).".
        format(obs_cube.attributes['dataset']), ['mean', 'corr', 'var'],
        ['scatter'], _get_ancestor_files(cfg, obs_cube.attributes['dataset']))

    # Plot
    if cfg['write_plots']:
        obs_mean = np.mean(obs_cube.data)
        obs_std = np.std(obs_cube.data)

        # Calculate regression line
        lines = ec.regression_surface(psi_cube.data, ecs_cube.data,
                                      n_points=1000)
        logger.info("Found emergent relationship with slope %.2f (R2 = %.2f)",
                    lines['coef'], lines['R2'])

        # Plot points
        for model in psi_cube.coord('dataset').points:
            _plot_model_point(model, psi_cube, ecs_cube, lambda_cube)

        # Plot lines
        AXES.set_xlim(auto=False)
        AXES.set_ylim(auto=False)
        AXES.plot(lines['x'],
                  lines['y'],
                  color='black',
                  linestyle='dashdot',
                  label='Linear regression')
        AXES.plot(lines['x'],
                  lines['y_minus_err'],
                  color='black',
                  linestyle='dashed')
        AXES.plot(lines['x'],
                  lines['y_plus_err'],
                  color='black',
                  linestyle='dashed')
        AXES.axvline(obs_mean,
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
        provenance_record['plot_file'] = _save_fig(cfg, filename, legend)

    # Write provenance
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)


def plot_pdf(cfg, ecs_lin, ecs_pdf, ecs_cube, obs_name):
    """Plot probability density function of ECS."""
    filename = 'pdf_{}'.format(obs_name)
    netcdf_path = get_diagnostic_filename(filename, cfg)
    cube = iris.cube.Cube(ecs_pdf,
                          var_name='pdf',
                          long_name='Probability density function',
                          units='K-1')
    cube.add_aux_coord(
        iris.coords.AuxCoord(ecs_lin, **ih.convert_to_iris(ECS_ATTRS)), 0)
    io.iris_save(cube, netcdf_path)
    project = _get_project(cfg)
    provenance_record = get_provenance_record(
        "The PDF for ECS. The orange histograms show the prior distributions "
        "that arise from equal weighting of the {} models in 0.5 K bins.".
        format(project), ['mean'], ['other'],
        _get_ancestor_files(cfg, obs_name))

    # Plot
    if cfg['write_plots']:
        AXES.plot(ecs_lin,
                  ecs_pdf,
                  color='black',
                  linewidth=2.0,
                  label='Emergent constraint')
        AXES.hist(ecs_cube.data,
                  bins=6,
                  range=(2.0, 5.0),
                  density=True,
                  color='orange',
                  label='{} models'.format(project))

        # Plot appearance
        AXES.set_title('PDF of emergent constraint')
        AXES.set_xlabel('ECS / K')
        AXES.set_ylabel('Probability density')
        legend = AXES.legend(loc='upper left')

        # Save plot
        provenance_record['plot_file'] = _save_fig(cfg, filename, legend)

    # Write provenance
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)


def plot_cdf(cfg, ecs_lin, ecs_pdf, ecs_cube, obs_name):
    """Plot cumulative distribution function of ECS."""
    confidence_level = cfg.get('confidence_level', 0.66)
    ecs_cdf = ec.cdf(ecs_lin, ecs_pdf)

    # Provenance
    filename = 'cdf_{}'.format(obs_name)
    netcdf_path = get_diagnostic_filename(filename, cfg)
    cube = iris.cube.Cube(ecs_cdf,
                          var_name='cdf',
                          long_name='Cumulative distribution function',
                          units='1')
    cube.add_aux_coord(
        iris.coords.AuxCoord(ecs_lin, **ih.convert_to_iris(ECS_ATTRS)), 0)
    io.iris_save(cube, netcdf_path)
    project = _get_project(cfg)
    provenance_record = get_provenance_record(
        "The CDF for ECS. The horizontal dot-dashed lines show the {}% "
        "confidence limits. The orange histograms show the prior "
        "distributions that arise from equal weighting of the {} models in "
        "0.5 K bins.".format(int(confidence_level * 100), project), ['mean'],
        ['other'], _get_ancestor_files(cfg, obs_name))

    # Plot
    if cfg['write_plots']:
        AXES.plot(ecs_lin,
                  ecs_cdf,
                  color='black',
                  linewidth=2.0,
                  label='Emergent constraint')
        AXES.hist(ecs_cube.data,
                  bins=6,
                  range=(2.0, 5.0),
                  cumulative=True,
                  density=True,
                  color='orange',
                  label='{} models'.format(project))
        AXES.axhline((1.0 - confidence_level) / 2.0,
                     color='black',
                     linestyle='dashdot')
        AXES.axhline((1.0 + confidence_level) / 2.0,
                     color='black',
                     linestyle='dashdot')

        # Plot appearance
        AXES.set_title('CDF of emergent constraint')
        AXES.set_xlabel('ECS / K')
        AXES.set_ylabel('CDF')
        legend = AXES.legend(loc='upper left')

        # Save plot
        provenance_record['plot_file'] = _save_fig(cfg, filename, legend)

    # Write provenance
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)


def get_ecs_range(cfg, ecs_lin, ecs_pdf):
    """Get constrained ecs range."""
    confidence_level = cfg.get('confidence_level', 0.66)
    conf_low = (1.0 - confidence_level) / 2.0
    conf_high = (1.0 + confidence_level) / 2.0

    # Calculate CDF
    ecs_cdf = ec.cdf(ecs_lin, ecs_pdf)

    # Calculate constrained ECS range
    ecs_mean = ecs_lin[np.argmax(ecs_pdf)]
    ecs_index_range = np.where((ecs_cdf >= conf_low)
                               & (ecs_cdf <= conf_high))[0]
    ecs_range = ecs_lin[ecs_index_range]
    ecs_low = min(ecs_range)
    ecs_high = max(ecs_range)
    return (ecs_mean, ecs_low, ecs_high)


def main(cfg):
    """Run the diagnostic."""
    input_data = (
        select_metadata(cfg['input_data'].values(), short_name='tas') +
        select_metadata(cfg['input_data'].values(), short_name='tasa'))
    if not input_data:
        raise ValueError("This diagnostics needs 'tas' or 'tasa' variable")

    # Get time-dependent data
    (tas_cubes, tas_obs) = get_tas(input_data)
    (psi_cubes, psi_obs) = get_psi(cfg)

    # Get scalar psi, ECS and climate feedback parameter for models
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
        (ecs_lin, ecs_pdf) = ec.gaussian_pdf(psi_cube.data, ecs_cube.data,
                                             np.mean(obs_cube.data),
                                             np.var(obs_cube.data))
        plot_pdf(cfg, ecs_lin, ecs_pdf, ecs_cube, obs_name)
        plot_cdf(cfg, ecs_lin, ecs_pdf, ecs_cube, obs_name)

        # Print ECS range
        ecs_range = get_ecs_range(cfg, ecs_lin, ecs_pdf)
        logger.info("Observational constraint: Ψ = (%.2f ± %.2f) K",
                    np.mean(obs_cube.data), np.std(obs_cube.data))
        logger.info(
            "Constrained ECS range: (%.2f - %.2f) K with best "
            "estimate %.2f K", ecs_range[1], ecs_range[2], ecs_range[0])


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
        plt.close()
