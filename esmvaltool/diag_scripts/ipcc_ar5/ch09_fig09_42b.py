#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to plot figure 9.42b of IPCC AR5 chapter 9.

Description
-----------
Calculate and plot the transient climate response (TCR) vs. the equilibrium
climate sensitivity (ECS) (see IPCC AR5 WG1 ch.9, fig. 9.42b).

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recipe
-------------------------------
dataset_style : str, optional
    Dataset style file (located in
    :mod:`esmvaltool.diag_scripts.shared.plot.styles_python`). The entry
    ``marker`` is ignored when ``marker_file`` is given.
log_x : bool, optional (default: False)
    Apply logarithm to X axis (ECS).
log_y : bool, optional (default: False)
    Apply logarithm to Y axis (TCR).
marker_column : str, optional (default: 'marker')
    Name of the column to look up markers in ``marker_file``.
marker_file : str, optional
    CSV file with markers (can also be integers). Must have the columns
    ``dataset`` and ``marker`` (or the column specified by ``marker_column``).
    If a relative path is given, assumes that this is a pattern to search for
    ancestor files.
savefig_kwargs : dict, optional
    Keyword arguments for :func:`matplotlib.pyplot.savefig`.
seaborn_settings : dict, optional
    Options for :func:`seaborn.set` (affects all plots).
x_lim : list of float, optional (default: [1.5, 6.0])
    Plot limits for X axis (ECS).
y_lim : list of float, optional (default: [0.5, 3.5])
    Plot limits for Y axis (TCR).

"""

import logging
import os
from copy import deepcopy

import iris
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    io,
    iris_helpers,
    plot,
    run_diagnostic,
)

logger = logging.getLogger(os.path.basename(__file__))

COLORS = sns.color_palette()


def _get_reg_line(x_cube, y_cube, cfg, n_points=100):
    """Get regression line and means for two cubes."""
    x_lim = [cfg['x_lim'][0] - 1.0, cfg['x_lim'][1] + 1.0]
    x_mean = np.mean(x_cube.data)
    y_mean = np.mean(y_cube.data)

    # Apply logarithms if desired
    if cfg.get('log_x'):
        x_cube = x_cube.copy(np.ma.log(x_cube.data))
    if cfg.get('log_y'):
        y_cube = y_cube.copy(np.ma.log(y_cube.data))

    # Regression
    reg = stats.linregress(x_cube.data, y_cube.data)
    logger.info("Regression stats")
    logger.info("Slope: %.2f", reg.slope)
    logger.info("Intercept: %.2f", reg.intercept)
    logger.info("R2: %.2f", reg.rvalue**2)

    # Regression line
    x_reg = np.linspace(x_lim[0], x_lim[1], n_points)
    y_reg = reg.slope * x_reg + reg.intercept
    if cfg.get('log_x'):
        x_reg = np.exp(x_reg)
    if cfg.get('log_y'):
        y_reg = np.exp(y_reg)

    return ((x_reg, y_reg), (x_mean, y_mean), reg.rvalue)


def _get_style(dataset_name, cfg):
    """Get style for individual data points."""
    style = plot.get_dataset_style(dataset_name, cfg.get('dataset_style'))
    if not cfg.get('marker_file'):
        return style
    marker_file = os.path.expanduser(cfg['marker_file'])
    if not os.path.isabs(marker_file):
        marker_file = io.get_ancestor_file(cfg, marker_file)
    data_frame = pd.read_csv(marker_file)
    marker_column = cfg['marker_column']
    for column in ('dataset', marker_column):
        if column not in data_frame.columns:
            raise ValueError(
                f"Marker file '{marker_file}' does not contain necessary "
                f"column '{column}'")
    marker = data_frame[marker_column][data_frame['dataset'] == dataset_name]
    if len(marker) != 1:
        raise ValueError(
            f"Expected exactly one entry for marker of '{dataset_name}' in "
            f"file '{marker_file}', got {len(marker):d}")
    style['mark'] = marker.values[0]
    return style


def get_provenance_record(project, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption':
        (f'Transient climate response (TCR) against equilibrium climate '
         f'sensitivity (ECS) for {project} models.'),
        'statistics': ['mean'],
        'domains': ['global'],
        'authors': ['schlund_manuel'],
        'references': ['flato13ipcc'],
        'realms': ['atmos'],
        'themes': ['phys'],
        'ancestors':
        ancestor_files,
    }
    return record


def plot_data(cfg, ecs_cube, tcr_cube):
    """Plot data."""
    logger.debug("Plotting Fig. 9.42b of IPCC AR5")
    (_, axes) = plt.subplots()
    project = ecs_cube.attributes['project']

    # Plot scatterplot
    plot_types = []
    for dataset_name in ecs_cube.coord('dataset').points:
        if dataset_name == 'MultiModelMean':
            continue
        style = _get_style(dataset_name, cfg)

        # Plot single point
        if isinstance(style['mark'], str):
            axes.plot(
                ecs_cube.extract(iris.Constraint(dataset=dataset_name)).data,
                tcr_cube.extract(iris.Constraint(dataset=dataset_name)).data,
                marker=style['mark'],
                linestyle='none',
                markeredgecolor=style['color'],
                markerfacecolor=style['facecolor'],
                label=dataset_name,
            )
            plot_types.append(0)
        else:
            axes.text(
                ecs_cube.extract(iris.Constraint(dataset=dataset_name)).data,
                tcr_cube.extract(iris.Constraint(dataset=dataset_name)).data,
                str(int(style['mark'])),
                size=7,
                ha='center',
                va='center',
            )
            plot_types.append(1)
    plot_type = np.mean(plot_types)

    # Plot regression line and MMM
    (reg_line, mmm, r_value) = _get_reg_line(ecs_cube, tcr_cube, cfg)
    if plot_type <= 0.5:
        axes.plot(reg_line[0], reg_line[1], 'k-')
        axes.plot(mmm[0], mmm[1], 'ro', label=f'{project} mean')
    else:
        axes.plot(reg_line[0], reg_line[1], linestyle='-', color=COLORS[0])

    # Plot appearance
    title = f"Linear TCR vs. ECS for {project} models"
    axes.set_xlim(cfg['x_lim'])
    axes.set_ylim(cfg['y_lim'])
    if cfg.get('log_x'):
        axes.set_xscale('log')
        axes.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
        axes.get_xaxis().set_minor_formatter(ticker.ScalarFormatter())
        axes.grid(True, which='both', axis='x')
        title = f"Non-linear TCR vs. ECS for {project} models"
    if cfg.get('log_y'):
        axes.set_yscale('log')
        axes.get_yaxis().set_major_formatter(ticker.ScalarFormatter())
        axes.get_yaxis().set_minor_formatter(ticker.ScalarFormatter())
        axes.grid(True, which='both', axis='y')
        title = f"Non-linear TCR vs. ECS for {project} models"
    axes.set_title(title)
    axes.set_xlabel('ECS [K]')
    axes.set_ylabel('TCR [K]')
    if plot_type <= 0.5:
        legend = axes.legend(loc='center left',
                             bbox_to_anchor=[1.05, 0.5],
                             borderaxespad=0.0,
                             ncol=2)
    else:
        legend = None
    axes.text(0.05, 0.9, f'$R^2$ = {r_value**2:.2f}', transform=axes.transAxes)

    # Save plot
    plot_path = get_plot_filename(f'{project}_ch09_fig09_42b', cfg)
    plt.savefig(plot_path,
                additional_artists=[legend],
                **cfg['savefig_kwargs'])
    logger.info("Wrote %s", plot_path)
    plt.close()
    return plot_path


def set_default_cfg(cfg):
    """Set default values for cfg."""
    cfg = deepcopy(cfg)
    cfg.setdefault('marker_column', 'marker')
    cfg.setdefault('savefig_kwargs', {
        'dpi': 300,
        'orientation': 'landscape',
        'bbox_inches': 'tight',
    })
    cfg.setdefault('x_lim', [1.5, 6.0])
    cfg.setdefault('y_lim', [0.5, 3.5])
    cfg['x_lim'] = np.array(cfg['x_lim'])
    cfg['y_lim'] = np.array(cfg['y_lim'])
    return cfg


def write_data(cfg, ecs_cube, tcr_cube):
    """Write netcdf file."""
    project = ecs_cube.attributes['project']
    ecs_attrs = {
        'var_name': ecs_cube.var_name,
        'long_name': ecs_cube.long_name,
        'units': ecs_cube.units,
    }

    # Write data to netcdf
    ecs_coord = iris.coords.AuxCoord(ecs_cube.data, **ecs_attrs)
    tcr_cube.add_aux_coord(ecs_coord, 0)
    tcr_cube.attributes.pop('provenance', None)
    netcdf_path = get_diagnostic_filename(f'{project}_ch09_fig09_42b', cfg)
    io.iris_save(tcr_cube, netcdf_path)
    return netcdf_path


def main(cfg):
    """Run the diagnostic."""
    cfg = set_default_cfg(cfg)
    sns.set(**cfg.get('seaborn_settings', {}))
    ecs_file = io.get_ancestor_file(cfg, 'ecs.nc')
    tcr_file = io.get_ancestor_file(cfg, 'tcr.nc')
    ecs_cube = iris.load_cube(ecs_file)
    tcr_cube = iris.load_cube(tcr_file)

    # Project
    if (ecs_cube.attributes.get('project', 'a') != tcr_cube.attributes.get(
            'project', 'b')):
        raise ValueError(
            "ECS and TCR input files have either no 'project' attribute or "
            "differ in it")
    project = ecs_cube.attributes['project']

    # Remove missing data and use equal coordinate
    [ecs_cube, tcr_cube
     ] = iris_helpers.intersect_dataset_coordinates([ecs_cube, tcr_cube])

    # Create plot
    plot_path = plot_data(cfg, ecs_cube, tcr_cube)

    # Write netcdf file
    netcdf_path = write_data(cfg, ecs_cube, tcr_cube)

    # Provenance
    ancestor_files = [ecs_file, tcr_file]
    provenance_record = get_provenance_record(project, ancestor_files)
    provenance_record.update({
        'plot_types': ['scatter'],
    })
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)
        provenance_logger.log(plot_path, provenance_record)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
