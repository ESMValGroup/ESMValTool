#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to plot figure 9.42b of IPCC AR5 chapter 9.

Description
-----------
Calculate and plot the effective climate sensitivity (ECS) vs. the transient
climate response (TCR) (see IPCC AR5 WG1 ch.9, fig. 9.42b).

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
    :mod:`esmvaltool.diag_scripts.shared.plot.styles_python`).
seaborn_settings : dict, optional
    Options for seaborn's `set()` method (affects all plots), see
    <https://seaborn.pydata.org/generated/seaborn.set.html>.

"""

import logging
import os

import iris
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats

from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            get_plot_filename, io,
                                            iris_helpers, plot, run_diagnostic)

logger = logging.getLogger(os.path.basename(__file__))


def _get_reg_line(x_cube, y_cube, xlim=(1.5, 5.0), n_points=100):
    """Get regression line and means for two cubes."""
    reg = stats.linregress(x_cube.data, y_cube.data)
    x_reg = np.linspace(xlim[0], xlim[1], n_points)
    y_reg = reg.slope * x_reg + reg.intercept
    x_mean = np.mean(x_cube.data)
    y_mean = np.mean(y_cube.data)
    return ((x_reg, y_reg), (x_mean, y_mean))


def get_provenance_record(project, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption':
        (f'Effective climate sensitivity (ECS) against transient climate '
         f'response of {project} models.'),
        'statistics': ['mean'],
        'domains': ['global'],
        'authors': ['schl_ma'],
        'references': ['flato13ipcc'],
        'realms': ['atmos'],
        'themes': ['phys'],
        'ancestors':
        ancestor_files,
    }
    return record


def plot_data(cfg, ecs_cube, tcr_cube):
    """Plot data."""
    if not cfg['write_plots']:
        return None
    logger.debug("Plotting Fig. 9.42b of IPCC AR5")
    (_, axes) = plt.subplots()
    project = ecs_cube.attributes['project']

    # Plot scatterplot
    for dataset_name in ecs_cube.coord('dataset').points:
        style = plot.get_dataset_style(dataset_name, cfg.get('dataset_style'))
        ecs = ecs_cube.extract(iris.Constraint(dataset=dataset_name)).data
        tcr = tcr_cube.extract(iris.Constraint(dataset=dataset_name)).data
        axes.plot(ecs,
                  tcr,
                  marker=style['mark'],
                  linestyle='none',
                  markeredgecolor=style['color'],
                  markerfacecolor=style['facecolor'],
                  label=dataset_name)

    # Plot regression line and MMM
    (reg_line, mmm) = _get_reg_line(ecs_cube, tcr_cube)
    axes.plot(reg_line[0], reg_line[1], 'k-')
    axes.plot(mmm[0], mmm[1], 'ro', label=f'{project} mean')

    # Plot appearance
    axes.set_title(f"TCR vs. ECS for {project} models")
    axes.set_xlabel("ECS / K")
    axes.set_ylabel("TCR / K")
    axes.set_xlim([1.5, 5.0])
    axes.set_ylim([0.5, 3.0])
    legend = axes.legend(loc='center left',
                         bbox_to_anchor=[1.05, 0.5],
                         borderaxespad=0.0,
                         ncol=2)

    # Save plot
    plot_path = get_plot_filename('ch09_fig09_42b', cfg)
    plt.savefig(plot_path,
                orientation='landscape',
                bbox_inches='tight',
                additional_artists=[legend])
    logger.info("Wrote %s", plot_path)
    plt.close()
    return plot_path


def write_data(cfg, ecs_cube, tcr_cube):
    """Write netcdf file."""
    ecs_attrs = {
        'var_name': ecs_cube.var_name,
        'long_name': ecs_cube.long_name,
        'units': ecs_cube.units,
    }

    # Write data
    ecs_coord = iris.coords.AuxCoord(ecs_cube.data, **ecs_attrs)
    tcr_cube.add_aux_coord(ecs_coord, 0)
    tcr_cube.attributes.pop('provenance', None)
    netcdf_path = get_diagnostic_filename('ch09_fig09_42b', cfg)
    io.iris_save(tcr_cube, netcdf_path)
    return netcdf_path


def main(cfg):
    """Run the diagnostic."""
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
    if plot_path is not None:
        provenance_record.update({
            'plot_file': plot_path,
            'plot_types': ['scatter'],
        })
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
