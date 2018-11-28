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
import pprint

import iris
from iris import Constraint
import numpy as np

from esmvaltool.diag_scripts.shared import (extract_variables, group_metadata,
                                            get_all_ancestor_files, plot,
                                            run_diagnostic, save_iris_cube,
                                            variables_available,
                                            get_file_from_ancestors)
from esmvaltool.diag_scripts.shared.emergent_constraints import (
    standard_prediction_error)

import matplotlib # noqa
matplotlib.use('Agg') # noqa
import matplotlib.pyplot as plt # noqa

logger = logging.getLogger(os.path.basename(__file__))
plt.style.use(plot.get_path_to_mpl_style())


def constraint_no_obs(cfg):
    """Create `iris.Constraint` to remove OBS data from cubes."""
    datasets = []
    grouped_data = group_metadata(cfg['input_data'].values(), 'project')
    for data in grouped_data.get('OBS', []):
        datasets.append(data['dataset'])

    # Constraint function
    def no_obs(cell):
        return cell not in datasets

    return iris.Constraint(dataset=no_obs)


def calculate_emergent_relationsip(cfg, psi_cube, ecs_cube):
    """Calculate emergent relationship."""
    lines = {}
    psi_cube_models = psi_cube.extract(constraint_no_obs(cfg))
    ecs_cube_models = ecs_cube.extract(constraint_no_obs(cfg))
    if (set(psi_cube_models.coord('dataset').points) !=
            set(ecs_cube_models.coord('dataset').points)):
        raise ValueError("Expected identical climate models for psi and ECS")

    # Statistical calculations
    spe = standard_prediction_error(psi_cube_models.data, ecs_cube_models.data)






def plot_emergent_relationship(cfg, psi_cube, ecs_cube, lambda_cube, lines):
    """Plot emergent relationship."""
    if not cfg['write_plots']:
        return
    (fig, axes) = plt.subplots()


def main(cfg):
    """Run the diagnostic."""
    input_data = cfg['input_data'].values()

    # Check if tas is available
    if not variables_available(cfg, ['tas']):
        raise ValueError("This diagnostic needs 'tas' variable")

    # Get external data (ECS, lambda and psi)
    ecs_filepath = get_file_from_ancestors(cfg, 'ecs.nc')
    lambda_filepath = get_file_from_ancestors(cfg, 'lambda.nc')
    psi_filepath = get_file_from_ancestors(cfg, 'psi.nc')
    psi_files = get_all_ancestor_files(cfg, pattern='psi_*.nc')
    ecs_cube = iris.load_cube(ecs_filepath)
    lambda_cube = iris.load_cube(lambda_filepath)
    psi_cube = iris.load_cube(psi_filepath)

    # Calculation
    calculate_emergent_relationsip(cfg, psi_cube, ecs_cube)

    # Plots
    plot_emergent_relationship(cfg, psi_cube, ecs_cube, lambda_cube)

    # # Create iris cubes for each dataset
    # hist_cubes = {}
    # pi_cubes = {}
    # for data in input_data:
    #     name = data['dataset']
    #     logger.info("Processing %s", name)
    #     cube = iris.load_cube(data['filename'])

    #     # Preprocess cubes
    #     cube.convert_units(cfg.get('tas_units', 'celsius'))
    #     cube = cube.collapsed(['time'], iris.analysis.MEAN)

    #     # Save cubes
    #     if data.get('exp') == 'historical':
    #         hist_cubes[name] = cube
    #     elif data.get('exp') == 'piControl':
    #         pi_cubes[name] = cube
    #     else:
    #         pass

    # # Plot data
    # plot_data(cfg, hist_cubes, pi_cubes, ecs_cube)

    # # Write netcdf file
    # write_data(cfg, hist_cubes, pi_cubes, ecs_cube)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
