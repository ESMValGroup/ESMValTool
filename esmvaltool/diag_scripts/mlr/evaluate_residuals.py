#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Simple evaluation of residuals (coming from MLR model output).

Description
-----------
This diagnostic evaluates residuals created by MLR models.

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recipe
-------------------------------
ignore: list of dict, optional
    Ignore specific datasets by specifying multiple :obj:`dict` s of metadata.
mse_plot: dict, optional
    Additional options for plotting the mean square errors (MSE). Specify
    additional keyword arguments for :func:`seaborn.boxplot` by ``plot_kwargs``
    and plot appearance options by ``pyplot_kwargs`` (processed as functions of
    :mod:`matplotlib.pyplot`).
pattern: str, optional
    Pattern matched against ancestor file names.
rmse_plot: dict, optional
    Additional options for plotting the root mean square errors (RMSE).
    Specify additional keyword arguments for :func:`seaborn.boxplot` by
    ``plot_kwargs`` and plot appearance options by ``pyplot_kwargs`` (processed
    as functions of :mod:`matplotlib.pyplot`).
savefig_kwargs: dict, optional
    Keyword arguments for :func:`matplotlib.pyplot.savefig`.
seaborn_settings: dict, optional
    Options for :func:`seaborn.set` (affects all plots).
weighted_samples: dict
    If specified, use weighted root mean square error. The given keyword
    arguments are directly passed to
    :func:`esmvaltool.diag_scripts.mlr.get_all_weights` to calculate the sample
    weights. By default, area weights and time weights are used.

"""

import logging
import os
from copy import deepcopy

import iris
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import esmvaltool.diag_scripts.emergent_constraints as ec
import esmvaltool.diag_scripts.mlr.plot as mlr_plot
from esmvaltool.diag_scripts import mlr
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_plot_filename,
    group_metadata,
    io,
    run_diagnostic,
    select_metadata,
)

logger = logging.getLogger(os.path.basename(__file__))


def _plot_boxplot(cfg, data_frame, plot_name):
    """Plot boxplot."""
    boxplot_kwargs = {
        'color': 'b',
        'data': data_frame,
        'showfliers': False,
        'showmeans': True,
        'meanprops': {
            'marker': 'x',
            'markeredgecolor': 'k',
            'markerfacecolor': 'k',
            'markersize': 8,
        },
        'whis': [0, 100],
    }
    boxplot_kwargs.update(mlr_plot.get_plot_kwargs(cfg, plot_name))
    sns.boxplot(**boxplot_kwargs)
    sns.swarmplot(data=data_frame, color='k', alpha=0.6)

    # Plot appearance
    plt.ylim(0.0, plt.ylim()[1])
    mlr_plot.process_pyplot_kwargs(cfg, plot_name)

    # Save plot
    plot_path = get_plot_filename(plot_name, cfg)
    plt.savefig(plot_path, **mlr_plot.get_savefig_kwargs(cfg))
    logger.info("Wrote %s", plot_path)
    plt.close()
    return plot_path


def _write_provenance(cfg, data_frame, plot_path, title, ancestors,
                      **cube_kwargs):
    """Write provenance information."""
    cube = ec.pandas_object_to_cube(data_frame, **cube_kwargs)
    netcdf_path = mlr.get_new_path(cfg, plot_path)
    io.iris_save(cube, netcdf_path)
    record = {
        'ancestors': ancestors,
        'authors': ['schlund_manuel'],
        'caption': f"Boxplot of {title}.",
        'plot_file': plot_path,
        'plot_types': ['box'],
        'references': ['schlund20jgr'],
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, record)


def get_residual_data(cfg):
    """Get residual data."""
    input_data = mlr_plot.get_input_datasets(cfg)
    residual_data = select_metadata(input_data, var_type='prediction_residual')
    if not residual_data:
        raise ValueError("No 'prediction_residual' data found")
    return group_metadata(residual_data, 'mlr_model_name')


def plot_mse(cfg, residual_data):
    """Plot distribution of mean square error (MSE)."""
    logger.info("Plotting mean square error (MSE) distribution")
    mlr_models_mse = []

    # Collect data for every statistical model
    ancestors = []
    for (model_name, datasets) in residual_data.items():
        mse_data = []
        for dataset in datasets:
            cube = iris.load_cube(dataset['filename'])
            ancestors.append(dataset['filename'])
            weights = mlr.get_all_weights(cube, **cfg['weighted_samples'])
            mse = np.ma.average(cube.data**2, weights=weights)
            mse_data.append(mse)
        data_frame = pd.DataFrame(mse_data, columns=[model_name])
        mlr_models_mse.append(data_frame)
    boxplot_data_frame = pd.concat(mlr_models_mse, axis=1)
    boxplot_data_frame.columns.name = 'mlr_model'

    # Plot
    plot_path = _plot_boxplot(cfg, boxplot_data_frame, 'mse_plot')
    logger.info("MSEs:\n%s", boxplot_data_frame.describe())

    # Provenance
    _write_provenance(cfg, boxplot_data_frame, plot_path, 'MSE', ancestors,
                      var_name='mse', long_name='Mean Square Error')


def plot_rmse(cfg, residual_data):
    """Plot distribution of root mean square error (RMSE)."""
    logger.info("Plotting root mean square error (RMSE) distribution")
    mlr_models_rmse = []

    # Collect data for every statistical model
    ancestors = []
    for (model_name, datasets) in residual_data.items():
        rmse_data = []
        for dataset in datasets:
            cube = iris.load_cube(dataset['filename'])
            ancestors.append(dataset['filename'])
            weights = mlr.get_all_weights(cube, **cfg['weighted_samples'])
            mse = np.ma.average(cube.data**2, weights=weights)
            rmse_data.append(np.ma.sqrt(mse))
        data_frame = pd.DataFrame(rmse_data, columns=[model_name])
        mlr_models_rmse.append(data_frame)
    boxplot_data_frame = pd.concat(mlr_models_rmse, axis=1)
    boxplot_data_frame.columns.name = 'mlr_model'

    # Plot
    plot_path = _plot_boxplot(cfg, boxplot_data_frame, 'rmse_plot')
    logger.info("RMSEs:\n%s", boxplot_data_frame.describe())

    # Provenance
    _write_provenance(cfg, boxplot_data_frame, plot_path, 'RMSE', ancestors,
                      var_name='rmse', long_name='Root Mean Square Error')


def main(cfg):
    """Run the diagnostic."""
    cfg = deepcopy(cfg)
    cfg.setdefault('weighted_samples',
                   {'area_weighted': True, 'time_weighted': True})
    sns.set(**cfg.get('seaborn_settings', {}))

    # Extract data
    residual_data = get_residual_data(cfg)

    # Plots
    plot_mse(cfg, residual_data)
    plot_rmse(cfg, residual_data)


# Run main function when this script is called
if __name__ == '__main__':
    mlr.ignore_warnings()
    with run_diagnostic() as config:
        main(config)
