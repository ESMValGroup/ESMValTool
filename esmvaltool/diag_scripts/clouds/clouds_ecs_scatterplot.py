#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Plotting scatterplot

Description
-----------
This diagnostic creates a scatterplot.

Author
------
Lisa Bock (DLR, Germany)


Configuration options in recipe
-------------------------------
file1: contains the feedback parameters for all cmip models from the
       climate_metrics/feedback_parameters.py`` diagnostic
file2: contains the ECS values for all cmip models from the
       climate_metrics/ecs.py diagnostic
output_file_name: Set output file name.
xlabel: label on x-axis
ylabel: label on y-axis
x_range: Range for the x-axis of the plot
y_range: Range for the y-axis of the plot
print_corr: bool, optional (default: False)
    Print Pearson correlation coefficient between all datasets at the
    end.  Requires identical shapes for all datasets.
seaborn_settings: dict, optional
    Options for :func:`seaborn.set` (affects all plots).

"""

import itertools
import logging
import os
from copy import deepcopy
from pprint import pformat

import iris
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from cf_units import Unit
from scipy.stats import linregress
from scipy.stats import pearsonr

import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts import mlr
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    io,
    plot,
    run_diagnostic,
)

logger = logging.getLogger(os.path.basename(__file__))

def _write_xy_provenance(cfg, data_list, plot_path, title):
    """Write provenance information for X-Y plots."""
    cubes = cubes.copy()
    if isinstance(cubes, iris.cube.Cube):
        cubes = iris.cube.CubeList([cubes])
    #ancestors = []
    #    ancestors.extend(attr['filename'].split('|'))
    netcdf_path = mlr.get_new_path(cfg, plot_path)
    io.iris_save(cubes, netcdf_path)
    long_name = ' and '.join([cube.long_name for cube in cubes])
    caption = f"Scatter plot of cloud feedback and ECS"
    record = {
        'ancestors': ancestors,
        'authors': ['bock_lisa'],
        'caption': caption,
        'plot_types': ['scatter'],
        'references': ['bock24acp'],
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, record)
        provenance_logger.log(plot_path, record)


def _xy_plot(cfg, x_data, y_data, reg_line=False, **plot_kwargs):
    """Create single X-Y plot."""
    plot_kwargs = deepcopy(plot_kwargs)
    if reg_line:
        if plot_kwargs.get('linestyle', '-') == '-':
            plot_kwargs.setdefault('marker', 'o')
        else:
            plot_kwargs.setdefault('marker', 's')
        plot_kwargs['linestyle'] = 'none'
        plot_kwargs.setdefault('markersize', 3)
    if not reg_line:
        return
    plot_kwargs['linestyle'] = '-'
    plot_kwargs['marker'] = None
    plot_kwargs.pop('label', None)
    sns.regplot(x=x_data, y=y_data, ci=95, truncate=False,
                line_kws={"color":"r","alpha":0.7,"lw":5})

    if cfg['print_corr']:
        # Pearson correlation with SciPy:
        corr = pearsonr(x_data, y_data)
        print("Pearson correlation: r = ", corr[0], "p = ", corr[1])
        # Extracting the r-value and the p-value:
        corr = [np.round(c, 2) for c in corr]
        if corr[1] < 0.001:
            title = 'r=%s, p<<0.001' % (corr[0])
        else:
            title = 'r=%s, p=%s' % (corr[0], corr[1])
        # Adding the text to the Seaborn plot:
        plt.title(title, loc='right')


def plot_scatter(x_data, y_data, cfg):
    """Plot scattterplot from dictionary."""

    _xy_plot(cfg, x_data, y_data, True)

    plt.ylabel(cfg['ylabel'])
    plt.xlabel(cfg['xlabel'])
    plt.axhline(4.03, color='grey', linestyle='dashed')
    plt.axhline(2.87, color='grey', linestyle='dashed')

    plt.xlim(cfg.get('x_range'))
    plt.ylim(cfg.get('y_range'))

    plot_path = get_plot_filename(cfg['output_file_name'], cfg)
    plt.savefig(plot_path)
    logger.info("Wrote %s", plot_path)
    plt.close()
    #_write_xy_provenance(cfg, [x_data, y_data], plot_path, None)


def read_data_and_preprocess(cfg):
    """Read files and extract vars."""
    logger.debug("Reading data")
    filename1 = cfg['file1']
    filename2 = cfg['file2']

    input_files = io.get_all_ancestor_files(cfg, pattern=filename1)
    data1 = []
    datasets1 = []
    for input_file in input_files:
      logger.debug("Loading %s", input_file)
      cube = iris.load_cube(input_file)
      grid_areas = iris.analysis.cartography.area_weights(cube)

      cube.data = np.ma.masked_invalid(cube.data)
      new_cube = cube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN,
                                weights=grid_areas)
      data1.append(new_cube.data.item())
      datasets1.append(cube.attributes['dataset'])
      logger.debug("Reading %s", input_file)

    input_files = io.get_all_ancestor_files(cfg, pattern=filename2)
    data2 = []
    datasets2 = []
    for input_file in input_files:
      logger.debug("Loading %s", input_file)
      cube = iris.load_cube(input_file)
      datasets = cube.coord('dataset').points
      for item in datasets:
        datasets2.append(item)
      for item in cube.data[:]:
        data2.append(item)
      logger.debug("Reading %s", input_file)

    x_data = []
    y_data = []
    for idx, dataset in enumerate(datasets2):
        if dataset in datasets1:
            if not np.isnan(data1[datasets1.index(dataset)]):
                x_data.append(data1[datasets1.index(dataset)])
                y_data.append(data2[idx])
        else:
            logger.debug("Dataset %s is not part of file1", dataset)
    for idx, dataset in enumerate(datasets1):
        if dataset not in datasets2:
            logger.debug("Dataset %s is not part of file2", dataset)

    return x_data, y_data


def main(cfg):
    """Run the diagnostic."""
    sns.set(**cfg.get('seaborn_settings', {}))
    cfg = deepcopy(cfg)
    cfg.setdefault('print_corr', False)
    cfg.setdefault('output_file_name', None)

    x_data, y_data = read_data_and_preprocess(cfg)

    plot_scatter(x_data, y_data, cfg)


# Run main function when this script is called
if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
