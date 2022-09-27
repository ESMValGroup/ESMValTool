#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Plotting scatterplot

Description
-----------
This diagnostic creates a scatterplot.

Author
------
Lisa Bock (DLR, Germany)


Notes
-----
All configuration options starting with ``plot_`` specify keyword arguments for
a specific plot type. A certain plot type is only plotted if the corresponding
option is given in the recipe (if no additional keyword arguments are desired,
use ``{}``).

Configuration options in recipe
-------------------------------
additional_plot_kwargs_xy_plots: dict, optional
    Optional keyword arguments (values) for single datasets used in X-Y plots.
    They keys may include a ``var_type`` or values of the attribute given by
    ``group_by_attribute``.
alias: dict, optional
    :obj:`str` to :obj:`str` mapping for nicer plot labels (e.g.
    ``{'feature': 'Historical CMIP5 data'}``.
apply_common_mask: bool, optional (default: False)
    Apply common mask to all datasets prior to plotting. Requires identical
    shapes for all datasets.
group_attribute_as_default_alias: bool, optional (default: True)
    If ``True``, use value of attribute given by ``group_by_attribute`` as
    default alias if possible. If ``False``, use full group name (including
    ``var_type``) as default alias.
group_by_attribute: str, optional (default: 'mlr_model_name')
    By default, datasets are grouped using the ``var_type`` attribute. This
    option can be used to specify a further attribute to group datasets. This
    diagnostic expects a single dataset per group.
ignore: list of dict, optional
    Ignore specific datasets by specifying multiple :obj:`dict` s of metadata.
legend_kwargs: dict, optional
    Optional keyword arguments of :func:`matplotlib.pyplot.legend` (affects
    only plots with legends).
map_plot_type: str, optional (default: 'pcolormesh')
    Type of plot used for plotting maps. Must be one of ``'pcolormesh'`` or
    ``'contourf'``.
pattern: str, optional
    Pattern matched against ancestor file names.
plot_map: dict, optional
    Specify additional keyword arguments for plotting global maps showing
    datasets by ``plot_kwargs`` and plot appearance options by
    ``pyplot_kwargs`` (processed as functions of :mod:`matplotlib.pyplot`).
plot_map_abs_biases: dict, optional
    Specify additional keyword arguments for plotting global maps showing
    absolute biases by ``plot_kwargs`` and plot appearance options by
    ``pyplot_kwargs`` (processed as functions of :mod:`matplotlib.pyplot`).
plot_map_ratios: dict, optional
    Specify additional keyword arguments for plotting global maps showing
    ratios of datasets by ``plot_kwargs`` and plot appearance options by
    ``pyplot_kwargs`` (processed as functions of :mod:`matplotlib.pyplot`).
plot_map_rel_biases: dict, optional
    Specify additional keyword arguments for plotting global maps showing
    relative biases of datasets  by ``plot_kwargs`` and plot appearance options
    by ``pyplot_kwargs`` (processed as functions of :mod:`matplotlib.pyplot`).
plot_xy: dict, optional
    Specify additional keyword arguments for simple X-Y plots by
    ``plot_kwargs`` and plot appearance options by ``pyplot_kwargs`` (processed
    as functions of :mod:`matplotlib.pyplot`). By default, plots data against
    dimensional coordinate (if available). Use ``x_coord`` (:obj:`str`) to use
    another coordinate as X-axis. Use ``reg_line: True`` to additionally plot
    a linear regression line.
plot_xy_with_errors: dict, optional
    Specify additional keyword arguments for X-Y plots with error ranges
    ``plot_kwargs`` and plot appearance options by ``pyplot_kwargs`` (processed
    as functions of :mod:`matplotlib.pyplot`). By default, plots data against
    dimensional coordinate (if available). Use ``x_coord`` (:obj:`str`) to use
    another coordinate as X-axis.
print_corr: bool, optional (default: False)
    Print and save Pearson correlation coefficient between all datasets at the
    end.  Requires identical shapes for all datasets.
savefig_kwargs: dict, optional
    Keyword arguments for :func:`matplotlib.pyplot.savefig`.
seaborn_settings: dict, optional
    Options for :func:`seaborn.set` (affects all plots).
years_in_title: bool, optional (default: False)
    Print years in default title of plots.

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

ALL_CUBES = pd.DataFrame()
COLORS = sns.color_palette()
SEP = '___'


def _write_xy_provenance(cfg, cubes, plot_path, title, *attrs):
    """Write provenance information for X-Y plots."""
    cubes = cubes.copy()
    if isinstance(cubes, iris.cube.Cube):
        cubes = iris.cube.CubeList([cubes])
    ancestors = []
    for attr in attrs:
        ancestors.extend(attr['filename'].split('|'))
    netcdf_path = mlr.get_new_path(cfg, plot_path)
    io.iris_save(cubes, netcdf_path)
    long_name = ' and '.join([cube.long_name for cube in cubes])
    caption = f"Line plot of {long_name}"
    if title:
        caption += f" for {title}."
    else:
        caption += '.'
    record = {
        'ancestors': ancestors,
        'authors': ['schlund_manuel'],
        'caption': caption,
        'plot_types': ['line'],
        'references': ['schlund20jgr'],
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, record)
        provenance_logger.log(plot_path, record)


def _xy_plot(x_data, y_data, reg_line=False, **plot_kwargs):
    """Create single X-Y plot."""
    plot_kwargs = deepcopy(plot_kwargs)
    if reg_line:
        if plot_kwargs.get('linestyle', '-') == '-':
            plot_kwargs.setdefault('marker', 'o')
        else:
            plot_kwargs.setdefault('marker', 's')
        plot_kwargs['linestyle'] = 'none'
        plot_kwargs.setdefault('markersize', 3)
    #plt.plot(x_data, y_data, **plot_kwargs)
    if not reg_line:
        return
    plot_kwargs['linestyle'] = '-'
    plot_kwargs['marker'] = None
    plot_kwargs.pop('label', None)
    print(x_data)
    print(y_data)
    #reg = linregress(x_data, y_data)
    #y_reg = reg.slope * np.array(x_data) + reg.intercept
    #plt.plot(x_data, y_reg, **plot_kwargs)
    x, y = pd.Series(x_data, name="x_var"), pd.Series(y_data, name="y_var")
    sns.regplot(x=x_data, y=y_data, ci=95, truncate=False,
                line_kws={"color":"r","alpha":0.7,"lw":5})

    # Pearson correlation with SciPy:
    corr = pearsonr(x_data, y_data)
    corr = [np.round(c, 2) for c in corr]
    # Extracting the r-value and the p-value:
    title = 'r=%s, p=%s' % (corr[0], corr[1])
    # Adding the text to the Seaborn plot:
    plt.title(title, loc='right')


def get_plot_kwargs(cfg, option, key=None):
    """Get keyword arguments for desired plot function and key."""
    plot_kwargs = cfg.get(option, {}).get('plot_kwargs', {})
    print(plot_kwargs)
    if key is None:
        return plot_kwargs
    if '_xy' in option:
        additional_plot_kwargs = cfg.get('additional_plot_kwargs_xy_plots', {})
        if key in additional_plot_kwargs:
            return {**plot_kwargs, **additional_plot_kwargs[key]}
        subkey = key.split(SEP)[-1]
        if subkey in additional_plot_kwargs:
            return {**plot_kwargs, **additional_plot_kwargs[subkey]}
    return deepcopy(plot_kwargs)


def get_savefig_kwargs(cfg):
    """Get keyword arguments for :func:`matplotlib.pyplot.savefig`."""
    if 'savefig_kwargs' in cfg:
        return cfg['savefig_kwargs']
    savefig_kwargs = {
        'bbox_inches': 'tight',
        'dpi': 300,
        'orientation': 'landscape',
    }
    return savefig_kwargs


def process_pyplot_kwargs(cfg, option):
    """Process functions for :mod:`matplotlib.pyplot`."""
    for (key, val) in cfg.get(option, {}).get('pyplot_kwargs', {}).items():
        getattr(plt, key)(val)


def plot_scatter(x_data, y_data, cfg):
    """Plot scattterplot from dictionary."""

    plot_kwargs = get_plot_kwargs(cfg, 'plot_xy')

    _xy_plot(x_data, y_data, True, **plot_kwargs)

    #plt.title(title)
    plt.ylabel(cfg['ylabel'])
    plt.xlabel(cfg['xlabel'])
    plt.axvline(4.03, color='grey', linestyle='dashed')
    plt.axvline(2.87, color='grey', linestyle='dashed')

    plt.ylim(cfg.get('y_range'))

    process_pyplot_kwargs(cfg, 'plot_xy')
    plt.legend(**cfg.get('legend_kwargs'))
    plot_path = get_plot_filename(cfg['output_file_name'], cfg)
    savefig_kwargs = get_savefig_kwargs(cfg)
    plt.savefig(plot_path, **savefig_kwargs)
    logger.info("Wrote %s", plot_path)
    plt.close()
    #_write_xy_provenance(cfg, cubes, plot_path, None, *all_attrs)


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
      datasets = cube.coord('dataset').points
      for item in datasets:
        datasets1.append(item)
      for item in cube.data[:]:
        data1.append(item)
      logger.debug("Reading %s", input_file)

    input_files = io.get_all_ancestor_files(cfg, pattern=filename2)
    data2 = []
    datasets2 = []
    for input_file in input_files:
      logger.debug("Loading %s", input_file)
      cube = iris.load_cube(input_file)
      grid_areas = iris.analysis.cartography.area_weights(cube)
      new_cube = cube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN,
                                weights=grid_areas)
      data2.append(new_cube.data.item())
      datasets2.append(cube.attributes['dataset'])
      logger.debug("Reading %s", input_file)

    x_data = []
    y_data = []
    for idx, dataset in enumerate(datasets1):
        if dataset in datasets2:
            x_data.append(data1[idx])
            y_data.append(data2[datasets2.index(dataset)])
        else:
            logger.debug("Dataset %s is not part of file2", dataset)
    for idx, dataset in enumerate(datasets2):
        if dataset not in datasets1:
            logger.debug("Dataset %s is not part of file1", dataset)

    return x_data, y_data


def main(cfg):
    """Run the diagnostic."""
    sns.set(**cfg.get('seaborn_settings', {}))
    cfg = deepcopy(cfg)
    cfg.setdefault('legend_kwargs', {})
    cfg.setdefault('map_plot_type', 'pcolormesh')
    cfg.setdefault('print_corr', False)
    cfg.setdefault('years_in_title', False)
    cfg.setdefault('output_file_name', None)

    x_data, y_data = read_data_and_preprocess(cfg)

    plot_scatter(x_data, y_data, cfg)


# Run main function when this script is called
if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
