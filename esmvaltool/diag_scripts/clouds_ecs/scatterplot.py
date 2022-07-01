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


def _add_correlation_information(cfg, title, cube):
    """Add data from cube to :class:`pandas.DataFrame` holding all data."""
    if not cfg['print_corr']:
        return
    if not ALL_CUBES.empty and len(ALL_CUBES.index) != cube.data.size:
        raise ValueError(
            "Expected datasets with identical shapes when 'print_corr' is set")
    ALL_CUBES[title] = np.ma.filled(cube.data.ravel(), np.nan)


def _get_alias(cfg, name):
    """Get alias for given ``name``."""
    aliases = cfg.get('aliases', {})
    if name in aliases:
        return aliases[name]
    if cfg['group_attribute_as_default_alias']:
        return name.split(SEP)[-1]
    return name


def _get_cube(var_type, group_by_attribute, attr, datasets):
    """Get single cube for datasets of type ``key``."""
    key = _get_key(var_type, attr)
    logger.info("Found the following datasets for '%s':\n%s", key,
                pformat([d['filename'] for d in datasets]))
    if 'error' in var_type:
        logger.debug("Calculating cube for '%s' by squared error aggregation",
                     key)
        ref_cube = iris.load_cube(datasets[0]['filename'])
        cube = mlr.get_squared_error_cube(ref_cube, datasets)
        mlr.square_root_metadata(cube)
        cube.data = np.ma.sqrt(cube.data)
    else:
        if len(datasets) != 1:
            raise ValueError(f"Expected exactly one dataset for '{key}' got "
                             f"{len(datasets):d}:\n"
                             f"{pformat([d['filename'] for d in datasets])}")
        cube = iris.load_cube(datasets[0]['filename'])
    dataset_names = sorted(list({d['dataset'] for d in datasets}))
    end_years = list({d['end_year'] for d in datasets})
    filenames = sorted(list({d['filename'] for d in datasets}))
    projects = sorted(list({d['project'] for d in datasets}))
    start_years = list({d['start_year'] for d in datasets})
    cube.attributes.update({
        'dataset': '|'.join(dataset_names),
        'end_year': max(end_years),
        'filename': '|'.join(filenames),
        'project': '|'.join(projects),
        'start_year': min(start_years),
        'tag': datasets[0]['tag'],
        'var_type': var_type,
    })
    if attr is not None:
        cube.attributes[group_by_attribute] = attr
    _unify_time_coord(cube)
    return cube


def _get_key(var_type, attr):
    """Get dictionary key for specific dataset."""
    if attr is None:
        return var_type
    return f'{var_type}{SEP}{attr}'


def _get_levels(val_min, val_max):
    """Get symmetric levels for contour plot.

    This function might be changed to consider cube in the future.

    """
    n_levels_per_sign = 50
    #val_max = 8.0
    #val_min = -8.0
    max_range = max([val_max, -val_min])
    range_ = np.linspace(0.0, max_range, n_levels_per_sign + 1)[1:]
    levels = list(-range_[::-1]) + [0.0] + list(range_)
    return levels


def _get_map_plot_func(cfg):
    """Get function used for plotting maps."""
    allowed_funcs = {
        'contourf': plot.global_contourf,
        'pcolormesh': plot.global_pcolormesh,
    }
    if cfg['map_plot_type'] not in allowed_funcs:
        raise ValueError(
            f"Expected one of {list(allowed_funcs.keys())} for "
            f"'map_plot_type', got '{cfg['map_plot_type']}'")
    return allowed_funcs[cfg['map_plot_type']]


def _get_title(cfg, alias_1, attrs_1, alias_2=None, attrs_2=None,
               op_type='-'):
    """Get title for plots."""
    if alias_2 is None:
        title = alias_1
        if cfg['years_in_title']:
            title += f" ({attrs_1['start_year']}-{attrs_1['end_year']})"
        return title
    if attrs_2 is None:
        raise ValueError(
            "'attrs_2' needs to be given when 'alias_2' is not None")
    if op_type == 'rel_bias':
        if not cfg['years_in_title']:
            title = f"({alias_1} - {alias_2}) / {alias_2}"
            return title
        if (attrs_1['start_year'] == attrs_2['start_year']
                and attrs_1['end_year'] == attrs_2['end_year']):
            title = (f"({alias_1} - {alias_2}) / {alias_2} "
                     f"({attrs_1['start_year']}-{attrs_1['end_year']})")
        else:
            title = (f"({alias_1} ({attrs_1['start_year']}-"
                     f"{attrs_1['end_year']}) - {alias_2} ("
                     f"{attrs_2['start_year']}-{attrs_2['end_year']})) / "
                     f"{alias_2} ({attrs_2['start_year']}-"
                     f"{attrs_2['end_year']})")
        return title
    if not cfg['years_in_title']:
        title = f"{alias_1} {op_type} {alias_2}"
        return title
    if (attrs_1['start_year'] == attrs_2['start_year']
            and attrs_1['end_year'] == attrs_2['end_year']):
        title = (f"{alias_1} {op_type} {alias_2} ({attrs_1['start_year']}-"
                 f"{attrs_1['end_year']})")
    else:
        title = (f"{alias_1} ({attrs_1['start_year']}-{attrs_1['end_year']}) "
                 f"{op_type} {alias_2} ({attrs_2['start_year']}-"
                 f"{attrs_2['end_year']})")
    return title


def _mask_cube(cube):
    """Mask cube to avoid divisions by zero."""
    cube = cube.copy()
    val_range = np.ma.max(cube.data) - np.ma.min(cube.data)
    threshold = val_range * 5e-2
    cube.data = np.ma.masked_inside(cube.data, -threshold, threshold)
    return cube


def _unify_time_coord(cube):
    """Unify time coordinate of cube."""
    if not cube.coords('time', dim_coords=True):
        return
    time_coord = cube.coord('time')
    dates_points = time_coord.units.num2date(time_coord.points)
    dates_bounds = time_coord.units.num2date(time_coord.bounds)
    new_units = Unit('days since 1850-01-01 00:00:00')
    new_time_coord = iris.coords.DimCoord(
        new_units.date2num(dates_points),
        bounds=new_units.date2num(dates_bounds),
        var_name='time',
        standard_name='time',
        long_name='time',
        units=new_units,
    )
    coord_dims = cube.coord_dims('time')
    cube.remove_coord('time')
    cube.add_dim_coord(new_time_coord, coord_dims)


def _write_map_provenance(cfg, cube, plot_path, title, *attrs):
    """Write provenance information for map plots."""
    cube = cube.copy()
    ancestors = []
    for attr in attrs:
        ancestors.extend(attr['filename'].split('|'))
    netcdf_path = mlr.get_new_path(cfg, plot_path)
    io.iris_save(cube, netcdf_path)
    record = {
        'ancestors': ancestors,
        'authors': ['schlund_manuel'],
        'caption': f"Geographical distribution of {cube.long_name} for "
                   f"{title}.",
        'plot_types': ['geo'],
        'references': ['schlund20jgr'],
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, record)
        provenance_logger.log(plot_path, record)


def _write_xy_error_provenance(cfg, cubes, plot_path, title, ancestors):
    """Write provenance information for X-Y plots with error range."""
    cubes = cubes.copy()
    if isinstance(cubes, iris.cube.Cube):
        cubes = iris.cube.CubeList([cubes])
    netcdf_path = mlr.get_new_path(cfg, plot_path)
    io.iris_save(cubes, netcdf_path)
    long_name = ' and '.join([cube.long_name for cube in cubes])
    caption = f"Line plot with error bars of {long_name}"
    if title:
        caption += f" for {title}."
    else:
        caption += '.'
    record = {
        'ancestors': ancestors,
        'authors': ['schlund_manuel'],
        'caption': caption,
        'plot_types': ['line', 'errorbar'],
        'references': ['schlund20jgr'],
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, record)
        provenance_logger.log(plot_path, record)


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


def _xy_plot_with_errors(cfg, cube_dict, split_key, **plot_kwargs):
    """Create single X-Y plot with error ranges."""
    ancestors = []
    plot_kwargs = deepcopy(plot_kwargs)
    key = SEP.join(split_key)
    error_key = split_key[0] + '_error'
    if len(split_key) > 1:
        error_key = SEP.join([error_key, *split_key[1:]])
    if error_key not in cube_dict:
        raise ValueError(
            f"Corresponding error '{error_key}' for '{key}' not available")
    x_coord = cfg['plot_xy_with_errors'].get('x_coord')

    # Extract data
    cube = cube_dict[key]
    error_cube = cube_dict[error_key]
    ancestors.extend(cube.attributes['filename'].split('|'))
    ancestors.extend(error_cube.attributes['filename'].split('|'))
    if cube.ndim != 1 or error_cube.ndim != 1:
        raise ValueError(
            f"Expected 1D cube for X-Y plots with error range, got "
            f"{cube.ndim:d}D and {error_cube.ndim:d} (error) cubes")
    if x_coord is not None:
        coord = cube.coord(x_coord)
        x_data = coord.points
    elif cube.coords(dim_coords=True):
        coord = cube.coord(dim_coords=True)
        x_data = coord.points
    else:
        coord = None
        x_data = np.arange(cube.shape[0])
    if coord is not None:
        if not error_cube.coord(coord):
            raise iris.exceptions.CoordinateNotFoundError(
                f"Coordinate '{coord.name()}' of '{key}' not found in "
                f"corresponding error '{error_key}'")

    # Plot
    alias = _get_alias(cfg, key)
    plot_kwargs.setdefault('color', COLORS[0])
    plot_kwargs.setdefault('label', alias)
    plt.plot(x_data, cube.data, **plot_kwargs)
    plot_kwargs.pop('label')
    plot_kwargs['alpha'] = 0.12
    plt.fill_between(x_data, cube.data - error_cube.data,
                     cube.data + error_cube.data, **plot_kwargs)

    return (cube, error_cube, coord, ancestors)


def get_cube_dict(cfg, group_by_attribute):
    """Get dictionary of mean cubes (values) with ``var_type`` (keys)."""
    logger.info("Grouping datasets by 'var_type' and '%s'", group_by_attribute)
    input_data = get_input_datasets(cfg)
    cube_dict = {}
    masks = []
    for (var_type, datasets) in group_metadata(input_data, 'var_type').items():
        print(datasets)
        grouped_datasets = group_metadata(datasets, group_by_attribute)
        for (attr, attr_datasets) in grouped_datasets.items():
            key = _get_key(var_type, attr)
            cube = _get_cube(var_type, group_by_attribute, attr, attr_datasets)
            logger.info("Found cube for '%s'", key)
            cube_dict[key] = cube
            masks.append(np.ma.getmaskarray(cube.data))
    if cfg.get('apply_common_mask'):
        mask = masks[0]
        for new_mask in masks[1:]:
            if new_mask.shape != mask.shape:
                raise ValueError(
                    "Expected datasets with identical shapes when "
                    "'apply_common_mask' is set")
            mask |= new_mask
        for cube in cube_dict.values():
            cube.data = np.ma.array(cube.data, mask=mask)
    return cube_dict


def get_input_datasets(cfg):
    """Get grouped datasets (by tag)."""
    print(cfg.get('patterns')[1])
    input_data = mlr.get_input_data(cfg,
                                    pattern=cfg.get('patterns')[1],
                                    ignore=cfg.get('ignore'))
    tags = list(group_metadata(input_data, 'tag').keys())
    if len(tags) != 1:
        raise ValueError(
            f"Expected unique 'tag' for all input datasets, got {len(tags):d} "
            f"different ones ({tags})")
    return input_data


def get_plot_kwargs(cfg, option, key=None):
    """Get keyword arguments for desired plot function and key."""
    plot_kwargs = cfg.get(option, {}).get('plot_kwargs', {})
    print(plot_kwargs)
    #if 'plot_map' in option:
    #    if 'vmin' in plot_kwargs and 'vmax' in plot_kwargs:
    #        plot_kwargs['levels'] = _get_levels(plot_kwargs['vmin'], plot_kwargs['vmax'])
    #print(plot_kwargs)
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


def plot_xy(cfg, cube_dict):
    """Plot X-Y plots."""
    logger.info("Creating X-Y plots")
    x_coord = cfg['plot_xy'].get('x_coord')
    all_attrs = []

    # Individual plots
    for (key, cube) in cube_dict.items():
        logger.debug("Plotting '%s'", key)
        if cube.ndim != 1:
            raise ValueError(
                f"Expected 1D cube for X-Y plots, got {cube.ndim:d}D cube")
        alias = _get_alias(cfg, key)
        plot_kwargs = get_plot_kwargs(cfg, 'plot_xy', key=key)
        plot_kwargs.setdefault('label', alias)
        _xy_plot(cube, x_coord=x_coord,
                 reg_line=cfg['plot_xy'].get('reg_line', False), **plot_kwargs)
        attrs = cube.attributes
        all_attrs.append(attrs)
        title = _get_title(cfg, alias, attrs)
        plt.title(title)
        plt.ylabel(f'{cube.var_name} / {cube.units}')
        if x_coord is not None:
            coord = cube.coord(x_coord)
            plt.xlabel(f'{coord.var_name} / {coord.units}')
        elif cube.coords(dim_coords=True):
            coord = cube.coord(dim_coords=True)
            plt.xlabel(f'{coord.var_name} / {coord.units}')
        process_pyplot_kwargs(cfg, 'plot_xy')
        plt.legend(**cfg['legend_kwargs'])

        # Save plot
        plot_path = get_plot_filename(f'xy_{key}', cfg)
        savefig_kwargs = get_savefig_kwargs(cfg)
        plt.savefig(plot_path, **savefig_kwargs)
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Provenance
        _write_xy_provenance(cfg, cube, plot_path, title, attrs)

        # Add to global DataFrame
        _add_correlation_information(cfg, title, cube)

    # Merged plot
    logger.debug("Plotting merged plot")
    cubes = iris.cube.CubeList()
    for (key, cube) in cube_dict.items():
        alias = _get_alias(cfg, key)
        plot_kwargs = get_plot_kwargs(cfg, 'plot_xy', key=key)
        plot_kwargs.setdefault('label', alias)
        _xy_plot(cube, x_coord=x_coord,
                 reg_line=cfg['plot_xy'].get('reg_line', False), **plot_kwargs)
        cube = cube.copy()
        ih.prepare_cube_for_merging(cube, key)
        cubes.append(cube)
    cubes = cubes.merge()
    process_pyplot_kwargs(cfg, 'plot_xy')
    plt.legend(**cfg.get('legend_kwargs'))
    plot_path = get_plot_filename('merged_xy', cfg)
    savefig_kwargs = get_savefig_kwargs(cfg)
    plt.savefig(plot_path, **savefig_kwargs)
    logger.info("Wrote %s", plot_path)
    plt.close()
    _write_xy_provenance(cfg, cubes, plot_path, None, *all_attrs)


def plot_scatter_with_errors(cfg, cube_dict):
    """Plot X-Y plots with error range."""
    logger.info("Creating X-Y plots with error ranges")
    keys = {key: key.split(SEP) for key in cube_dict if 'error' not in key}

    # Individual plots
    for (key, split_key) in keys.items():
        cubes = iris.cube.CubeList()
        logger.debug("Plotting '%s'", key)
        plot_kwargs = get_plot_kwargs(cfg, 'plot_xy_with_errors', key=key)
        (cube, error_cube, coord, ancestors) = _xy_plot_with_errors(
            cfg, cube_dict, split_key, **plot_kwargs)

        # Plot appearance
        alias = _get_alias(cfg, key)
        attrs = cube.attributes
        title = _get_title(cfg, alias, attrs)
        plt.title(title)
        plt.ylabel(f'{cube.var_name} / {cube.units}')
        if coord is not None:
            plt.xlabel(f'{coord.var_name} / {coord.units}')
        process_pyplot_kwargs(cfg, 'plot_xy_with_errors')
        plt.legend(**cfg['legend_kwargs'])

        # Save plot
        plot_path = get_plot_filename(f'xy_with_errors_{key}', cfg)
        savefig_kwargs = get_savefig_kwargs(cfg)
        plt.savefig(plot_path, **savefig_kwargs)
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Add to global DataFrame
        _add_correlation_information(cfg, title, cube)

        # Provenance
        cube = cube.copy()
        error_cube = error_cube.copy()
        ih.prepare_cube_for_merging(cube, key)
        ih.prepare_cube_for_merging(error_cube, f'{key}{SEP}error')
        cubes.append(cube)
        cubes.append(error_cube)
        cubes = cubes.merge()
        _write_xy_error_provenance(cfg, cubes, plot_path, title, ancestors)

    # Merged plot
    all_ancestors = []
    cubes = iris.cube.CubeList()
    logger.debug("Plotting merged plot")
    plot_kwargs.pop('color', None)
    for (idx, (key, split_key)) in enumerate(keys.items()):
        plot_kwargs = get_plot_kwargs(cfg, 'plot_xy_with_errors', key=key)
        plot_kwargs['color'] = COLORS[idx]
        (cube, error_cube, _, ancestors) = _xy_plot_with_errors(cfg, cube_dict,
                                                                split_key,
                                                                **plot_kwargs)
        all_ancestors.extend(ancestors)
        cube = cube.copy()
        error_cube = error_cube.copy()
        ih.prepare_cube_for_merging(cube, key)
        ih.prepare_cube_for_merging(error_cube, f'{key}{SEP}error')
        cubes.append(cube)
        cubes.append(error_cube)
    cubes = cubes.merge()
    process_pyplot_kwargs(cfg, 'plot_xy_with_errors')
    plt.legend(**cfg['legend_kwargs'])
    plot_path = get_plot_filename('merged_xy_with_errors', cfg)
    savefig_kwargs = get_savefig_kwargs(cfg)
    plt.savefig(plot_path, **savefig_kwargs)
    logger.info("Wrote %s", plot_path)
    plt.close()
    _write_xy_error_provenance(cfg, cubes, plot_path, None, all_ancestors)


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


def read_data(cfg):
    """Read files and extract vars."""
    logger.debug("Reading data")
    filename1 = cfg['filename1']
    filename2 = cfg['filename2']
    var1 = cfg['var1']
    var2 = cfg['var2']
    logger.debug("Loading %s", filename1)
    cube1 = iris.load_cube(filename)
    logger.debug("Reading %s", filename1)
    cube1 = iris.util.squeeze(cube)
    data1 = {}
    for dataset in cube1['dataset']:
        data1[dataset] = cube1[dataset].data

    logger.debug("Loading %s", filename2)
    cube2 = iris.load_cube(filename)
    logger.debug("Reading %s", filename2)
    cube2 = iris.util.squeeze(cube)
    data2 = {}
    for dataset in cube2['dataset']:
        data2[dataset] = cube2[dataset].data

    for key, value in data2.items():
        if key in data1:
            data1[key].extend(value)
        else:
            logger.debug("Dataset %s is ot part of file %s", key, filename1)
    for key, value in data1.items():
        if len(value) == 1:
            logger.debug("Dataset %s is ot part of file %s", key, filename2)

    return data1


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
      data2.append(np.asscalar(new_cube.data))
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
    #cfg.setdefault('group_by_attribute', 'mlr_model_name')
    #cfg.setdefault('group_attribute_as_default_alias', True)
    cfg.setdefault('legend_kwargs', {})
    cfg.setdefault('map_plot_type', 'pcolormesh')
    cfg.setdefault('print_corr', False)
    cfg.setdefault('years_in_title', False)
    cfg.setdefault('output_file_name', None)

    #cube_dict = get_cube_dict(cfg, cfg['group_by_attribute'])
    x_data, y_data = read_data_and_preprocess(cfg)

    plot_scatter(x_data, y_data, cfg)


# Run main function when this script is called
if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
