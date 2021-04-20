#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Rescale label data using a single emergent constraint.

Description
-----------
This diagnostic uses an emergent relationship between data marked as
``var_type=label`` (Y axis) and ``var_type=feature`` (X axis) together with an
observation of the X axis (``var_type=prediction_input`` and
``var_type=prediction_input_error``) to calculate factors that are necessary to
rescale each input point so that it matches the constraint. The rescaling is
applied to data marked as ``var_type=label_to_rescale``. All data needs the
attribute ``tag`` which needs to be identical for ``label``,
``prediction_input``, ``prediction_input_error`` and ``label_to_rescale``. Only
a single ``tag`` for ``feature`` is possible.

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recipe
-------------------------------
group_by_attributes: list of str, optional (default: ['dataset'])
    List of attributes used to separate different input points.
ignore: list of dict, optional
    Ignore specific datasets by specifying multiple :obj:`dict` s of metadata.
legend_kwargs: dict, optional
    Optional keyword arguments of :func:`matplotlib.pyplot.legend` (affects
    only plots with legends).
pattern: str, optional
    Pattern matched against ancestor file names.
plot_emergent_relationship: dict, optional
    If given, plot emergent relationship between X and Y data. Specify
    additional keyword arguments by ``plot_kwargs`` and plot appearance options
    by ``pyplot_kwargs`` (processed as functions of :mod:`matplotlib.pyplot`).
    Use ``{}`` to plot with default settings.
plot_kwargs_for_groups: dict, optional
    Specify additional keyword arguments (values) for the different points
    defined by ``group_by_attributes`` (keys) used in plots.
savefig_kwargs: dict, optional
    Keyword arguments for :func:`matplotlib.pyplot.savefig`.
seaborn_settings: dict, optional
    Options for :func:`seaborn.set` (affects all plots).

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
import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts import mlr
from esmvaltool.diag_scripts.mlr.plot import get_savefig_kwargs
from esmvaltool.diag_scripts.mlr.preprocess import cube_to_aux_coord
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    io,
    run_diagnostic,
    select_metadata,
)

logger = logging.getLogger(os.path.basename(__file__))


GROUPS_SEP = '__'
UNITS_SEP = '___'


def _check_datasets(datasets, name):
    """Check input data."""
    if not datasets:
        raise ValueError(f"No data with var_type '{name}' given")
    keys_to_check = ['tag', 'short_name', 'long_name', 'units']
    for key in keys_to_check:
        vals = list(group_metadata(datasets, key).keys())
        if len(vals) != 1:
            raise ValueError(
                f"Expected data with unique '{key}' for var_type "
                f"'{name}', got {vals}")


def _get_data(dataset):
    """Get data from netcdf file."""
    cube = iris.load_cube(dataset['filename'])
    if cube.data.size != 1:
        raise ValueError(
            f"Expected scalar data for dataset {dataset}, got data with shape "
            f"{cube.shape}")
    return cube.data.squeeze()


def _get_data_frame(datasets, group_by_attributes, var_type):
    """Extract :class:`pandas.DataFrame` from :obj:`list` of datasets."""
    tag = f"{datasets[0]['tag']}{UNITS_SEP}{datasets[0]['units']}"
    data_frame = pd.DataFrame(columns=[tag])
    for dataset in datasets:
        new_group = _get_group(dataset, group_by_attributes)
        if new_group in data_frame.index:
            raise ValueError(
                f"Got duplicate data for group '{new_group}' of var_type "
                f"'{var_type}', consider extending list of attributes that "
                f"is used for grouping (currently: {group_by_attributes})")
        new_data = _get_data(dataset)
        data_frame.loc[new_group] = new_data
    data_frame = data_frame.sort_index()
    data_frame.index.name = 'dataset'
    return data_frame


def _get_datasets_for_ec(input_data):
    """Check input data."""
    features = select_metadata(input_data, var_type='feature')
    labels = select_metadata(input_data, var_type='label')
    pred_input = select_metadata(input_data, var_type='prediction_input')
    pred_input_err = select_metadata(input_data,
                                     var_type='prediction_input_error')
    data_to_check = {
        'feature': features,
        'label': labels,
        'prediction_input': pred_input,
        'prediction_input_error': pred_input_err,
    }
    for (name, data) in data_to_check.items():
        _check_datasets(data, name)
    return (features, labels, pred_input, pred_input_err)


def _get_ec_ancestors(cfg):
    """Get ancestor files for emergent constraint."""
    input_data = _get_input_data(cfg)
    ancestors = []
    for var_type in ('feature', 'label', 'prediction_input',
                     'prediction_input_error'):
        datasets = select_metadata(input_data, var_type=var_type)
        ancestors.extend([d['filename'] for d in datasets])
    return ancestors


def _get_ec_cube(x_data, y_data):
    """Get :class:`iris.cube.Cube` representing emergent relationship."""
    (feature, feature_units, label, label_units) = _get_tags(x_data, y_data)
    x_cube = ec.pandas_object_to_cube(x_data, var_name=feature,
                                      units=feature_units)[:, 0]
    x_coord = cube_to_aux_coord(x_cube)
    y_cube = ec.pandas_object_to_cube(y_data, var_name=label,
                                      units=label_units)[:, 0]
    y_cube.add_aux_coord(x_coord, 0)
    y_cube.remove_coord('columns')
    y_cube.attributes['dataset'] = ''
    y_cube.attributes['project'] = ''
    return y_cube


def _get_error_dataset(cfg, datasets):
    """Get error dataset."""
    error_dataset = {}
    for key in mlr.NECESSARY_KEYS:
        vals = sorted(list({str(d[key]) for d in datasets}))
        error_dataset[key] = '|'.join(vals)
    start_years = list({d['start_year'] for d in datasets})
    end_years = list({d['end_year'] for d in datasets})
    error_dataset['start_year'] = min(start_years)
    error_dataset['end_year'] = max(end_years)
    error_dataset['standard_name'] = None
    error_dataset['short_name'] += '_standard_error'
    error_dataset['long_name'] += ' (Standard Error)'
    error_dataset['var_type'] = 'prediction_output_error'
    error_dataset['error'] = 'due to rescaling using emergent relationship'
    error_dataset['filename'] = get_diagnostic_filename(
        'standard_error_of_emergent_constraint', cfg)
    return error_dataset


def _get_group(dataset, group_by_attributes):
    """Get name of group."""
    values = []
    for attr in group_by_attributes:
        if attr not in dataset:
            raise KeyError(
                f"Attribute '{attr}' not available in dataset {dataset}")
        values.append(dataset[attr])
    return GROUPS_SEP.join(values)


def _get_input_data(cfg):
    """Get input data."""
    input_data = mlr.get_input_data(cfg,
                                    pattern=cfg.get('pattern'),
                                    ignore=cfg.get('ignore'))
    return input_data


def _get_plot_kwargs(cfg, option, group=None):
    """Get plot keyword arguments for a group."""
    plot_kwargs = cfg.get(option, {}).get('plot_kwargs', {})
    plot_kwargs = deepcopy(plot_kwargs)
    if group is None:
        return plot_kwargs
    group_plot_kwargs = cfg.get('plot_kwargs_for_groups', {}).get(group, {})
    plot_kwargs.update(group_plot_kwargs)
    plot_kwargs.setdefault('linestyle', '-')
    if plot_kwargs['linestyle'] == '-':
        plot_kwargs.setdefault('marker', 'o')
    else:
        plot_kwargs.setdefault('marker', 's')
    return plot_kwargs


def _get_ref_cube(datasets):
    """Get (unique) shape of datasets."""
    cubes = [iris.load_cube(d['filename']) for d in datasets]
    shapes = list({cube.shape for cube in cubes})
    if len(shapes) != 1:
        raise ValueError(
            f"Expected unique shape for 'label_to_rescale' data, got {shapes}")
    ref_cube = cubes[0]
    ref_cube.attributes = {}
    return ref_cube


def _get_mmm_cube(datasets):
    """Extract data."""
    cubes = iris.cube.CubeList()
    cube_labels = []
    ref_cube = iris.load_cube(datasets[0]['filename'])
    for (idx, dataset) in enumerate(datasets):
        path = dataset['filename']
        cube = iris.load_cube(path)
        ih.prepare_cube_for_merging(cube, str(idx))
        cubes.append(cube)
        cube_labels.append(str(idx))
    mmm_cube = cubes.merge_cube()
    if len(cube_labels) > 1:
        mmm_cube = mmm_cube.collapsed(['cube_label'], iris.analysis.MEAN)
    for aux_coord in ref_cube.coords(dim_coords=False):
        mmm_cube.add_aux_coord(aux_coord, ref_cube.coord_dims(aux_coord))
    mmm_cube.remove_coord('cube_label')
    return mmm_cube


def _get_tags(x_data, y_data):
    """Extract tags from X and Y data."""
    feature = x_data.columns[0].split(UNITS_SEP)[0]
    feature_units = x_data.columns[0].split(UNITS_SEP)[1]
    label = y_data.columns[0].split(UNITS_SEP)[0]
    label_units = y_data.columns[0].split(UNITS_SEP)[1]
    return (feature, feature_units, label, label_units)


def _process_pyplot_kwargs(cfg, option):
    """Process functions for :mod:`matplotlib.pyplot`."""
    for (key, val) in cfg.get(option, {}).get('pyplot_kwargs', {}).items():
        getattr(plt, key)(val)


def get_constraint(x_data, y_data, x_ref, x_ref_err):
    """Print constraint value for Y axis."""
    (feature, feature_units, label, label_units) = _get_tags(x_data, y_data)
    x_data = x_data.values.squeeze()
    y_data = y_data.values.squeeze()
    x_ref = x_ref.values.squeeze()
    x_ref_err = x_ref_err.values.squeeze()
    (y_data_lin, y_pdf) = ec.target_pdf(x_data, y_data, x_ref, x_ref_err)
    y_mean = np.sum(y_data_lin * y_pdf) / np.sum(y_pdf)
    y_var = np.sum((y_data_lin - y_mean)**2 * y_pdf) / np.sum(y_pdf)
    y_std = np.sqrt(y_var)
    lines = ec.regression_line(x_data, y_data)
    logger.info("Observational constraint on '%s': (%.3f ± %.3f) %s", feature,
                x_ref, x_ref_err, feature_units)
    logger.info("Constraint on target variable '%s': (%.3f ± %.3f) %s", label,
                y_mean, y_std, label_units)
    logger.info("R2 of emergent relationship: %.3f (p = %.4f)",
                lines['rvalue']**2, lines['pvalue'])
    return (y_mean, y_std)


def get_emergent_constraint_data(cfg):
    """Get :class:`pandas.DataFrame` that contains the data."""
    input_data = _get_input_data(cfg)
    (features, labels, pred_input,
     pred_input_err) = _get_datasets_for_ec(input_data)

    # Extract data frames
    x_data = _get_data_frame(features, cfg['group_by_attributes'], 'feature')
    y_data = _get_data_frame(labels, cfg['group_by_attributes'], 'label')
    x_ref = _get_data_frame(pred_input, cfg['group_by_attributes'],
                            'prediction_input')
    x_ref_err = _get_data_frame(pred_input_err, cfg['group_by_attributes'],
                                'prediction_input_error')

    # Check data frames
    if len(x_data.index) < 2:
        raise ValueError(
            f"Expected at least two input points for X data, got "
            f"{len(x_data.index):d}")
    if not x_data.index.equals(y_data.index):
        raise ValueError(
            f"Expected identical input points for X and Y data, got\nX: "
            f"{x_data.index.values}\nY: {y_data.index.values}")
    if len(x_ref.index) != 1:
        raise ValueError(
            f"Expected exactly one prediction input point for X data, got "
            f"{len(x_ref.index):d}")
    if not x_ref.index.equals(x_ref_err.index):
        raise ValueError(
            f"Expected identical input points for prediction input and its "
            f"corresponding errors for X data, got {x_ref.index.values} and "
            f"{x_ref_err.index.values}, respectively")
    logger.info("Found X data:\n%s", x_data)
    logger.info("Found Y data:\n%s", y_data)
    logger.info("Found X reference data:\n%s", x_ref)
    logger.info("Found X reference error data:\n%s", x_ref_err)
    return (x_data, y_data, x_ref, x_ref_err)


def plot_emergent_relationship(cfg, x_data, y_data, x_ref, x_ref_err, y_mean):
    """Plot emergent relationship."""
    (feature, feature_units, label, label_units) = _get_tags(x_data, y_data)
    logger.info("Plotting emergent relationship between '%s' and '%s'",
                label, feature)
    (_, axes) = plt.subplots()

    # Plot data points
    for group in x_data.index:
        plot_kwargs = _get_plot_kwargs(cfg, 'plot_emergent_relationship',
                                       group=group)
        plot_kwargs['linestyle'] = 'none'
        plot_kwargs['label'] = group
        axes.plot(x_data.loc[group], y_data.loc[group], **plot_kwargs)

    # Plot regression lines
    axes.set_xlim(auto=False)
    axes.set_ylim(auto=False)
    lines = ec.regression_line(x_data.values.squeeze(),
                               y_data.values.squeeze())
    lines['x'] = np.squeeze(lines['x'])
    axes.plot(lines['x'],
              lines['y'],
              color='orange',
              linestyle='-',
              label='Linear regression')
    axes.fill_between(lines['x'],
                      lines['y_minus_err'],
                      lines['y_plus_err'],
                      color='orange',
                      alpha=0.2)

    # Plot reference
    x_ref = x_ref.values.squeeze()
    x_ref_err = x_ref_err.values.squeeze()
    axes.axvline(x_ref,
                 color='k',
                 linestyle=':',
                 label='Observational constraint')
    axes.axvspan(x_ref - x_ref_err, x_ref + x_ref_err, color='k', alpha=0.1)
    axes.axhline(y_mean, color='k', linestyle=':')

    # Plot appearance
    axes.set_title(f"Emergent relationship between {label} and {feature}")
    axes.set_xlabel(f"{feature} [{feature_units}]")
    axes.set_ylabel(f"{label} [{label_units}]")
    _process_pyplot_kwargs(cfg, 'plot_emergent_relationship')
    plt.legend(**cfg['legend_kwargs'])
    text = rf"$R^2$ = {lines['rvalue']**2:.2f}, p = {lines['pvalue']:.3f}"
    if lines['rvalue'] > 0.0:
        axes.text(0.6, 0.05, text, transform=axes.transAxes)
    else:
        axes.text(0.6, 0.95, text, transform=axes.transAxes)

    # Save plot
    plot_path = get_plot_filename(f'{label}_vs_{feature}', cfg)
    savefig_kwargs = get_savefig_kwargs(cfg)
    plt.savefig(plot_path, **savefig_kwargs)
    logger.info("Wrote %s", plot_path)
    plt.close()

    # Provenance
    cube = _get_ec_cube(x_data, y_data)
    netcdf_path = get_diagnostic_filename(f'{label}_vs_{feature}', cfg)
    io.iris_save(cube, netcdf_path)
    record = {
        'ancestors': _get_ec_ancestors(cfg),
        'authors': ['schlund_manuel'],
        'caption': f"Emergent relationship between {label} and {feature}.",
        'plot_file': plot_path,
        'plot_types': ['scatter'],
        'references': ['schlund20jgr'],
        'themes': ['EC'],
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, record)


def rescale_labels(cfg, y_data, y_mean, y_std):
    """Rescale labels."""
    input_data = _get_input_data(cfg)
    labels_to_rescale = select_metadata(input_data,
                                        var_type='label_to_rescale')
    _check_datasets(labels_to_rescale, 'label_to_rescale')

    # Get groups
    groups = []
    for dataset in labels_to_rescale:
        group = _get_group(dataset, cfg['group_by_attributes'])
        groups.append(group)
        dataset['group__for__rescaling'] = group

    groups.sort()
    if set(groups) != set(y_data.index):
        raise ValueError(
            f"Expected identical groups for 'label' and 'label_to_rescale' "
            f"data, got\n'label': {y_data.index.values}\n'label_to_rescale': "
            f"{np.array(groups)}")

    # Rescale data
    ref_cube = _get_ref_cube(labels_to_rescale)
    for dataset in labels_to_rescale:
        cube = iris.load_cube(dataset['filename'])
        rescaling_factor = (
            y_mean / y_data.loc[dataset['group__for__rescaling']].values)
        logger.info("Rescaling '%s' with factor %.2f",
                    dataset['group__for__rescaling'], rescaling_factor)
        rescaled_cube = cube.copy(cube.data * rescaling_factor)

        # Adapt metadata
        rescaled_dataset = deepcopy(dataset)
        rescaled_dataset['var_type'] = 'label'
        rescaled_dataset['rescaled'] = 'using emergent relationship'
        if '_label' in dataset['filename']:
            rescaled_dataset['filename'] = dataset['filename'].replace(
                '_label_to_rescale', '_rescaled_label')
        else:
            rescaled_dataset['filename'] = dataset['filename'].replace(
                '.nc', '_rescaled_label.nc')

        # Save data
        rescaled_dataset['filename'] = mlr.get_new_path(
            cfg, rescaled_dataset['filename'])
        io.metadata_to_netcdf(rescaled_cube, rescaled_dataset)

        # Provenance
        record = {
            'ancestors': [dataset['filename']] + _get_ec_ancestors(cfg),
            'authors': ['schlund_manuel'],
            'caption': f"Rescaled {rescaled_cube.long_name} for "
                       f"{mlr.get_alias(rescaled_dataset)} using emergent "
                       f"relationship.",
            'references': ['schlund20jgr'],
            'themes': ['EC'],
        }
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(rescaled_dataset['filename'], record)

    # Rescale MMM to estimate error
    logger.debug("Estimating error using global error %e", y_std)
    mmm_cube = _get_mmm_cube(labels_to_rescale)
    error_cube = ref_cube.copy(mmm_cube.data * y_std / y_data.mean().values)
    error_dataset = _get_error_dataset(cfg, labels_to_rescale)
    io.metadata_to_netcdf(error_cube, error_dataset)

    # Provenance
    record = {
        'ancestors': ([d['filename'] for d in labels_to_rescale] +
                      _get_ec_ancestors(cfg)),
        'authors': ['schlund_manuel'],
        'caption': f"Rescaled {error_cube.long_name} using emergent "
                   f"relationship.",
        'references': ['schlund20jgr'],
        'themes': ['EC'],
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(error_dataset['filename'], record)


def main(cfg):
    """Run the diagnostic."""
    sns.set(**cfg.get('seaborn_settings', {}))
    cfg = deepcopy(cfg)
    cfg.setdefault('group_by_attributes', ['dataset'])
    cfg.setdefault('legend_kwargs', {})
    logger.info("Using attributes %s to group input data",
                cfg['group_by_attributes'])

    # Extract data
    (x_data, y_data, x_ref, x_ref_err) = get_emergent_constraint_data(cfg)

    # Get constraint
    (y_mean, y_std) = get_constraint(x_data, y_data, x_ref, x_ref_err)

    # Plots
    if 'plot_emergent_relationship' in cfg:
        plot_emergent_relationship(cfg, x_data, y_data, x_ref, x_ref_err,
                                   y_mean)

    # Rescale labels
    rescale_labels(cfg, y_data, y_mean, y_std)


# Run main function when this script is called
if __name__ == '__main__':
    mlr.ignore_warnings()
    with run_diagnostic() as config:
        main(config)
