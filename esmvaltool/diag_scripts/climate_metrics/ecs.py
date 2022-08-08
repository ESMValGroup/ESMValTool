#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to calculate ECS following Gregory et al. (2004).

Description
-----------
Calculate the equilibrium climate sensitivity (ECS) using the regression method
proposed by Gregory et al. (2004). Further plots related to ECS can be found
in the script ``climate_metrics/feedback_parameters.py``.

If datasets with different numbers of years are given, assume that all data
starts with year 1 in the MMM calculation.

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recipe
-------------------------------
calculate_mmm : bool, optional (default: True)
    Calculate multi-model mean ECS.
complex_gregory_plot : bool, optional (default: False)
    Plot complex Gregory plot (also add response for first ``sep_year`` years
    and last 150 - ``sep_year`` years, default: ``sep_year=20``) if ``True``.
output_attributes : dict, optional
    Write additional attributes to netcdf files.
read_external_file : str, optional
    Read ECS and feedback parameters from external file. The path can be given
    relative to this diagnostic script or as absolute path.
savefig_kwargs : dict, optional
    Keyword arguments for :func:`matplotlib.pyplot.savefig`.
seaborn_settings : dict, optional
    Options for :func:`seaborn.set` (affects all plots).
sep_year : int, optional (default: 20)
    Year to separate regressions of complex Gregory plot. Only effective if
    ``complex_gregory_plot`` is ``True``.
x_lim : list of float, optional (default: [1.5, 6.0])
    Plot limits for X axis of Gregory regression plot (T).
y_lim : list of float, optional (default: [0.5, 3.5])
    Plot limits for Y axis of Gregory regression plot (N).

"""

import logging
import os
from copy import deepcopy
from functools import partial
from pprint import pformat

import cf_units
import iris
import iris.coord_categorisation
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import yaml
from scipy import stats

from esmvaltool.diag_scripts.climate_metrics.feedback_parameters import (
    calculate_anomaly,
)
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    io,
    run_diagnostic,
    select_metadata,
    sorted_metadata,
    variables_available,
)

logger = logging.getLogger(os.path.basename(__file__))

COLORS = sns.color_palette()
EXP_4XCO2 = {
    'CMIP5': 'abrupt4xCO2',
    'CMIP6': 'abrupt-4xCO2',
}
RTMT_DATASETS = set()


def _get_anomaly_data(input_data):
    """Calculate anomaly data for all variables."""
    logger.info("Calculating anomaly data")
    project = input_data[0]['project']
    new_input_data = []
    for (var, var_data) in group_metadata(input_data, 'short_name').items():
        grouped_data = group_metadata(var_data, 'dataset')
        for (dataset_name, datasets) in grouped_data.items():
            logger.debug("Calculating '%s' anomaly for dataset '%s'", var,
                         dataset_name)
            data_4x = select_metadata(datasets, exp=EXP_4XCO2[project])
            data_pic = select_metadata(datasets, exp='piControl')

            # Check if all experiments are available
            if not data_4x:
                raise ValueError(
                    f"No '{EXP_4XCO2[project]}' data available for '{var}' of "
                    f"'{dataset_name}'")
            if not data_pic:
                raise ValueError(
                    f"No 'piControl' data available for '{var}' of "
                    f"'{dataset_name}'")

            # Calculate anomaly, extract correct years and save it
            cube = calculate_anomaly(data_4x, data_pic)
            cube.attributes['project'] = data_4x[0]['project']
            cube.attributes['dataset'] = dataset_name
            if cube.ndim != 1:
                raise ValueError(
                    f"This diagnostic supports only 1D (time), input data, "
                    f"got {cube.ndim}D data")
            new_input_data.append({
                **data_4x[0],
                'ancestors': [data_4x[0]['filename'], data_pic[0]['filename']],
                'cube':
                cube,
            })
    return new_input_data


def _get_data_time_last(cube):
    """Get data of ``cube`` with time axis as last dimension."""
    return np.moveaxis(cube.data, cube.coord_dims('time')[0], -1)


@partial(np.vectorize, excluded=['x_arr'], signature='(n),(n)->()')
def _get_intercept(x_arr, y_arr):
    """Get intercept of linear regression of two (masked) arrays."""
    if np.ma.is_masked(y_arr):
        x_arr = x_arr[~y_arr.mask]
        y_arr = y_arr[~y_arr.mask]
    if len(y_arr) < 2:
        return np.nan
    reg = stats.linregress(x_arr, y_arr)
    return reg.intercept


def _get_multi_model_mean(input_data):
    """Get multi-model mean for all variables."""
    logger.info("Calculating multi-model means")
    project = input_data[0]['project']
    mmm_data = []

    # Iterate over all variables
    for (var, datasets) in group_metadata(input_data, 'short_name').items():
        logger.debug("Calculating multi-model mean for variable '%s'", var)
        ancestors = []
        dataset_names = []
        mmm = {}

        # Read data from every datasets
        for dataset in datasets:
            try:
                cube = dataset['cube']
            except KeyError:
                raise KeyError(
                    f"No data for '{var}' of dataset '{dataset['dataset']}' "
                    f"for multi-model mean calculation")
            if cube.ndim > 1:
                raise ValueError(
                    f"Calculation of multi-model mean not supported for input "
                    f"data with more than one dimension (which should be "
                    f"time), got {cube.ndim:d}-dimensional cube")
            ancestors.extend(dataset['ancestors'])
            dataset_names.append(dataset['dataset'])
            mmm[dataset['dataset']] = cube.data

        # Adapt shapes if necessary
        target_shape = max([d.shape[0] for d in mmm.values()])
        for (dataset_name, dataset_data) in mmm.items():
            if dataset_data.shape[0] != target_shape:
                dataset_data = np.pad(
                    dataset_data, (0, target_shape - dataset_data.shape[0]),
                    constant_values=np.nan)
                mmm[dataset_name] = dataset_data

        # Calculate MMM
        mmm = np.ma.masked_invalid(list(mmm.values()))
        mmm_cube = cube.copy(data=np.ma.mean(mmm, axis=0))
        attributes = {
            'ancestors': ancestors,
            'dataset': 'MultiModelMean',
            'datasets': '|'.join(dataset_names),
            'project': project,
            'short_name': var,
        }
        mmm_cube.attributes = attributes
        mmm_data.append({**attributes, 'cube': mmm_cube})
    input_data.extend(mmm_data)
    return input_data


@partial(np.vectorize, excluded=['x_arr'], signature='(n),(n)->()')
def _get_slope(x_arr, y_arr):
    """Get slope of linear regression of two (masked) arrays."""
    if np.ma.is_masked(y_arr):
        x_arr = x_arr[~y_arr.mask]
        y_arr = y_arr[~y_arr.mask]
    if len(y_arr) < 2:
        return np.nan
    reg = stats.linregress(x_arr, y_arr)
    return reg.slope


def _plot_complex_gregroy_plot(cfg, axes, tas_cube, rtnt_cube, reg_all):
    """Plot complex Gregory plot."""
    sep = cfg['sep_year']

    # Regressions
    x_reg = np.linspace(cfg['x_lim'][0] - 1.0, cfg['x_lim'][1] + 1.0, 2)
    reg_first = stats.linregress(tas_cube.data[:sep], rtnt_cube.data[:sep])
    reg_last = stats.linregress(tas_cube.data[sep:], rtnt_cube.data[sep:])
    y_reg_first = reg_first.slope * x_reg + reg_first.intercept
    y_reg_last = reg_last.slope * x_reg + reg_last.intercept
    y_reg_all = reg_all.slope * x_reg + reg_all.intercept
    ecs_first = -reg_first.intercept / (2.0 * reg_first.slope)
    ecs_last = -reg_last.intercept / (2.0 * reg_last.slope)
    ecs_all = -reg_all.intercept / (2.0 * reg_all.slope)

    # Plots
    axes.scatter(tas_cube.data[:sep],
                 rtnt_cube.data[:sep],
                 color=COLORS[0],
                 marker='o',
                 s=8,
                 alpha=0.7,
                 label=f'first {sep:d} years: ECS = {ecs_first:.2f} K')
    axes.scatter(tas_cube.data[sep:],
                 rtnt_cube.data[sep:],
                 color=COLORS[1],
                 marker='o',
                 s=8,
                 alpha=0.7,
                 label=f'last {tas_cube.shape[0] - sep:d} years: ECS = '
                 f'{ecs_last:.2f} K')
    axes.plot(x_reg, y_reg_first, color=COLORS[0], linestyle='-', alpha=0.6)
    axes.plot(x_reg, y_reg_last, color=COLORS[1], linestyle='-', alpha=0.6)
    axes.plot(x_reg,
              y_reg_all,
              color='k',
              linestyle='-',
              alpha=0.9,
              label=r'all years: ECS = {:.2f} K ($R^2$ = {:.2f})'.format(
                  ecs_all, reg_all.rvalue**2))

    # Legend
    return axes.legend(loc='upper right')


def _write_ecs_regression(cfg, tas_cube, rtnt_cube, reg_stats, dataset_name):
    """Write Gregory regression cube."""
    ecs = -reg_stats.intercept / (2.0 * reg_stats.slope)
    attrs = {
        'anomaly': 'relative to piControl run',
        'regression_r_value': reg_stats.rvalue,
        'regression_slope': reg_stats.slope,
        'regression_interception': reg_stats.intercept,
        'Climate Feedback Parameter': reg_stats.slope,
        'ECS': ecs,
    }
    attrs.update(cfg.get('output_attributes', {}))
    cubes = iris.cube.CubeList()
    for cube in [tas_cube, rtnt_cube]:
        cube.var_name += '_anomaly'
        cube.long_name += ' Anomaly'
        cube.attributes = attrs
        cubes.append(cube)
    netcdf_path = get_diagnostic_filename('ecs_regression_' + dataset_name,
                                          cfg)
    io.iris_save(cubes, netcdf_path)
    return netcdf_path


def check_input_data(cfg):
    """Check input data."""
    if not variables_available(cfg, ['tas']):
        raise ValueError(
            "This diagnostic needs variable 'tas' if 'read_external_file' is "
            "not given")
    if not (variables_available(cfg, ['rtnt'])
            or variables_available(cfg, ['rtmt'])):
        raise ValueError(
            "This diagnostic needs the variable 'rtnt' or 'rtmt' if "
            "'read_external_file' is not given")
    input_data = cfg['input_data'].values()
    project_group = group_metadata(input_data, 'project')
    projects = list(project_group.keys())
    if len(projects) > 1:
        raise ValueError(
            f"This diagnostic supports only unique 'project' attributes, got "
            f"{projects}")
    project = projects[0]
    if project not in EXP_4XCO2:
        raise ValueError(f"Project '{project}' not supported yet")
    exp_group = group_metadata(input_data, 'exp')
    exps = set(exp_group.keys())
    if exps != {'piControl', EXP_4XCO2[project]}:
        raise ValueError(
            f"This diagnostic needs 'piControl' and '{EXP_4XCO2[project]}' "
            f"experiments, got {exps}")


def preprocess_data(cfg):
    """Extract input data."""
    input_data = deepcopy(list(cfg['input_data'].values()))
    input_data = sorted_metadata(input_data, ['short_name', 'exp', 'dataset'])
    if not input_data:
        return ([], [])

    # Use 'rtmt' instead of 'rtmt' if necessary
    for dataset in input_data:
        if dataset['short_name'] == 'rtmt':
            RTMT_DATASETS.add(dataset['dataset'])
            dataset['short_name'] = 'rtnt'
    if RTMT_DATASETS:
        logger.info("Using 'rtmt' instead of 'rtnt' for datasets '%s'",
                    RTMT_DATASETS)

    # Calculate anomalies for every dataset
    input_data = _get_anomaly_data(input_data)

    # Calculate multi-model mean
    if cfg.get('calculate_mmm', True):
        input_data = _get_multi_model_mean(input_data)

    # Group data in terms of dataset
    tas_data = select_metadata(input_data, short_name='tas')
    rtnt_data = select_metadata(input_data, short_name='rtnt')
    tas_data = group_metadata(tas_data, 'dataset')
    rtnt_data = group_metadata(rtnt_data, 'dataset')
    return (tas_data, rtnt_data)


def get_provenance_record(caption):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption': caption,
        'statistics': ['mean', 'diff'],
        'domains': ['global'],
        'authors': ['schlund_manuel'],
        'references': ['gregory04grl'],
        'realms': ['atmos'],
        'themes': ['phys'],
    }
    return record


def read_external_file(cfg):
    """Read external file to get ECS."""
    filepath = os.path.expanduser(os.path.expandvars(
        cfg['read_external_file']))
    if not os.path.isabs(filepath):
        filepath = os.path.join(os.path.dirname(__file__), filepath)
    if not os.path.isfile(filepath):
        raise FileNotFoundError(
            f"Desired external file '{filepath}' does not exist")
    with open(filepath, 'r') as infile:
        external_data = yaml.safe_load(infile)
    ecs = external_data.get('ecs', {})
    feedback_parameter = external_data.get('feedback_parameter', {})
    logger.info("Reading external file '%s'", filepath)
    logger.info("Found ECS (K):")
    logger.info("%s", pformat(ecs))
    logger.info("Found climate feedback parameters (W m-2 K-1):")
    logger.info("%s", pformat(feedback_parameter))
    return (ecs, feedback_parameter, filepath)


def plot_gregory_plot(cfg, dataset_name, tas_cube, rtnt_cube, reg_stats):
    """Plot linear regression used to calculate ECS."""
    (_, axes) = plt.subplots()
    ecs = -reg_stats.intercept / (2 * reg_stats.slope)
    project = tas_cube.attributes['project']

    # Regression line
    x_reg = np.linspace(cfg['x_lim'][0] - 1.0, cfg['x_lim'][1] + 1.0, 2)
    y_reg = reg_stats.slope * x_reg + reg_stats.intercept

    # Plot data
    if cfg.get('complex_gregory_plot'):
        legend = _plot_complex_gregroy_plot(cfg, axes, tas_cube, rtnt_cube,
                                            reg_stats)
    else:
        axes.scatter(tas_cube.data,
                     rtnt_cube.data,
                     color=COLORS[0],
                     marker='o',
                     s=8,
                     alpha=0.7)
        legend = None
        axes.plot(x_reg, y_reg, color='k', linestyle='-', alpha=0.8)
        axes.text(
            0.05,
            0.9,
            r'R$^2$ = {:.2f}, ECS = {:.2f} K'.format(reg_stats.rvalue**2, ecs),
            transform=axes.transAxes,
        )
    axes.axhline(0.0, color='gray', linestyle=':')

    # Plot appearance
    axes.set_title(f"Gregory regression for {dataset_name} ({project})")
    axes.set_xlabel("ΔT [K]")
    axes.set_ylabel(r"N [W m$^{-2}$]")
    axes.set_xlim(cfg['x_lim'])
    axes.set_ylim(cfg['y_lim'])

    # Save plot
    plot_path = get_plot_filename(dataset_name, cfg)
    plt.savefig(plot_path,
                additional_artists=[legend],
                **cfg['savefig_kwargs'])
    logger.info("Wrote %s", plot_path)
    plt.close()

    # Write netcdf file for every plot
    netcdf_path = _write_ecs_regression(cfg, tas_cube, rtnt_cube, reg_stats,
                                        dataset_name)

    # Provenance
    provenance_record = get_provenance_record(
        f"Scatterplot between TOA radiance and global mean surface "
        f"temperature anomaly for 150 years of the abrupt 4x CO2 experiment "
        f"including linear regression to calculate ECS for {dataset_name} "
        f"({project}).")
    provenance_record.update({
        'plot_types': ['scatter'],
    })

    return (netcdf_path, plot_path, provenance_record)


def set_default_cfg(cfg):
    """Set default values for cfg."""
    cfg = deepcopy(cfg)
    cfg.setdefault('savefig_kwargs', {
        'dpi': 300,
        'orientation': 'landscape',
        'bbox_inches': 'tight',
    })
    cfg.setdefault('sep_year', 20)
    cfg.setdefault('x_lim', [0.0, 12.0])
    cfg.setdefault('y_lim', [-2.0, 10.0])
    return cfg


def write_data(cfg, ecs_data, feedback_parameter_data, ancestor_files):
    """Write netcdf files."""
    data = [ecs_data, feedback_parameter_data]
    var_attrs = [
        {
            'short_name': 'ecs',
            'long_name': 'Equilibrium Climate Sensitivity (Gregory method)',
            'units': cf_units.Unit('K'),
        },
        {
            'short_name': 'lambda',
            'long_name': 'Climate Feedback Parameter',
            'units': cf_units.Unit('W m-2 K-1'),
        },
    ]
    input_data = list(cfg['input_data'].values())
    if input_data:
        attrs = {
            'project': input_data[0]['project'],
        }
    else:
        attrs = {}
    if RTMT_DATASETS:
        rtmt_datasets = sorted(list(RTMT_DATASETS))
        attrs['net_toa_radiation'] = (
            f"For datasets {rtmt_datasets}, 'rtmt' (net top of model "
            f"radiation) instead of 'rtnt' (net top of atmosphere radiation) "
            f"is used due to lack of data. These two variables might differ.")
    attrs.update(cfg.get('output_attributes', {}))
    data_available = False
    for (idx, var_attr) in enumerate(var_attrs):
        if not data[idx]:
            logger.info(
                "Skipping writing of '%s' for all models, no data available",
                var_attr['short_name'])
            continue
        data_available = True
        path = get_diagnostic_filename(var_attr['short_name'], cfg)
        io.save_scalar_data(data[idx], path, var_attr, attributes=attrs)
        caption = "{long_name} for multiple climate models.".format(**var_attr)
        provenance_record = get_provenance_record(caption)
        provenance_record['ancestors'] = ancestor_files
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(path, provenance_record)
    if not data_available:
        raise ValueError("No input data given")


def main(cfg):
    """Run the diagnostic."""
    cfg = set_default_cfg(cfg)
    sns.set(**cfg.get('seaborn_settings', {}))

    # Read external file if desired
    if cfg.get('read_external_file'):
        (ecs, feedback_parameter, external_file) = read_external_file(cfg)
    else:
        check_input_data(cfg)
        ecs = {}
        feedback_parameter = {}
        external_file = None

    # Read and preprocess data
    all_ancestors = []
    (tas_data, rtnt_data) = preprocess_data(cfg)

    # Iterate over all datasets and save ECS and feedback parameter
    for dataset_name in tas_data:
        logger.info("Processing '%s'", dataset_name)
        if dataset_name not in rtnt_data:
            raise ValueError(
                f"No 'rtmt' or 'rtnt' data for '{dataset_name}' available")
        tas_cube = tas_data[dataset_name][0]['cube']
        rtnt_cube = rtnt_data[dataset_name][0]['cube']
        ancestor_files = (tas_data[dataset_name][0]['ancestors'] +
                          rtnt_data[dataset_name][0]['ancestors'])

        # Perform linear regression
        reg = stats.linregress(tas_cube.data, rtnt_cube.data)

        # Plot Gregory plots
        (path, plot_path, provenance_record) = plot_gregory_plot(
            cfg, dataset_name, tas_cube, rtnt_cube, reg)

        # Provenance
        provenance_record['ancestors'] = ancestor_files
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(path, provenance_record)
            provenance_logger.log(plot_path, provenance_record)

        # Save data
        if cfg.get('read_external_file') and dataset_name in ecs:
            logger.warning(
                "Overwriting externally given ECS and climate feedback "
                "parameter from file '%s' for '%s'", external_file,
                dataset_name)
        ecs[dataset_name] = -reg.intercept / (2 * reg.slope)
        feedback_parameter[dataset_name] = reg.slope
        all_ancestors.extend(ancestor_files)

    # Write data
    if external_file is not None:
        all_ancestors.append(external_file)
    all_ancestors = list(set(all_ancestors))
    write_data(cfg, ecs, feedback_parameter, all_ancestors)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
