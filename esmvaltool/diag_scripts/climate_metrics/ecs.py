#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to calculate ECS following Gregory et al. (2004).

Description
-----------
Calculate the equilibrium climate sensitivity (ECS) using the regression method
proposed by Gregory et al. (2004).

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
output_attributes : dict, optional
    Write additional attributes to netcdf files.
read_external_file : str, optional
    Read ECS and feedback parameters from external file. The path can be given
    relative to this diagnostic script or as absolute path.
seaborn_settings : dict, optional
    Options for :func:`seaborn.set` (affects all plots), see
    <https://seaborn.pydata.org/generated/seaborn.set.html>.

"""

import logging
import os
from copy import deepcopy
from functools import partial
from pprint import pformat

import cf_units
import iris
import iris.coord_categorisation
import numpy as np
import seaborn as sns
import yaml
from scipy import stats

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger, extract_variables, get_diagnostic_filename,
    get_plot_filename, group_metadata, io, plot, run_diagnostic,
    select_metadata, variables_available)

logger = logging.getLogger(os.path.basename(__file__))

EXP_4XCO2 = {
    'CMIP5': 'abrupt4xCO2',
    'CMIP6': 'abrupt-4xCO2',
}
RTMT_DATASETS = set()


def _calculate_anomaly(data_4x, data_pic):
    """Calculate anomaly cube for a dataset."""
    cube_4x = iris.load_cube(data_4x[0]['filename'])
    iris.coord_categorisation.add_year(cube_4x, 'time')
    cube_4x = cube_4x.aggregated_by('year', iris.analysis.MEAN)

    cube_pic = iris.load_cube(data_pic[0]['filename'])
    iris.coord_categorisation.add_year(cube_pic, 'time')
    cube_pic = cube_pic.aggregated_by('year', iris.analysis.MEAN)

    x_data = cube_pic.coord('year').points
    y_data = _get_data_time_last(cube_pic)
    slope = _get_slope(x_data, y_data)
    intercept = _get_intercept(x_data, y_data)
    for _ in range(cube_pic.ndim - 1):
        x_data = np.expand_dims(x_data, -1)
    new_x_data = np.broadcast_to(x_data, cube_pic.shape)
    new_data = slope * new_x_data + intercept
    cube_4x.data -= np.ma.masked_invalid(new_data)
    return cube_4x


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
            cube = _calculate_anomaly(data_4x, data_pic)
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
    for (var, datasets) in group_metadata(input_data, 'short_name').items():
        logger.debug("Calculating multi-model mean for variable '%s'", var)
        ancestors = []
        dataset_names = []
        mmm = []
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
            mmm.append(cube.data)
        mmm = np.ma.array(mmm)
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


def plot_ecs_regression(cfg, dataset_name, tas_cube, rtnt_cube, reg_stats):
    """Plot linear regression used to calculate ECS."""
    if not cfg['write_plots']:
        return (None, None)
    ecs = -reg_stats.intercept / (2 * reg_stats.slope)

    # Regression line
    x_reg = np.linspace(-1.0, 9.0, 2)
    y_reg = reg_stats.slope * x_reg + reg_stats.intercept

    # Plot data
    text = r'r = {:.2f}, $\lambda$ = {:.2f}, F = {:.2f}, ECS = {:.2f}'.format(
        reg_stats.rvalue, -reg_stats.slope, reg_stats.intercept, ecs)
    plot_path = get_plot_filename(dataset_name, cfg)
    plot.scatterplot(
        [tas_cube.data, x_reg],
        [rtnt_cube.data, y_reg],
        plot_path,
        plot_kwargs=[{
            'linestyle': 'none',
            'markeredgecolor': 'b',
            'markerfacecolor': 'none',
            'marker': 's',
        }, {
            'color': 'k',
            'linestyle': '-',
        }],
        save_kwargs={
            'bbox_inches': 'tight',
            'orientation': 'landscape',
        },
        axes_functions={
            'set_title': dataset_name,
            'set_xlabel': 'tas / ' + tas_cube.units.origin,
            'set_ylabel': 'rtnt / ' + rtnt_cube.units.origin,
            'set_xlim': [0.0, 8.0],
            'set_ylim': [-2.0, 10.0],
            'text': {
                'args': [0.05, 0.9, text],
                'kwargs': {
                    'transform': 'transAxes'
                },
            },
        },
    )

    # Write netcdf file for every plot
    tas_coord = iris.coords.AuxCoord(
        tas_cube.data,
        **extract_variables(cfg, as_iris=True)['tas'])
    attrs = {
        'model': dataset_name,
        'regression_r_value': reg_stats.rvalue,
        'regression_slope': reg_stats.slope,
        'regression_interception': reg_stats.intercept,
        'Climate Feedback Parameter': -reg_stats.slope,
        'ECS': ecs,
    }
    attrs.update(cfg.get('output_attributes', {}))
    cube = iris.cube.Cube(rtnt_cube.data,
                          attributes=attrs,
                          aux_coords_and_dims=[(tas_coord, 0)],
                          **extract_variables(cfg, as_iris=True)['rtnt'])
    netcdf_path = get_diagnostic_filename('ecs_regression_' + dataset_name,
                                          cfg)
    io.iris_save(cube, netcdf_path)

    # Provenance
    provenance_record = get_provenance_record(
        f"Scatterplot between TOA radiance and global mean surface "
        f"temperature anomaly for 150 years of the abrupt 4x CO2 experiment "
        f"including linear regression to calculate ECS for {dataset_name}.")
    provenance_record.update({
        'plot_file': plot_path,
        'plot_types': ['scatter'],
    })

    return (netcdf_path, provenance_record)


def write_data(ecs_data, feedback_parameter_data, ancestor_files, cfg):
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
        attrs['net_toa_radiation'] = (
            f"For datasets {RTMT_DATASETS}, 'rtmt' (net top of model "
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
            raise ValueError(f"No 'rtnt' data for '{dataset_name}' available")
        tas_cube = tas_data[dataset_name][0]['cube']
        rtnt_cube = rtnt_data[dataset_name][0]['cube']
        ancestor_files = (tas_data[dataset_name][0]['ancestors'] +
                          rtnt_data[dataset_name][0]['ancestors'])

        # Perform linear regression
        reg = stats.linregress(tas_cube.data, rtnt_cube.data)

        # Plot ECS regression if desired
        (path,
         provenance_record) = plot_ecs_regression(cfg, dataset_name, tas_cube,
                                                  rtnt_cube, reg)

        # Provenance
        if path is not None:
            provenance_record['ancestors'] = ancestor_files
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(path, provenance_record)

        # Save data
        if cfg.get('read_external_file') and dataset_name in ecs:
            logger.warning(
                "Overwriting externally given ECS and climate feedback "
                "parameter from file '%s' for '%s'", external_file,
                dataset_name)
        ecs[dataset_name] = -reg.intercept / (2 * reg.slope)
        feedback_parameter[dataset_name] = -reg.slope
        all_ancestors.extend(ancestor_files)

    # Write data
    if external_file is not None:
        all_ancestors.append(external_file)
    all_ancestors = list(set(all_ancestors))
    write_data(ecs, feedback_parameter, all_ancestors, cfg)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
