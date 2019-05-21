#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to calculate ECS following Andrews et al. (2012).

Description
-----------
Calculate the effective climate sensitivity (ECS) using the regression method
proposed by Andrews et al. (2012).

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recipe
-------------------------------
read_external_file : str, optional
    Read ECS from external file.

"""

import logging
import os
from pprint import pformat

import cf_units
import iris
import numpy as np
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


def check_input_data(cfg):
    """Check input data."""
    if not variables_available(cfg, ['tas', 'rtnt']):
        raise ValueError("This diagnostic needs 'tas' and 'rtnt' "
                         "variables if 'read_external_file' is not given")
    input_data = cfg['input_data'].values()
    project_group = group_metadata(input_data, 'project')
    projects = list(project_group.keys())
    if len(projects) > 1:
        raise ValueError("This diagnostic supports only unique 'project' "
                         "attributes, got {}".format(projects))
    project = projects[0]
    if project not in EXP_4XCO2:
        raise ValueError("Project '{}' not supported yet".format(project))
    exp_group = group_metadata(input_data, 'exp')
    exps = set(exp_group.keys())
    if exps != {'piControl', EXP_4XCO2[project]}:
        raise ValueError("This diagnostic needs 'piControl' and '{}' "
                         "experiments, got {}".format(EXP_4XCO2[project],
                                                      exps))


def get_anomaly_data(tas_data, rtnt_data, dataset):
    """Calculate anomaly data for both variables."""
    project = tas_data[0]['project']
    exp_4xco2 = EXP_4XCO2[project]
    paths = {
        'tas_4x': select_metadata(tas_data, dataset=dataset, exp=exp_4xco2),
        'tas_pi': select_metadata(tas_data, dataset=dataset, exp='piControl'),
        'rtnt_4x': select_metadata(rtnt_data, dataset=dataset, exp=exp_4xco2),
        'rtnt_pi': select_metadata(
            rtnt_data, dataset=dataset, exp='piControl'),
    }
    ancestor_files = []
    cubes = {}
    for (key, [path]) in paths.items():
        ancestor_files.append(path['filename'])
        cube = iris.load_cube(path['filename'])
        cube = cube.aggregated_by('year', iris.analysis.MEAN)
        cubes[key] = cube

    # Substract linear fit of piControl run from abrupt4xCO2 experiment
    shape = None
    for cube in cubes.values():
        if shape is None:
            shape = cube.shape
        else:
            if cube.shape != shape:
                raise ValueError(
                    "Expected all cubes of dataset '{}' to have identical "
                    "shapes, got {} and {}".format(dataset, shape, cube.shape))
    tas_pi_reg = stats.linregress(cubes['tas_pi'].coord('year').points,
                                  cubes['tas_pi'].data)
    rtnt_pi_reg = stats.linregress(cubes['rtnt_pi'].coord('year').points,
                                   cubes['rtnt_pi'].data)
    cubes['tas_4x'].data -= (
        tas_pi_reg.slope * cubes['tas_pi'].coord('year').points +
        tas_pi_reg.intercept)
    cubes['rtnt_4x'].data -= (
        rtnt_pi_reg.slope * cubes['rtnt_pi'].coord('year').points +
        rtnt_pi_reg.intercept)
    return (cubes['tas_4x'], cubes['rtnt_4x'], ancestor_files)


def get_provenance_record(caption):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption': caption,
        'statistics': ['mean', 'diff'],
        'domains': ['global'],
        'authors': ['schl_ma'],
        'references': ['andrews12grl'],
        'realms': ['atmos'],
        'themes': ['phys'],
    }
    return record


def read_external_file(cfg):
    """Read external file to get ECS."""
    ecs = {}
    feedback_parameter = {}
    if not cfg.get('read_external_file'):
        return (ecs, feedback_parameter)
    base_dir = os.path.dirname(__file__)
    filepath = os.path.join(base_dir, cfg['read_external_file'])
    if os.path.isfile(filepath):
        with open(filepath, 'r') as infile:
            external_data = yaml.safe_load(infile)
    else:
        logger.error("Desired external file %s does not exist", filepath)
        return (ecs, feedback_parameter)
    ecs = external_data.get('ecs', {})
    feedback_parameter = external_data.get('feedback_parameter', {})
    logger.info("External file %s", filepath)
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
    cube = iris.cube.Cube(
        rtnt_cube.data,
        attributes=attrs,
        aux_coords_and_dims=[(tas_coord, 0)],
        **extract_variables(cfg, as_iris=True)['rtnt'])
    netcdf_path = get_diagnostic_filename('ecs_regression_' + dataset_name,
                                          cfg)
    io.iris_save(cube, netcdf_path)

    # Provenance
    provenance_record = get_provenance_record(
        "Scatterplot between TOA radiance and global mean surface temperature "
        "anomaly for 150 years of the abrupt 4x CO2 experiment including "
        "linear regression to calculate ECS for {}.".format(dataset_name))
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
            'long_name': 'Effective Climate Sensitivity (ECS)',
            'units': cf_units.Unit('K'),
        },
        {
            'short_name': 'lambda',
            'long_name': 'Climate Feedback Parameter',
            'units': cf_units.Unit('W m-2 K-1'),
        },
    ]
    for (idx, var_attr) in enumerate(var_attrs):
        path = get_diagnostic_filename(var_attr['short_name'], cfg)
        io.save_scalar_data(data[idx], path, var_attr)
        caption = "{long_name} for multiple climate models.".format(**var_attr)
        provenance_record = get_provenance_record(caption)
        provenance_record['ancestors'] = ancestor_files
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(path, provenance_record)


def main(cfg):
    """Run the diagnostic."""
    input_data = cfg['input_data'].values()

    # Read external file if desired
    if cfg.get('read_external_file'):
        (ecs, feedback_parameter, external_file) = read_external_file(cfg)
    else:
        check_input_data(cfg)
        ecs = {}
        feedback_parameter = {}
        external_file = None

    # Read data
    tas_data = select_metadata(input_data, short_name='tas')
    rtnt_data = select_metadata(input_data, short_name='rtnt')

    # Iterate over all datasets and save ECS and feedback parameter
    for dataset in group_metadata(tas_data, 'dataset'):
        logger.info("Processing %s", dataset)
        (tas_cube, rtnt_cube, ancestor_files) = get_anomaly_data(
            tas_data, rtnt_data, dataset)

        # Perform linear regression
        reg = stats.linregress(tas_cube.data, rtnt_cube.data)

        # Plot ECS regression if desired
        (path, provenance_record) = plot_ecs_regression(
            cfg, dataset, tas_cube, rtnt_cube, reg)

        # Provenance
        if path is not None:
            provenance_record['ancestors'] = ancestor_files
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(path, provenance_record)

        # Save data
        if cfg.get('read_external_file') and dataset in ecs:
            logger.info(
                "Overwriting external given ECS and climate feedback "
                "parameter for %s", dataset)
        ecs[dataset] = -reg.intercept / (2 * reg.slope)
        feedback_parameter[dataset] = -reg.slope

    # Write data
    ancestor_files = [d['filename'] for d in tas_data + rtnt_data]
    if external_file is not None:
        ancestor_files.append(external_file)
    write_data(ecs, feedback_parameter, ancestor_files, cfg)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
