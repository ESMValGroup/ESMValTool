#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to calculate ECS following Andrews et al. (2012).

Description
-----------
Calculate the equilibrium climate sensitivity (ECS) using the regression method
proposed by Andrews et al. (2012).

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recipe
-------------------------------
plot_ecs_regression : bool, optional (default: False)
    Plot the linear regression graph.
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
    get_plot_filename, group_metadata, plot, run_diagnostic, save_iris_cube,
    save_scalar_data, select_metadata, variables_available)

logger = logging.getLogger(os.path.basename(__file__))

EXP_4XCO2 = {
    'CMIP5': 'abrupt4xCO2',
    'CMIP6': 'abrupt-4xCO2',
}


def check_input_data(cfg):
    """Check input data."""
    if not variables_available(cfg, ['tas', 'rtmt']):
        raise ValueError("This diagnostic needs 'tas' and 'rtmt' "
                         "variables if 'read_external_file' is not given")
    input_data = cfg['input_data'].values()
    project_group = group_metadata(input_data, 'project')
    projects = list(project_group.keys())
    if len(projects) > 1:
        raise ValueError("This diagnostic supports only unique 'project' "
                         "attributes, got {}".format(projects))
    exp_group = group_metadata(input_data, 'exp')
    exps = set(exp_group.keys())
    if exps != {'piControl', EXP_4XCO2[projects[0]]}:
        raise ValueError("This diagnostic needs 'piControl' and '{}' "
                         "experiments, got {}".format(EXP_4XCO2[projects[0]],
                                                      exps))


def get_anomaly_data(tas_data, rtmt_data, dataset):
    """Calculate anomaly data for both variables."""
    project = tas_data[0]['project']
    exp_4xco2 = EXP_4XCO2[project]
    paths = {
        'tas_4x': select_metadata(tas_data, dataset=dataset, exp=exp_4xco2),
        'tas_pi': select_metadata(tas_data, dataset=dataset, exp='piControl'),
        'rtmt_4x': select_metadata(rtmt_data, dataset=dataset, exp=exp_4xco2),
        'rtmt_pi': select_metadata(
            rtmt_data, dataset=dataset, exp='piControl'),
    }
    ancestor_files = []
    cubes = {}
    for (key, [path]) in paths.items():
        ancestor_files.append(path['filename'])
        cube = iris.load_cube(path['filename'])
        cube = cube.aggregated_by('year', iris.analysis.MEAN)
        cubes[key] = cube

    # Substract piControl run from abrupt4xCO2 experiment
    shape = None
    for cube in cubes.values():
        if shape is None:
            shape = cube.shape
        else:
            if cube.shape != shape:
                raise ValueError(
                    "Expected all cubes of dataset '{}' to have identical "
                    "shapes, got {} and {}".format(dataset, shape, cube.shape))
    cubes['tas_4x'].data -= cubes['tas_pi'].data
    cubes['rtmt_4x'].data -= cubes['rtmt_pi'].data
    return (cubes['tas_4x'], cubes['rtmt_4x'], ancestor_files)


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
    clim_sens = {}
    if not cfg.get('read_external_file'):
        return (ecs, clim_sens)
    base_dir = os.path.dirname(__file__)
    filepath = os.path.join(base_dir, cfg['read_external_file'])
    if os.path.isfile(filepath):
        with open(filepath, 'r') as infile:
            external_data = yaml.safe_load(infile)
    else:
        logger.error("Desired external file %s does not exist", filepath)
        return (ecs, clim_sens)
    ecs = external_data.get('ecs', {})
    clim_sens = external_data.get('climate_sensitivity', {})
    logger.info("External file %s", filepath)
    logger.info("Found ECS (K):")
    logger.info("%s", pformat(ecs))
    logger.info("Found climate sensitivities (W m-2 K-1):")
    logger.info("%s", pformat(clim_sens))
    return (ecs, clim_sens, filepath)


def plot_ecs_regression(cfg, dataset_name, tas_cube, rtmt_cube, reg_stats):
    """Plot linear regression used to calculate ECS."""
    if not (cfg['write_plots'] and cfg.get('plot_ecs_regression')):
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
        [rtmt_cube.data, y_reg],
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
            'set_ylabel': 'rtmt / ' + rtmt_cube.units.origin,
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
        'climate_sensitivity': -reg_stats.slope,
        'ECS': ecs,
    }
    cube = iris.cube.Cube(
        rtmt_cube.data,
        attributes=attrs,
        aux_coords_and_dims=[(tas_coord, 0)],
        **extract_variables(cfg, as_iris=True)['rtmt'])
    netcdf_path = get_diagnostic_filename('ecs_regression_' + dataset_name,
                                          cfg)
    save_iris_cube(cube, netcdf_path)

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


def write_data(ecs_data, clim_sens_data, ancestor_files, cfg):
    """Write netcdf files."""
    data = [ecs_data, clim_sens_data]
    var_attrs = [
        {
            'short_name': 'ecs',
            'long_name': 'Equilibrium Climate Sensitivity (ECS)',
            'units': cf_units.Unit('K'),
        },
        {
            'short_name': 'lambda',
            'long_name': 'Climate Sensitivity',
            'units': cf_units.Unit('W m-2 K-1'),
        },
    ]
    for (idx, var_attr) in enumerate(var_attrs):
        path = get_diagnostic_filename(var_attr['short_name'], cfg)
        save_scalar_data(data[idx], path, var_attr)
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
        (ecs, clim_sens, external_file) = read_external_file(cfg)
    else:
        check_input_data(cfg)
        ecs = {}
        clim_sens = {}
        external_file = None

    # Read data
    tas_data = select_metadata(input_data, short_name='tas')
    rtmt_data = select_metadata(input_data, short_name='rtmt')

    # Iterate over all datasets and save ECS and climate sensitivity
    for dataset in group_metadata(tas_data, 'dataset'):
        logger.info("Processing %s", dataset)
        (tas_cube, rtmt_cube, ancestor_files) = get_anomaly_data(
            tas_data, rtmt_data, dataset)

        # Perform linear regression
        reg = stats.linregress(tas_cube.data, rtmt_cube.data)

        # Plot ECS regression if desired
        (path, provenance_record) = plot_ecs_regression(
            cfg, dataset, tas_cube, rtmt_cube, reg)

        # Provenance
        if path is not None:
            provenance_record['ancestors'] = ancestor_files
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(path, provenance_record)

        # Save data
        if cfg.get('read_external_file') and dataset in ecs:
            logger.info(
                "Overwriting external given ECS and climate "
                "sensitivity for %s", dataset)
        ecs[dataset] = -reg.intercept / (2 * reg.slope)
        clim_sens[dataset] = -reg.slope

    # Write data
    ancestor_files = [d['filename'] for d in tas_data + rtmt_data]
    if external_file is not None:
        ancestor_files.append(external_file)
    write_data(ecs, clim_sens, ancestor_files, cfg)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
