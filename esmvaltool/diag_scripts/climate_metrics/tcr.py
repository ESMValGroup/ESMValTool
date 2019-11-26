#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to calculate Transient Climate Response (TCR).

Description
-----------
Calculate the transient climate response (see e.g. Gregory and Forster, 2008).

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recipe
-------------------------------
plot : bool, optional (default: True)
    Plot temperature vs. time.
read_external_file : str, optional
    Read TCR from external file. The path can be given relative to this
    diagnostic script or as absolute path.
seaborn_settings : dict, optional
    Options for :func:`seaborn.set` (affects all plots), see
    <https://seaborn.pydata.org/generated/seaborn.set.html>.

"""

import logging
import os
from pprint import pformat

import cf_units
import iris
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import yaml
from scipy import stats

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger, get_diagnostic_filename, get_plot_filename,
    group_metadata, io, run_diagnostic, select_metadata, variables_available)

logger = logging.getLogger(os.path.basename(__file__))

START_YEAR_IDX = 60
END_YEAR_IDX = 80


def _get_anomaly_cube(onepct_cube, pi_cube):
    """Get anomaly cube."""
    onepct_cube = onepct_cube.aggregated_by('year', iris.analysis.MEAN)
    pi_cube = pi_cube.aggregated_by('year', iris.analysis.MEAN)

    # Check cube
    if onepct_cube.ndim != 1:
        raise ValueError(
            f"This diagnostics needs 1D cubes, got {onepct_cube.ndim:d}D cube "
            f"for '1pctCO2' experiment")
    if pi_cube.ndim != 1:
        raise ValueError(
            f"This diagnostics needs 1D cubes, got {pi_cube.ndim:d}D cube for "
            f"'piControl' experiment")
    if onepct_cube.shape != pi_cube.shape:
        raise ValueError(
            f"Cube shapes of '1pctCO2' and 'piControl' are not identical, got "
            f"{onepct_cube.shape} and {pi_cube.shape}")
    if onepct_cube.shape[0] < END_YEAR_IDX:
        raise ValueError(
            f"Cubes need at least {END_YEAR_IDX:d} points for TCR "
            f"calculation, got only {onepct_cube.shape[0]:d}")

    # Calculate anomaly
    reg = stats.linregress(pi_cube.coord('year').points, pi_cube.data)
    onepct_cube.data -= (reg.slope * pi_cube.coord('year').points +
                         reg.intercept)

    # Adapt metadata
    onepct_cube.standard_name = None
    onepct_cube.var_name += '_anomaly'
    onepct_cube.long_name += ' (Anomaly)'
    onepct_cube.attributes['anomaly'] = ('relative to linear fit of piControl '
                                         'run')
    onepct_cube.convert_units('K')
    return onepct_cube


def _plot(cfg, cube, dataset_name, tcr):
    """Create scatterplot of temperature anomaly vs. time."""
    if not cfg['write_plots'] or not cfg.get('plot', True):
        return (None, None)
    logger.debug("Plotting temperature anomaly vs. time for '%s'",
                 dataset_name)
    (_, axes) = plt.subplots()

    # Plot data
    x_data = np.arange(cube.shape[0])
    y_data = cube.data
    axes.scatter(x_data, y_data, color='b', marker='o')

    # Plot lines
    line_kwargs = {'color': 'k', 'linewidth': 1.0, 'linestyle': '--'}
    axes.axhline(tcr, **line_kwargs)
    axes.axvline(START_YEAR_IDX, **line_kwargs)
    axes.axvline(END_YEAR_IDX, **line_kwargs)

    # Appearance
    units_str = (cube.units.symbol
                 if cube.units.origin is None else cube.units.origin)
    axes.set_title(dataset_name)
    axes.set_xlabel('Years after experiment start')
    axes.set_ylabel(f'Temperature anomaly / {units_str}')
    axes.set_ylim([x_data[0] - 1, x_data[-1] + 1])
    axes.set_ylim([-1.0, 7.0])
    axes.text(0.0, tcr + 0.1, 'TCR = {:.1f} {}'.format(tcr, units_str))

    # Save cube
    netcdf_path = get_diagnostic_filename(dataset_name, cfg)
    io.iris_save(cube, netcdf_path)

    # Save plot
    plot_path = get_plot_filename(dataset_name, cfg)
    plt.savefig(plot_path, orientation='landscape', bbox_inches='tight')
    logger.info("Wrote %s", plot_path)
    plt.close()

    # Provenance
    provenance_record = get_provenance_record(
        f"Time series of the global mean surface air temperature anomaly "
        f"(relative to the linear fit of the pre-industrial control run) of "
        f"{dataset_name} for the 1% CO2 increase per year experiment. The "
        f"horizontal dashed line indicates the transient climate response "
        f"(TCR) defined as the 20 year average temperature anomaly centered "
        f"at the time of CO2 doubling (vertical dashed lines).")
    provenance_record.update({
        'plot_file': plot_path,
        'plot_types': ['times'],
    })

    return (netcdf_path, provenance_record)


def calculate_tcr(cfg):
    """Calculate transient climate response (TCR)."""
    tcr = {}

    # Get data
    input_data = cfg['input_data'].values()
    onepct_data = select_metadata(input_data, short_name='tas', exp='1pctCO2')

    # Iterate over all datasets
    for dataset in onepct_data:
        pi_data = select_metadata(input_data,
                                  short_name='tas',
                                  exp='piControl',
                                  dataset=dataset['dataset'])
        if not pi_data:
            raise ValueError(f"No 'piControl' data available for dataset "
                             f"'{dataset['dataset']}'")

        onepct_cube = iris.load_cube(dataset['filename'])
        pi_cube = iris.load_cube(pi_data[0]['filename'])

        # Get anomaly cube
        anomaly_cube = _get_anomaly_cube(onepct_cube, pi_cube)

        # Calculate TCR
        tas_2x = anomaly_cube[START_YEAR_IDX:END_YEAR_IDX].collapsed(
            'time', iris.analysis.MEAN).data
        new_tcr = tas_2x
        tcr[dataset['dataset']] = new_tcr
        logger.info("TCR (%s) = %.2f %s", dataset['dataset'], new_tcr,
                    anomaly_cube.units)

        # Plot
        (path, provenance_record) = _plot(cfg, anomaly_cube,
                                          dataset['dataset'], new_tcr)
        if path is not None:
            provenance_record['ancestors'] = [
                dataset['filename'],
                pi_data[0]['filename'],
            ]
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(path, provenance_record)

    return tcr


def check_input_data(cfg):
    """Check input data."""
    if not variables_available(cfg, ['tas']):
        raise ValueError(
            "This diagnostic needs variable 'tas' if 'read_external_file' is "
            "not given")
    input_data = cfg['input_data'].values()
    project_group = group_metadata(input_data, 'project')
    projects = list(project_group.keys())
    if len(projects) > 1:
        raise ValueError(
            f"This diagnostic supports only unique 'project' attributes, got "
            f"{projects}")
    exp_group = group_metadata(input_data, 'exp')
    exps = set(exp_group.keys())
    if exps != {'piControl', '1pctCO2'}:
        raise ValueError(
            f"This diagnostic needs '1pctCO2' and 'piControl' experiment, got "
            f"{exps}")


def get_provenance_record(caption):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption': caption,
        'statistics': ['mean', 'diff'],
        'domains': ['global'],
        'authors': ['schlund_manuel'],
        'references': ['gregory08jgr'],
        'realms': ['atmos'],
        'themes': ['phys'],
    }
    return record


def read_external_file(cfg):
    """Read external file to get TCR."""
    filepath = os.path.expanduser(os.path.expandvars(
        cfg['read_external_file']))
    if not os.path.isabs(filepath):
        filepath = os.path.join(os.path.dirname(__file__), filepath)
    if not os.path.isfile(filepath):
        raise FileNotFoundError(
            f"Desired external file '{filepath}' does not exist")
    with open(filepath, 'r') as infile:
        external_data = yaml.safe_load(infile)
    tcr = external_data.get('tcr', {})
    logger.info("Reading external file '%s'", filepath)
    logger.info("Found TCR (K):")
    logger.info("%s", pformat(tcr))
    return (tcr, filepath)


def write_data(cfg, tcr, external_file=None):
    """Write netcdf files."""
    var_attr = {
        'short_name': 'tcr',
        'long_name': 'Transient Climate Response (TCR)',
        'units': cf_units.Unit('K'),
    }
    path = get_diagnostic_filename(var_attr['short_name'], cfg)
    project = list(cfg['input_data'].values())[0]['project']
    io.save_scalar_data(tcr, path, var_attr, attributes={'project': project})
    caption = "{long_name} for multiple climate models.".format(**var_attr)
    provenance_record = get_provenance_record(caption)
    ancestor_files = []
    for dataset_name in tcr.keys():
        datasets = select_metadata(cfg['input_data'].values(),
                                   dataset=dataset_name)
        ancestor_files.extend([d['filename'] for d in datasets])
    if external_file is not None:
        ancestor_files.append(external_file)
    provenance_record['ancestors'] = ancestor_files
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(path, provenance_record)


def main(cfg):
    """Run the diagnostic."""
    sns.set(**cfg.get('seaborn_settings', {}))

    # Read external file if desired
    if cfg.get('read_external_file'):
        (tcr, external_file) = read_external_file(cfg)
    else:
        check_input_data(cfg)
        tcr = {}
        external_file = None

    # Calculate TCR directly
    new_tcr = calculate_tcr(cfg)
    for dataset_name in new_tcr:
        if dataset_name in tcr:
            logger.warning(
                "Overwriting externally given TCR from file '%s' for '%s'",
                external_file, dataset_name)
    tcr.update(new_tcr)

    # Write TCR
    write_data(cfg, tcr)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
