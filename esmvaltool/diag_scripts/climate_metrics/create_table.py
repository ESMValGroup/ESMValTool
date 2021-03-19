#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to create table of scalar data.

Description
-----------
Create CSV table for scalar data (per dataset). All input data must be 1D
arrays with single ``'dataset'`` coordinate.

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recipe
-------------------------------
calculate_mean : bool, optional (default: True)
    Calculate mean over all datasets and add it to table.
calculate_std : bool, optional (default: True)
    Calculate standard deviation over all datasets and add it to table.
exclude_datasets : list of str, optional (default: ['MultiModelMean'])
    Exclude certain datasets when calculating statistics over all datasets and
    for assigning an index.
patterns : list of str, optional
    Patterns to filter list of input data.
round_output : int, optional
    If given, round output to given number of decimals.

"""

import logging
import os
from pprint import pformat

import iris
import numpy as np
import pandas as pd

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    io,
    run_diagnostic,
)

logger = logging.getLogger(os.path.basename(__file__))

EXCLUDE_VAL = 0
PANDAS_PRINT_OPTIONS = ['display.max_rows', None, 'display.max_colwidth', -1]


def _add_numerical_index(data_frame, exclude_datasets):
    """Get numerical index."""
    data_frame.loc[:, 'idx'] = np.arange(len(data_frame.index)) + 1
    for exclude_dataset in exclude_datasets:
        idx_exclude_dataset = (data_frame.index.get_level_values('dataset') ==
                               exclude_dataset)
        data_frame.loc[idx_exclude_dataset, 'idx'] = EXCLUDE_VAL
    idx_exclude_vals = (data_frame['idx'] == EXCLUDE_VAL).to_numpy().nonzero()
    for idx in idx_exclude_vals[0]:
        rows_above_exclude = data_frame.index[idx + 1:]
        row_to_exclude = data_frame.index[idx]
        data_frame.loc[rows_above_exclude, 'idx'] -= 1
        data_frame.loc[row_to_exclude, 'idx'] = EXCLUDE_VAL
    data_frame = data_frame.set_index('idx', append=True)
    return data_frame


def _calculate_statistic(data_frame, stat_func, exclude_datasets):
    """Calculate statistic."""
    projects = data_frame.index.get_level_values('project')
    series_to_append = []
    for project in list(set(projects)):
        sub_data_frame = data_frame.loc[projects == project]
        datasets = sub_data_frame.index.get_level_values('dataset')
        for exclude_dataset in exclude_datasets:
            sub_data_frame = sub_data_frame.loc[datasets != exclude_dataset]
        stat = stat_func(sub_data_frame, axis=0)
        series = pd.Series(
            stat,
            index=data_frame.columns,
            name=(project, f'--{stat_func.__name__.upper()}--', EXCLUDE_VAL),
        )
        series_to_append.append(series)
    for series in series_to_append:
        data_frame = data_frame.append(series)
    data_frame = data_frame.sort_index()
    return data_frame


def calculate_statistics(data_frame, cfg):
    """Calculate statistics over all datasets."""
    exclude_datasets = cfg['exclude_datasets']
    if cfg.get('calculate_mean', True):
        logger.info("Calculating mean over all datasets excluding %s",
                    exclude_datasets)
        data_frame = _calculate_statistic(data_frame, np.mean,
                                          exclude_datasets)
    if cfg.get('calculate_std', True):
        logger.info("Calculating standard deviation over all datasets "
                    "excluding %s", exclude_datasets)
        data_frame = _calculate_statistic(data_frame, np.std, exclude_datasets)
    return data_frame


def check_cube(cube, filename):
    """Check properties of cube."""
    if cube.ndim != 1:
        raise ValueError(
            f"Expected 1D data in file '{filename}', got {cube.ndim:d} cube")
    try:
        cube.coord('dataset')
    except iris.exceptions.CoordinateNotFoundError:
        raise iris.exceptions.CoordinateNotFoundError(
            f"Necessary coordinate 'dataset' not found in file '{filename}'")


def create_data_frame(input_files, exclude_datasets):
    """Create data frame."""
    data_frame = pd.DataFrame()
    for input_file in input_files:
        cube = iris.load_cube(input_file)
        check_cube(cube, input_file)
        project = cube.attributes.get('project', 'unknown project')
        index = pd.MultiIndex.from_product([[project],
                                            cube.coord('dataset').points],
                                           names=['project', 'dataset'])
        series = pd.Series(data=cube.data, index=index)

        # Expand index
        for row in series.index.difference(data_frame.index):
            data_frame = data_frame.append(pd.Series(name=row,
                                                     dtype=cube.dtype))

        # Add new data
        if cube.var_name in data_frame.columns:
            for row in series.index:
                if np.isnan(data_frame.loc[row, cube.var_name]):
                    data_frame.loc[row, cube.var_name] = series.loc[row]
                else:
                    if not np.isclose(data_frame.loc[row, cube.var_name],
                                      series.loc[row]):
                        raise ValueError(
                            f"Got duplicate data for '{cube.var_name}' of "
                            f"'{row}': {series.loc[row]:e} and "
                            f"{data_frame.loc[row, cube.var_name]:e}")
        else:
            data_frame.loc[:, cube.var_name] = series

    # Sort and add numerical index for labels
    data_frame.index.names = ['project', 'dataset']
    data_frame = data_frame.sort_index()
    data_frame = _add_numerical_index(data_frame, exclude_datasets)

    return data_frame


def write_provenance(cfg, filename, data_frame, ancestors):
    """Write provenance information."""
    variables = ', '.join(data_frame.columns)
    projects = ', '.join(list(set(data_frame.index.get_level_values(0))))
    caption = (f"Table including variable(s) {variables} for datasets of "
               f"project(s) {projects}.")
    provenance_record = {
        'caption': caption,
        'authors': ['schlund_manuel'],
        'references': ['acknow_project'],
        'ancestors': ancestors,
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance_record)


def main(cfg):
    """Run the diagnostic."""
    cfg.setdefault('exclude_datasets', ['MultiModelMean'])

    # Get input files
    patterns = cfg.get('patterns')
    if patterns is None:
        input_files = io.get_all_ancestor_files(cfg)
    else:
        input_files = []
        for pattern in patterns:
            input_files.extend(io.get_all_ancestor_files(cfg, pattern=pattern))
    if not input_files:
        raise ValueError("No input files found")
    logger.info("Found input files:\n%s", pformat(input_files))

    # Create data frame
    data_frame = create_data_frame(input_files, cfg['exclude_datasets'])

    # Calculate statistics
    data_frame = calculate_statistics(data_frame, cfg)

    # Round output if desired
    if 'round_output' in cfg:
        data_frame = data_frame.round(decimals=cfg['round_output'])

    # Save file
    basename = '-'.join(data_frame.index.levels[0]) + '_'
    basename += '-'.join(data_frame.columns)
    csv_path = get_diagnostic_filename(basename, cfg).replace('.nc', '.csv')
    data_frame.to_csv(csv_path)
    logger.info("Wrote %s", csv_path)
    with pd.option_context(*PANDAS_PRINT_OPTIONS):
        logger.info("Data:\n%s", data_frame)

    # Provenance
    write_provenance(cfg, csv_path, data_frame, input_files)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
