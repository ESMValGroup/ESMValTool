#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to create simple multi-dataset barplots.

Description
-----------
Create barplot for different datasets of a single variable. This diagnostic
needs preprocessed 1D cubes with single dimension ``dataset``.

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recipe
-------------------------------
label_attribute : str, optional
    Cube attribute which is used as label for different input files.
patterns : list of str, optional
    Patterns to filter list of input data.
seaborn_settings : dict, optional
    Options for :func:`seaborn.set` (affects all plots).
sort_ascending : bool, optional (default: False)
    Sort bars in ascending order.
sort_descending : bool, optional (default: False)
    Sort bars in descending order.
value_labels : bool, optional (default: False)
    Label bars with value of that bar.
y_range : list of float, optional
    Range for the y axis in the plot.

"""

import logging
import os
from pprint import pformat

import iris
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    io,
    iris_helpers,
    run_diagnostic,
)

logger = logging.getLogger(os.path.basename(__file__))


def get_all_data(cfg, input_files):
    """Get all data."""
    metadata = None
    all_data = {}
    all_files = []
    for filename in input_files:
        all_files.append(filename)
        cube = iris.load_cube(filename)
        try:
            cube.coord('dataset')
        except iris.exceptions.CoordinateNotFoundError:
            raise iris.exceptions.CoordinateNotFoundError(
                f"File '{filename}' does not contain necessary coordinate "
                f"'dataset'")
        logger.info("Processing '%s'", filename)

        # Sort coordinate 'dataset'
        [cube] = iris_helpers.intersect_dataset_coordinates([cube])

        # Add to data dictionary
        if cfg.get('label_attribute') in cube.attributes:
            label = cube.attributes[cfg.get('label_attribute')]
        else:
            label = filename
        all_data[label] = (cube.coord('dataset').points, cube.data)

        # Check cube metadata
        new_metadata = {
            'long_name': cube.long_name,
            'units': cube.units,
            'var_name': cube.var_name.upper(),
        }
        if metadata is None:
            metadata = new_metadata
        else:
            if metadata != new_metadata:
                raise ValueError(
                    f"Got differing metadata for the different input files, "
                    f"{metadata} and {new_metadata}")
    return (all_data, all_files, metadata)


def get_provenance_record(caption, ancestor_files, **kwargs):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption': caption,
        'authors': ['schlund_manuel'],
        'references': ['acknow_project'],
        'ancestors': ancestor_files,
    }
    record.update(kwargs)
    return record


def plot_data(cfg, all_data, metadata):
    """Create barplot."""
    logger.debug("Plotting barplot")
    (_, axes) = plt.subplots(figsize=(8, 4))

    # Plot
    all_pos = []
    x_labels = []
    offset = 0.0
    for (label, xy_data) in all_data.items():
        if cfg.get('sort_ascending'):
            sort_idx = np.argsort(xy_data[1])
        elif cfg.get('sort_descending'):
            sort_idx = np.argsort(xy_data[1])[::-1]
        else:
            sort_idx = np.arange(len(xy_data[1]))
        xy_data = (xy_data[0][sort_idx], xy_data[1][sort_idx])
        pos = np.arange(len(xy_data[0])) + offset + 0.5
        axes.bar(pos, xy_data[1], align='center', label=label)
        all_pos.extend(pos)
        x_labels.extend(xy_data[0])
        offset += len(pos) + 1.0

    # Plot appearance
    axes.set_title(metadata['long_name'])
    axes.set_xticks(all_pos)
    axes.set_xticklabels(x_labels, rotation=45.0, ha='right', size=7.0)
    axes.set_ylabel(f"{metadata['var_name']} / {metadata['units']}")
    axes.set_ylim(cfg.get('y_range'))
    if 'label_attribute' in cfg:
        axes.legend(loc='upper right')
    if cfg.get('value_labels'):
        for rect in axes.patches:
            axes.text(rect.get_x() + rect.get_width() / 2.0,
                      rect.get_height(),
                      "{:.1f}".format(rect.get_height()),
                      ha='center',
                      va='bottom',
                      size=5.0)

    # Save plot
    plot_path = get_plot_filename(metadata['var_name'], cfg)
    plt.savefig(plot_path,
                orientation='landscape',
                bbox_inches='tight',
                dpi=300)
    logger.info("Wrote %s", plot_path)
    plt.close()
    return plot_path


def write_data(cfg, all_data, metadata):
    """Write netcdf file."""
    new_data = {}
    for (label, xy_data) in all_data.items():
        for (idx, dataset_name) in enumerate(xy_data[0]):
            key = f'{label}-{dataset_name}'
            value = xy_data[1][idx]
            new_data[key] = value
    netcdf_path = get_diagnostic_filename(metadata['var_name'], cfg)
    var_attrs = metadata.copy()
    var_attrs['short_name'] = var_attrs.pop('var_name')
    io.save_scalar_data(new_data, netcdf_path, var_attrs)
    return netcdf_path


def main(cfg):
    """Run the diagnostic."""
    sns.set(**cfg.get('seaborn_settings', {}))
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

    # Iterate over all files and extract data
    (all_data, all_files, metadata) = get_all_data(cfg, input_files)

    # Create plot and netcdf file
    plot_path = plot_data(cfg, all_data, metadata)
    netcdf_path = write_data(cfg, all_data, metadata)

    # Provenance
    caption = f"{metadata['long_name']} for multiple datasets."
    provenance_record = get_provenance_record(caption, all_files)
    provenance_record.update({
        'plot_file': plot_path,
        'plot_types': ['bar'],
    })
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
