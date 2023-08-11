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
add_mean : bool, optional (default: False)
    If ``True``, add a bar representing the mean for each group identified by
    ``group_by_attribute``. If no grouping  is used, add the mean of all
    values. If ``False``, do not add mean.
colors : str or dict, optional
    If ``group_by_attribute`` is not given, specify a single color for all
    bars.  Otherwise, this needs to be a dict with keys given by the group
    elements and values given as single colors. E.g., ``colors: {CMIP5: blue,
    CMIP6: orange}`` for ``group_by_attribute: project``.
group_by_attribute : str, optional
    Cube attribute which is used to group the input data. E.g.,. use
    ``group_by_attribute: project`` on CMIP5 and CMIP6 data to create two
    groups of data (CMIP5 + CMIP6). Note that the input files need to have the
    global attribute given by ``group_by_attribute`` for this to work.
group_mode : str, optional (default: 'location+color')
    Must be one of ``'location+color', 'color'``. If ``'location+color'``,
    separate different groups by their location in the barplot and their
    colors. If ``'color'``, only use colors to separate the different groups.
order : list of str, optional
    Specify the order of the different groups in the barplot. This is only
    relevant of ``group_by_attribute`` is given.  Elements in the list are
    groups, e.g., use ``order: [CMIP5, CMIP6]`` for ``group_by_attribute:
    project``.
patterns : list of str, optional
    Patterns to filter list of input data.
savefig_kwargs : dict, optional
    Keyword arguments for :func:`matplotlib.pyplot.savefig`.
seaborn_settings : dict, optional
    Options for :func:`seaborn.set_theme` (affects all plots).
sort_ascending : bool, optional (default: False)
    Sort bars in ascending order.
sort_descending : bool, optional (default: False)
    Sort bars in descending order.
subplots_kwargs : dict, optional
    Keyword arguments for :func:`matplotlib.pyplot.subplots`.
value_labels : bool, optional (default: False)
    Label bars with value of that bar.
y_range : list of float, optional
    Range for the Y axis in the plot.

"""

import colorsys
import logging
import os
from collections import OrderedDict
from copy import deepcopy
from pprint import pformat

import iris
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import to_rgb

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    io,
    run_diagnostic,
)

logger = logging.getLogger(os.path.basename(__file__))


def _adjust_lightness(rgb_color, amount=0.6):
    """Adjust lightness of given RGB color."""
    hls = colorsys.rgb_to_hls(*rgb_color)
    new_color = (hls[0], max(0.0, min(1.0, amount * hls[1])), hls[2])
    return colorsys.hls_to_rgb(*new_color)


def _get_mean_name(group):
    """Get name for mean bar of a given group."""
    if group is None:
        return 'Mean'
    return f'{group} Mean'


def _add_mean(cfg, all_data):
    """Add mean to datasets if desired."""
    for (group, vals) in all_data.items():
        datasets = vals[0]
        values = vals[1]
        if cfg.get('add_mean', False):
            logger.debug("Adding mean to group %s", group)
            datasets.append(_get_mean_name(group))
            values.append(np.ma.mean(values))
        all_data[group] = (datasets, values)
    return all_data


def _check_cfg(cfg):
    """Check recipe options."""
    # colors
    colors = cfg.get('colors')
    if 'group_by_attribute' in cfg and 'colors' in cfg:
        if not isinstance(colors, dict):
            raise ValueError(
                f"Expected dict for 'colors' when 'group_by_attribute' is "
                f"given, got {type(colors)}")
    elif 'colors' in cfg:
        if not isinstance(colors, str):
            raise ValueError(
                f"Expected str for 'colors' when 'group_by_attribute' is not "
                f"given, got {type(colors)}")
        cfg['colors'] = {None: colors}

    return cfg


def _get_colors_for_group(cfg, group, group_idx, datasets):
    """Get colors for a single group."""
    colors = cfg.get('colors')
    if colors is None:
        color = f'C{group_idx}'
    else:
        if group not in colors:
            raise ValueError(
                f"Group '{group}' not available in 'colors'")
        color = colors[group]
    new_colors = [color] * len(datasets)
    mean_name = _get_mean_name(group)
    datasets = np.array(datasets)
    if mean_name in datasets:
        mean_idx = np.nonzero(datasets == mean_name)[0][0]
        new_colors[mean_idx] = _adjust_lightness(to_rgb(color))
    return (new_colors, color)


def _sort_and_add_colors(cfg, all_data):
    """Sort data and add colors."""
    all_datasets = []
    all_x_pos = []
    all_y_pos = []
    all_colors = []
    all_handles = []

    # If grouped by location, simply put the groups next to each other
    group_mode = cfg.get('group_mode', 'location+color')
    if group_mode == 'location+color':
        logger.info("Using group_mode = '%s'", group_mode)
        offset = 0.0
        for (idx, (group, vals)) in enumerate(all_data.items()):
            datasets = vals[0]
            y_pos = vals[1]

            # Sort datasets and y values
            if cfg.get('sort_ascending'):
                sort_idx = np.argsort(y_pos)
            elif cfg.get('sort_descending'):
                sort_idx = np.argsort(y_pos)[::-1]
            else:
                sort_idx = np.argsort(datasets)
            datasets = np.array(datasets)[sort_idx]
            y_pos = np.array(y_pos)[sort_idx]
            all_datasets.extend(datasets)
            all_y_pos.extend(y_pos)

            # Derive x position
            x_pos = np.arange(len(datasets)) + offset + 0.5
            all_x_pos.extend(x_pos)
            offset += len(x_pos) + 1.0

            # Derive colors
            (new_colors, color) = _get_colors_for_group(cfg, group, idx,
                                                        datasets)
            all_colors.extend(new_colors)

            # Derive handles for legend
            all_handles.append(mpatches.Patch(
                edgecolor='white', facecolor=color, label=group))

    # If grouped by color, combine all groups
    elif group_mode == 'color':
        logger.info("Using group_mode = '%s'", group_mode)
        for (idx, (group, vals)) in enumerate(all_data.items()):
            datasets = vals[0]
            y_pos = vals[1]
            all_datasets.extend(datasets)
            all_y_pos.extend(y_pos)

            # Derive colors
            (new_colors, color) = _get_colors_for_group(cfg, group, idx,
                                                        datasets)
            all_colors.extend(new_colors)

            # Derive handles for legend
            all_handles.append(mpatches.Patch(
                edgecolor='white', facecolor=color, label=group))

        # Sort arrays
        if cfg.get('sort_ascending'):
            sort_idx = np.argsort(all_y_pos)
        elif cfg.get('sort_descending'):
            sort_idx = np.argsort(all_y_pos)[::-1]
        else:
            sort_idx = np.argsort(all_datasets)
        all_datasets = np.array(all_datasets)[sort_idx]
        all_y_pos = np.array(all_y_pos)[sort_idx]
        all_colors = np.array(all_colors, dtype='object')[sort_idx]

        # Derive x position
        all_x_pos = np.arange(len(all_datasets)) + 0.5

    # Other group_modes are not implemented
    else:
        raise NotImplementedError(
            f"group_mode '{group_mode}' not supported yet, expected one of "
            f"'location+color', 'color'")

    return (all_x_pos, all_y_pos, all_datasets, all_colors, all_handles)


def _get_ordered_dict(cfg, all_data):
    """Get desired order of data."""
    if ('group_by_attribute' not in cfg) or ('order' not in cfg):
        return all_data
    new_dict = []
    order = cfg['order']
    if len(order) != len(set(order)):
        raise ValueError(
            f"Expected unique elements for 'order' option, got {order}")
    logger.info("Using order %s for barplot", order)
    if len(order) != len(all_data):
        raise ValueError(
            f"Expected {len(all_data):d} unique elements for 'order' option "
            f"(number of different labels for the barplot), got "
            f"{len(order):d}")
    for label in order:
        if label not in all_data:
            raise ValueError(
                f"Got invalid label '{label}' in 'order' option, expected one "
                f"of {list(all_data.keys())}")
        new_dict.append((label, all_data[label]))
    return OrderedDict(new_dict)


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
        except iris.exceptions.CoordinateNotFoundError as exc:
            raise iris.exceptions.CoordinateNotFoundError(
                f"File '{filename}' does not contain necessary coordinate "
                f"'dataset'") from exc
        logger.info("Processing '%s'", filename)

        # Get group if desired and collect data in a single dictionary whose
        # keys are groups (group=None if no grouping is desired)
        if 'group_by_attribute' in cfg:
            if cfg['group_by_attribute'] not in cube.attributes:
                raise ValueError(
                    f"Grouping by attribute '{cfg['group_by_attribute']}' "
                    f"failed: file {filename} does not contain global "
                    f"attribute '{cfg['group_by_attribute']}'")
            group = cube.attributes[cfg['group_by_attribute']]
        else:
            group = None
        all_data.setdefault(group, ([], []))  # idx 0: datasets, idx 1: values
        new_datasets = cube.coord('dataset').points
        new_values = cube.data
        all_data[group][0].extend(new_datasets)
        all_data[group][1].extend(new_values)

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

    # Respect order given by user
    all_data = _get_ordered_dict(cfg, all_data)

    # Add mean to datasets
    all_data = _add_mean(cfg, all_data)

    # Sort datasets and add colors
    all_data = _sort_and_add_colors(cfg, all_data)

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
    (_, axes) = plt.subplots(**cfg.get('subplots_kwargs', {}))

    # Plot
    (x_pos, y_pos, datasets, colors, handles) = all_data
    axes.bar(x_pos, y_pos, color=colors, align='center')

    # Plot appearance
    axes.set_title(metadata['long_name'])
    axes.set_xticks(x_pos)
    axes.set_xticklabels(datasets, rotation=45.0, ha='right', size=4.0)
    axes.tick_params(axis='x', which='major', pad=-5.0)
    axes.set_ylabel(f"{metadata['var_name']} / {metadata['units']}")
    axes.set_ylim(cfg.get('y_range'))
    if 'group_by_attribute' in cfg:
        axes.legend(handles=handles, loc='upper right')
    if cfg.get('value_labels'):
        for rect in axes.patches:
            axes.text(rect.get_x() + rect.get_width() / 2.0,
                      rect.get_height() + 0.05,
                      f"{rect.get_height():.2f}",
                      rotation=90.0,
                      ha='center',
                      va='bottom',
                      size=5.0)

    # Save plot
    plot_path = get_plot_filename(metadata['var_name'], cfg)
    plt.savefig(plot_path, **cfg['savefig_kwargs'])
    logger.info("Wrote %s", plot_path)
    plt.close()
    return plot_path


def write_data(cfg, all_data, metadata):
    """Write netcdf file."""
    (_, y_pos, datasets, _, _) = all_data
    new_data = dict(zip(datasets, y_pos))
    netcdf_path = get_diagnostic_filename(metadata['var_name'], cfg)
    var_attrs = metadata.copy()
    var_attrs['short_name'] = var_attrs.pop('var_name')
    io.save_scalar_data(new_data, netcdf_path, var_attrs)
    return netcdf_path


def main(cfg):
    """Run the diagnostic."""
    cfg = deepcopy(cfg)
    cfg.setdefault('savefig_kwargs', {
        'dpi': 300,
        'orientation': 'landscape',
        'bbox_inches': 'tight',
    })
    # Check options
    cfg = _check_cfg(cfg)

    # Get input files
    sns.set_theme(**cfg.get('seaborn_settings', {}))
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
        'plot_types': ['bar'],
    })
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)
        provenance_logger.log(plot_path, provenance_record)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
