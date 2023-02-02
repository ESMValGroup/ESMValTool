#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create arbitrary Seaborn plots.

Description
-----------
This diagnostic provides a high-level interface to Seaborn. For this, the input
data is arranged into a single :class:`pandas.DataFrame`, which is then used as
input for the Seaborn function defined by the option `seaborn_func`.

Caveats
-------
All datasets of a given variable must have the same units (e.g., it is not
allowed to use datasets with units `K` and datasets with units `°C` for the
variable `tas`).

Author
------
Manuel Schlund (DLR, Germany)

Configuration options in recipe
-------------------------------
add_ancillary_variables: bool, optional (default: False)
    Add :meth:`~iris.cube.Cube.ancillary_variables` to the main data frame.
    Note that this will assume that ancillary variables are identical across
    cubes within a group (see option `groupby_facet`). This equality is not
    checked!
add_aux_coords: bool, optional (default: False)
    Add :meth:`~iris.cube.Cube.aux_coords` to the main data frame. Note that
    this will assume that auxiliary coordinates are identical across cubes
    within a group (see option `groupby_facet`). This equality is not checked!
add_cell_measures: bool, optional (default: False)
    Add :meth:`~iris.cube.Cube.cell_measures` to the main data frame. Note that
    this will assume that cell measures are identical across cubes within a
    group (see option `groupby_facet`). This equality is not checked!
data_frame_ops: dict, optional
    Perform additional operations on the main data frame. Allowed operations
    are :func:`pandas.DataFrame.query` (dict key `query`) and
    :func:`pandas.DataFrame.eval` (dict key `eval`). Operations are defined by
    strings (dict values). Examples: ``{'query': 'latitude > 80', 'eval':
    'longitude = longitude - 180.0'}``.
facets_as_columns: list of str, optional
    Facets that will will be added as a columns to the main data frame. Values
    for these facets must be identical across all datasets within a group (see
    option `groupby_facet`).
groupby_facet: str, optional (default: 'alias')
    Facet which is used to group input datasets when creating the main data
    frame. All datasets within a group are expected to have the same index
    after calling :func:`iris.pandas.as_data_frame` on them. These datasets
    within a group will then get merged into a single data frame per group.
    Finally, the data frames for all groups are concatenated into one main data
    frame. `groupby_facet` is also added as a column to this main data frame.
reset_index: bool, optional (default: False)
    Put coordinate information of datasets into columns instead of (multi-)
    indices. This avoids the deletion of coordinate information if different
    groups of datasets have different dimensions but increases the memory
    footprint of this diagnostic.
savefig_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.savefig`. By
    default, uses ``bbox_inches: tight, dpi: 300, orientation: landscape``.
seaborn_func: str
    Function used to plot the data. Must be a function of
    :mod:`seaborn`. An overview of seaborn's plotting functions is given `here
    <https://seaborn.pydata.org/tutorial/function_overview.html>`__.
seaborn_kwargs: dict, optional
    Optional keyword arguments for the plotting function given by
    `seaborn_func`. Must not be include an argument called `data`.
seaborn_settings: dict, optional
    Options for :func:`seaborn.set_theme` (affects all plots).
set_kwargs: dict, optional
    Optional keyword arguments for the ``set`` method of the object returned by
    the plotting function (`seaborn_func`). This will either be
    :meth:`matplotlib.axes.Axes.set` (e.g., :func:`~seaborn.scatterplot`,
    :func:`~seaborn.lineplot`), :meth:`seaborn.FacetGrid.set` (e.g.,
    :func:`~seaborn.relplot`, :func:`~seaborn.displot`),
    :meth:`seaborn.JointGrid.set` (e.g., :func:`~seaborn.jointplot`), or
    :meth:`seaborn.PairGrid.set` (e.g., :func:`~seaborn.pairplot`).
suptitle: str or None, optional (default: None)
    Suptitle for the plot (see :func:`matplotlib.pyplot.suptitle`). If
    ``None``, do not create a suptitle.

"""
from __future__ import annotations

import logging
from copy import deepcopy
from pathlib import Path
from pprint import pformat

import iris
import iris.pandas
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(Path(__file__).stem)

# Use the new behavior of :func:`iris.pandas.as_data_frame`
iris.FUTURE.pandas_ndim = True

# Save units of different variables
# Note: units must be unique across datasets of the same variable
UNITS: dict[str, str] = {}


def _create_plot(
    plot_func: callable,
    data_frame: pd.DataFrame,
    cfg: dict,
) -> None:
    """Create plot."""
    logger.debug(
        "Using main data frame as input for plotting:\n%s", data_frame
    )

    # Plot
    plot_kwargs = cfg['seaborn_kwargs']
    plot_func_str = cfg['seaborn_func']
    if 'data' in plot_kwargs:
        raise ValueError("'data' is an invalid argument for 'seaborn_kwargs'")
    logger.info(
        "Creating plot with\nseaborn.%s(\n    data=main_data_frame,\n%s\n)",
        plot_func_str,
        "\n".join(f"    {k}={v!r}," for (k, v) in plot_kwargs.items()),
    )
    plot_obj = plot_func(data=data_frame, **plot_kwargs)

    # Adjust plot appearance
    plot_obj.set(**cfg['set_kwargs'])
    if cfg['suptitle'] is not None:
        plt.suptitle(cfg['suptitle'])

    # Save plot
    plot_path = get_plot_filename(f"seaborn_{plot_func_str}", cfg)
    plt.savefig(plot_path, **cfg['savefig_kwargs'])
    logger.info("Wrote %s", plot_path)
    plt.close()

    # Provenance tracking
    caption = f"Seaborn {cfg['seaborn_func']} for one or more dataset(s)"
    ancestors = [d['filename'] for d in cfg['input_data'].values()]
    provenance_record = {
        'ancestors': ancestors,
        'authors': ['schlund_manuel'],
        'caption': caption,
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(plot_path, provenance_record)


def _get_grouped_data(cfg: dict) -> dict:
    """Get grouped input data."""
    groupby_facet = cfg['groupby_facet']
    input_data = list(cfg['input_data'].values())

    # Check if necessary facets are present
    for dataset in input_data:
        if groupby_facet not in dataset:
            raise ValueError(
                f"groupby_facet '{groupby_facet}' is not available for "
                f"dataset {dataset['filename']}"
            )
        for facet in cfg['facets_as_columns']:
            if facet not in dataset:
                raise ValueError(
                    f"Facet '{facet}' used for option 'facets_as_columns' is "
                    f"not available for dataset {dataset['filename']}"
                )

    # Group data accordingly
    grouped_data = group_metadata(
        input_data,
        groupby_facet,
        sort='filename',
    )

    return grouped_data


def _get_dataframe(cfg: dict) -> pd.DataFrame:
    """Get main :class:`pandas.DataFrame` used as input for plotting.

    Note
    ----
    Data is stored in long form, see also :func:`iris.pandas.as_data_frame`.

    """
    logger.info(
        "Grouping datasets by '%s' to create main data frame (data frames "
        "are merged within groups, then concatenated across groups)",
        cfg['groupby_facet'],
    )
    if cfg['add_aux_coords']:
        logger.info("Adding aux_coords as columns")
    if cfg['add_cell_measures']:
        logger.info("Adding cell_measures as columns")
    if cfg['add_ancillary_variables']:
        logger.info("Adding ancillary_variables as columns")
    if cfg['facets_as_columns']:
        logger.info("Adding facets as columns: %s", cfg['facets_as_columns'])

    grouped_data = _get_grouped_data(cfg)

    # Merge data frames within groups
    df_dict = {}
    for (group, datasets) in grouped_data.items():
        logger.info("Processing group '%s'", group)
        df_group = _get_df_for_group(cfg, group, datasets)
        df_dict[group] = df_group

    # Concatenate data frames across groups and use dtype 'category' for facet
    # columns to reduce memory usage and decrease computation times
    groupby_facet = cfg['groupby_facet']
    df_main = pd.concat(df_dict.values(), ignore_index=cfg['reset_index'])
    df_main = df_main.astype({
        f: 'category' for f in cfg['facets_as_columns'] + [groupby_facet]
    })

    logger.info("Successfully retrieved main data frame from input data")
    logger.debug("Got main data frame:\n%s", df_main)
    return df_main


def _get_df_for_group(
    cfg: dict,
    group: str,
    datasets: list[dict],
) -> pd.DataFrame:
    """Extract :class:`pandas.DataFrame` for a single group of datasets.

    This merges all data frames of individual datasets of a group.

    """
    df_group = pd.DataFrame()
    facets_as_columns: dict[str, str] = {}
    for dataset in datasets:
        filename = dataset['filename']
        logger.info("Reading %s", filename)
        cube = iris.load_cube(filename)

        # Update units
        short_name = dataset['short_name']
        units = dataset['units']
        if short_name in UNITS and UNITS[short_name] != units:
            raise ValueError(
                f"Got duplicate units for variable '{short_name}': '{units}' "
                f"and '{UNITS[short_name]}'"
            )
        UNITS.setdefault(short_name, units)

        # Get data frame for individual dataset with proper name
        df_dataset = iris.pandas.as_data_frame(
            cube,
            add_aux_coords=cfg['add_aux_coords'],
            add_cell_measures=cfg['add_cell_measures'],
            add_ancillary_variables=cfg['add_ancillary_variables'],
        )
        df_dataset = df_dataset.rename(
            {cube.name(): short_name}, axis='columns'
        )

        # Merge
        if df_group.empty:
            df_group = df_dataset
            facets_as_columns = {
                f: dataset[f] for f in cfg['facets_as_columns']
            }
        else:
            # Make sure that dimensional coordinates match across cubes within
            # a group
            if not df_group.index.equals(df_dataset.index):
                raise ValueError(
                    f"Dimensions of cube {filename} differ from other cubes "
                    f"of group '{group}'. Cubes of that group:\n"
                    f"{pformat([d['filename'] for d in datasets])}"
                )

            # Make sure that facet values used as columns match across datasets
            # within a cube
            for (facet, val) in facets_as_columns.items():
                if dataset[facet] != val:
                    raise ValueError(
                        f"Facet value for facet '{facet}' (used by option "
                        f"'facets_as_columns') of dataset {filename} differs "
                        f"from value of other datasets of group '{group}': "
                        f"expected '{val}', got '{dataset[facet]}'. Datasets "
                        f"of that group:\n"
                        f"{pformat([d['filename'] for d in datasets])}"
                    )
            df_group = pd.merge(
                df_group,
                df_dataset,
                left_index=True,
                right_index=True,
                sort=False,
                suffixes=[None, '_DUPLICATE'],
            )

            # Assume that aux_coords, cell_measures, and ancillary_variables
            # (if requested) are equal across cubes within the group. Only add
            # them when they first appear.
            df_group = df_group.filter(regex='^(?!.*_DUPLICATE)')

    # Move dimensional coordinates from (multi-) index into columns if
    # requested
    if cfg['reset_index']:
        df_group = df_group.reset_index()

    # Add additional information as column and save the data frame
    for (facet, val) in facets_as_columns.items():
        df_group[facet] = val
    if cfg['groupby_facet'] not in df_group.columns:
        df_group[cfg['groupby_facet']] = group

    return df_group


def _get_default_cfg(cfg: dict) -> dict:
    """Get default options for configuration dictionary."""
    cfg = deepcopy(cfg)

    cfg.setdefault('add_ancillary_variables', False)
    cfg.setdefault('add_aux_coords', False)
    cfg.setdefault('add_cell_measures', False)
    cfg.setdefault('data_frame_ops', {})
    cfg.setdefault('facets_as_columns', [])
    cfg.setdefault('groupby_facet', 'alias')
    cfg.setdefault('reset_index', False)
    cfg.setdefault('savefig_kwargs', {
        'bbox_inches': 'tight',
        'dpi': 300,
        'orientation': 'landscape',
    })
    cfg.setdefault('seaborn_kwargs', {})
    cfg.setdefault('seaborn_settings', {})
    cfg.setdefault('set_kwargs', {})
    cfg.setdefault('suptitle', None)

    return cfg


def _get_plot_func(cfg):
    """Get seaborn plot function."""
    if 'seaborn_func' not in cfg:
        raise ValueError("Necessary option 'seaborn_func' missing")
    if not hasattr(sns, cfg['seaborn_func']):
        raise AttributeError(
            f"Invalid seaborn_func '{cfg['seaborn_func']}' (must be a "
            f"function of the module seaborn; an overview of seaborn plotting "
            f"functions is given here: https://seaborn.pydata.org/tutorial/"
            f"function_overview.html)"
        )
    logger.info("Using plotting function seaborn.%s", cfg['seaborn_func'])
    return getattr(sns, cfg['seaborn_func'])


def _modify_dataframe(data_frame: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Modify data frame according to the option ``data_frame_ops``."""
    allowed_funcs = ('query', 'eval')
    for (func, expr) in cfg['data_frame_ops'].items():
        if func not in allowed_funcs:
            raise ValueError(
                f"Got invalid operation '{func}' for option 'data_frame_ops', "
                f"expected one of {allowed_funcs}"
            )
        op_str = f"'{func}' with argument '{expr}'"
        logger.info("Modifying main data frame through operation %s", op_str)
        data_frame = getattr(data_frame, func)(expr)
        logger.debug(
            "Main data frame after operation %s:\n%s", op_str, data_frame
        )
    return data_frame


def main(cfg: dict) -> None:
    """Run diagnostic."""
    cfg = _get_default_cfg(cfg)

    sns.set_theme(**cfg['seaborn_settings'])
    plot_func = _get_plot_func(cfg)

    df_main = _get_dataframe(cfg)
    df_main = _modify_dataframe(df_main, cfg)

    _create_plot(plot_func, df_main, cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
