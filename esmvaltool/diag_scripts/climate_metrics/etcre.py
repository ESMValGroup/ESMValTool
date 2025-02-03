#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Calculate effective Transient Climate Response to Emissions (eTCRE).

Description
-----------
This diagnostics calucates the effective Transient Climate Response to
Emissions (eTCRE) and produces relevante plots.

Author
------
Manuel Schlund (DLR, Germany)

Configuration options in recipe
-------------------------------
figure_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.figure`. By
    default, uses ``constrained_layout: true``.
gridline_kwargs: dict, optional
    Optional keyword arguments for grid lines. By default, ``color: lightgrey,
    alpha: 0.5`` are used. Use ``gridline_kwargs: false`` to not show grid
    lines.
groupby_facet: str, optional (default: 'alias')
    Facet used to group datasets. Each group must contain exactly two datasets,
    one for the variable ``tas``, and on for the variable ``fco2antt``.
legend_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.legend`. Use
    ``legend_kwargs: false`` to not show legends.
matplotlib_rc_params: dict, optional
    Optional :class:`matplotlib.RcParams` used to customize matplotlib plots.
    Options given here will be passed to :func:`matplotlib.rc_context` and used
    for all plots produced with this diagnostic.
plot_kwargs: dict, optional
    Optional keyword arguments for :func:`iris.plot.plot`. Dictionary keys are
    elements identified by ``groupby_facet`` or ``default``, e.g.,
    ``CMIP6`` if ``groupby_facet: project`` or ``historical`` if
    ``groupby_facet: exp``. Dictionary values are dictionaries used as
    keyword arguments for :func:`iris.plot.plot`.
pyplot_kwargs: dict, optional
    Optional calls to functions of :mod:`matplotlib.pyplot`. Dictionary keys
    are functions of :mod:`matplotlib.pyplot`. Dictionary values are used as
    argument(s) for these functions (if values are dictionaries, these are
    interpreted as keyword arguments; otherwise a single argument is assumed).
savefig_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.savefig`. By
    default, uses ``bbox_inches: tight, dpi: 300, orientation: landscape``.
seaborn_settings: dict, optional
    Options for :func:`seaborn.set_theme` (affects all plots). By default, uses
    ``style: ticks``.
var_emissions: str, optional (default: 'cumulative_fco2antt')
    Short name of the variable describing the cumulative anthropogenic CO2
    emissions.
var_temperature: str, optional (default: 'tas')
    Short name of the variable describing the temperature change.

"""
import logging
from copy import deepcopy
from pathlib import Path

import iris
import iris.plot
import matplotlib as mpl
import seaborn as sns
from iris.exceptions import ConstraintMismatchError
import matplotlib.pyplot as plt


from esmvaltool.diag_scripts.shared._base import get_plot_filename, group_metadata, select_metadata
import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts.shared import (
    run_diagnostic,
)

logger = logging.getLogger(Path(__file__).stem)


def _get_data_for_group(
    cfg: dict,
    group: str,
    group_data: list[dict],
) -> tuple[dict, dict]:
    logger.info("Processing group '%s'", group)
    if len(group_data) != 2:
        raise ValueError(
                f"Expected exactly 2 datasets for group '{group}', got "
                f"{len(group_data)}"
            )
    data_temperature = select_metadata(
            group_data, short_name=cfg["var_temperature"]
        )
    if len(data_temperature) != 1:
        raise ValueError(
                f"Expected exactly 1 temperature dataset with variable "
                f"'{cfg["var_temperature"]}' for group '{group}', got "
                f"{len(data_temperature)}"
            )
    data_emissions = select_metadata(
            group_data, short_name=cfg["var_emissions"]
        )
    if len(data_emissions) != 1:
        raise ValueError(
                f"Expected exactly 1 temperature dataset with variable "
                f"'{cfg["var_emissions"]}' for group '{group}', got "
                f"{len(data_emissions)}"
            )
    temperature = data_temperature[0]["cube"]
    emissions = data_emissions[0]["cube"]
    if temperature.shape[0] != emissions.shape[0]:
        raise ValueError(
                f"Temperature and emissions data for group '{group}' need to "
                f"have identical dimensions, got {temperature.shape[0]} for "
                f"{data_temperature['filename']} and {emissions.shape[0]} for "
                f"{data_emissions['filename']}"
            )
    return data_temperature[0], data_emissions[0]


def _get_default_cfg(cfg: dict) -> dict:
    """Get default options for configuration dictionary."""
    cfg = deepcopy(cfg)

    cfg.setdefault('figure_kwargs', {'constrained_layout': True})
    cfg.setdefault('gridline_kwargs', {})
    cfg.setdefault('groupby_facet', 'alias')
    cfg.setdefault('legend_kwargs', {})
    cfg.setdefault('matplotlib_rc_params', {})
    cfg.setdefault('plot_kwargs', {})
    cfg.setdefault('pyplot_kwargs', {})
    cfg.setdefault('savefig_kwargs', {
        'bbox_inches': 'tight',
        'dpi': 300,
        'orientation': 'landscape',
    })
    cfg.setdefault("seaborn_settings", {"style": "ticks"})
    cfg.setdefault("var_emissions", "cumulative_fco2antt")
    cfg.setdefault("var_temperature", "tas")

    return cfg


def _get_plot_kwargs(all_plot_kwargs: dict, group: str) -> dict:
    """Get keyword arguments for plot functions."""
    all_plot_kwargs = deepcopy(all_plot_kwargs)

    # First get default kwargs, then overwrite them with group-specific
    # ones
    plot_kwargs = all_plot_kwargs.get('default', {})
    plot_kwargs.update(all_plot_kwargs.get(group, {}))
    plot_kwargs.setdefault("label", group)

    return deepcopy(plot_kwargs)


def _load_and_preprocess_data(cfg: dict) -> list[dict]:
    """Load and preprocess data."""
    input_data = list(cfg['input_data'].values())

    for dataset in input_data:
        filename = dataset['filename']
        logger.info("Loading %s", filename)
        cubes = iris.load(filename)
        if len(cubes) == 1:
            cube = cubes[0]
        else:
            var_name = dataset['short_name']
            try:
                cube = cubes.extract_cube(iris.NameConstraint(
                    var_name=var_name
                ))
            except ConstraintMismatchError as exc:
                var_names = [c.var_name for c in cubes]
                raise ValueError(
                    f"Cannot load data: multiple variables ({var_names}) "
                    f"are available in file {filename}, but not the "
                    f"requested '{var_name}'"
                ) from exc

        # Fix time coordinate if present
        if cube.coords('time', dim_coords=True):
            ih.unify_time_coord(cube)

        if cube.ndim != 1:
            raise ValueError(
                f"Expcected 1D data for '{filename}', got {cube.ndim}D data"
            )

        dataset['cube'] = cube

    return input_data


def _process_pyplot_kwargs(**pyplot_kwargs) -> None:
    """Process functions for :mod:`matplotlib.pyplot`."""
    for (func, arg) in pyplot_kwargs.items():
        if arg is None:
            getattr(plt, func)()
        elif isinstance(arg, dict):
            getattr(plt, func)(**arg)
        else:
            getattr(plt, func)(arg)


def main(cfg: dict) -> None:
    """Run diagnostic."""
    sns.set_theme(**cfg['seaborn_settings'])
    fig = plt.figure(**cfg["figure_kwargs"])
    axes = fig.add_subplot()

    # Load and group data
    input_data = _load_and_preprocess_data(cfg)
    grouped_data = group_metadata(input_data, cfg["groupby_facet"])
    logger.info("Grouped data by '%s'", cfg["groupby_facet"])

    # Process groups
    for group, group_data in grouped_data.items():
        data_temperature, data_emissions = _get_data_for_group(
            cfg, group, group_data
        )
        cube_temperature = data_temperature["cube"]
        cube_emissions = data_emissions["cube"]

        # Plot data from one group
        plot_kwargs = _get_plot_kwargs(cfg["plot_kwargs"], group)
        plot_kwargs["axes"] = axes
        iris.plot.plot(cube_emissions, cube_temperature, **plot_kwargs)

    # Plot appearance
    if cfg["gridline_kwargs"] is not False:
        axes.grid(**cfg["gridline_kwargs"])
    if cfg["legend_kwargs"] is not False:
        axes.legend(**cfg["legend_kwargs"])
    _process_pyplot_kwargs(**cfg["pyplot_kwargs"])

    # Save data
    plot_path = get_plot_filename("etcre", cfg)
    fig.savefig(plot_path, **cfg['savefig_kwargs'])
    logger.info("Wrote %s", plot_path)


if __name__ == '__main__':
    with run_diagnostic() as config:
        config = _get_default_cfg(config)
        with mpl.rc_context(config['matplotlib_rc_params']):
            main(config)
