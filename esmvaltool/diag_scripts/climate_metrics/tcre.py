#!/usr/bin/env python
"""Calculate Transient Climate Response to Cumulative CO2 Emissions (TCRE).

Description
-----------
This diagnostics calculates the Transient Climate Response to Cumulative CO2
Emissions (TCRE) and produces relevant plots.

Author
------
Manuel Schlund (DLR, Germany)

Configuration options in recipe
-------------------------------
calc_tcre_period: list of int, optional (default: [90, 110])
    Period considered to calculate TCRE. This needs to be a sequence of two
    integers, where the first integer defines the start and the second integer
    the end of a slice. TCRE is calculated by averaging over a subarray of the
    temperature change array which is determined by this slice. For example, if
    the input data are annual means, the values here correspond to the years
    (measured from simulation start) over which the temperature change is
    averaged (by default from years 90 to 110).
exp_control: str, optional (default: 'esm-piControl')
    Name of the control experiment.
exp_target: str, optional (default: 'esm-flat10')
    Name of the target experiment.
figure_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.figure`. By
    default, uses ``constrained_layout: true``.
gridline_kwargs: dict, optional
    Optional keyword arguments for grid lines. By default, ``color: lightgrey,
    alpha: 0.5`` are used. Use ``gridline_kwargs: false`` to not show grid
    lines.
groupby_facet: str, optional (default: 'dataset')
    Facet used to group datasets. TCRE is calculated for each group element
    individually.
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
    ``CMIP6`` if ``groupby_facet: project`` or ``CESM2`` if ``groupby_facet:
    dataset``. Dictionary values are dictionaries used as keyword arguments for
    :func:`iris.plot.plot`.
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
    emissions. This must be the name of the variable given in the recipe. Note
    that using the :func:`~esmvalcore.preprocessor.cumulative_sum` preprends
    ``cumulative_`` to this name.
var_temperature: str, optional (default: 'tas')
    Short name of the variable describing the temperature change. This must be
    the name of the variable given in the recipe.

"""

import logging
from collections.abc import Sequence
from copy import deepcopy
from pathlib import Path

import iris
import iris.plot
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from iris.cube import Cube
from iris.exceptions import ConstraintMismatchError
from scipy.stats import linregress

import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts.shared import io, run_diagnostic
from esmvaltool.diag_scripts.shared._base import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    select_metadata,
)

logger = logging.getLogger(Path(__file__).stem)


def _calculate_tcre(
    cfg: dict,
    grouped_anomaly_data: dict[str, tuple[Cube, Cube]],
) -> str:
    """Calculate TCRE."""
    # Get center index at interval for emissions calculations
    # Note: For the default flat10 experiment, emissions reach 1000 Pg at the
    # END of year 99 (starting with year 0). For the default calculation period
    # for TCRE, we use 10 years before and after that, which gives us the
    # default [90, 110]. The following formula basically calculates the center
    # (by default: 99) for the calculation period (by default: [90, 110]).
    idx_start = cfg["calc_tcre_period"][0]
    idx_end = cfg["calc_tcre_period"][1]
    idx_center = int((idx_start + idx_end) / 2.0) - 1

    tcre = {}
    target_units = "K/Eg"
    for group, (cube_e, cube_t) in grouped_anomaly_data.items():
        emissions = cube_e[idx_center].data
        temperature = cube_t[idx_start:idx_end].data.mean()
        tcre_cube = Cube(
            np.array(temperature / emissions, dtype=cube_t.dtype),
            units=cube_t.units / cube_e.units,
        )
        tcre_cube.convert_units(target_units)
        tcre[group] = tcre_cube.data
        logger.info("TCRE of %s: %.2f %s", group, tcre_cube.data, target_units)

    netcdf_path = get_diagnostic_filename("tcre", cfg)
    var_attrs = {
        "short_name": "tcre",
        "long_name": (
            "Transient Climate Response to Cumulative CO2 Emissions (TCRE)"
        ),
        "units": target_units,
    }
    attrs = {
        "comment": (
            "Expressed relative to mass of carbon (not CO2), i.e., units are "
            "K/EgC.",
        ),
    }
    io.save_scalar_data(tcre, netcdf_path, var_attrs, attributes=attrs)

    return netcdf_path


def _get_default_cfg(cfg: dict) -> dict:
    """Get default options for configuration dictionary."""
    cfg = deepcopy(cfg)

    cfg.setdefault("calc_tcre_period", [90, 110])
    cfg.setdefault("exp_control", "esm-piControl")
    cfg.setdefault("exp_target", "esm-flat10")
    cfg.setdefault("figure_kwargs", {"constrained_layout": True})
    cfg.setdefault("gridline_kwargs", {})
    cfg.setdefault("groupby_facet", "dataset")
    cfg.setdefault("legend_kwargs", {})
    cfg.setdefault("matplotlib_rc_params", {})
    cfg.setdefault("plot_kwargs", {})
    cfg.setdefault("pyplot_kwargs", {})
    cfg.setdefault(
        "savefig_kwargs",
        {
            "bbox_inches": "tight",
            "dpi": 300,
            "orientation": "landscape",
        },
    )
    cfg.setdefault("seaborn_settings", {"style": "ticks"})
    cfg.setdefault("var_emissions", "cumulative_fco2antt")
    cfg.setdefault("var_temperature", "tas")

    if not isinstance(cfg["calc_tcre_period"], Sequence):
        raise ValueError(
            f"Option 'calc_tcre_period' needs to be a sequence of exactly 2 "
            f"integers, got '{cfg['calc_tcre_period']}' of type "
            f"{type(cfg['calc_tcre_period'])}"
        )
    if len(cfg["calc_tcre_period"]) != 2:
        raise ValueError(
            f"Option 'calc_tcre_period' needs to be a sequence of exactly 2 "
            f"integers, got '{cfg['calc_tcre_period']}'"
        )
    if any(not isinstance(i, int) for i in cfg["calc_tcre_period"]):
        raise ValueError(
            f"Option 'calc_tcre_period' needs to be a sequence of exactly 2 "
            f"integers, got '{cfg['calc_tcre_period']}'"
        )
    if cfg["calc_tcre_period"][0] >= cfg["calc_tcre_period"][1]:
        raise ValueError(
            f"Invalid value for option 'calc_tcre_period': the first integer "
            f"needs to be smaller than the second, got "
            f"'{cfg['calc_tcre_period']}'"
        )

    return cfg


def _get_grouped_anomaly_data(
    cfg: dict,
    input_data: list[dict],
) -> dict[str, tuple[Cube, Cube]]:
    """Get grouped anomaly data."""
    grouped_anomaly_data: dict[str, tuple[Cube, Cube]] = {}
    grouped_data = group_metadata(input_data, cfg["groupby_facet"])
    logger.info("Grouping input data by '%s'", cfg["groupby_facet"])
    for group, group_data in grouped_data.items():
        var_t = cfg["var_temperature"]
        var_e = cfg["var_emissions"]
        exp_control = cfg["exp_control"]
        exp_target = cfg["exp_target"]
        data_e_target = select_metadata(
            group_data, short_name=var_e, exp=exp_target
        )
        data_t_target = select_metadata(
            group_data, short_name=var_t, exp=exp_target
        )
        data_t_control = select_metadata(
            group_data, short_name=var_t, exp=exp_control
        )
        if len(data_e_target) != 1:
            raise ValueError(
                f"Expected exactly 1 dataset for emissions (variable name "
                f"'{var_e}') of experiment '{exp_target}' for group "
                f"'{group}', got {len(data_e_target)}"
            )
        if len(data_t_target) != 1:
            raise ValueError(
                f"Expected exactly 1 dataset for temperature (variable name "
                f"'{var_t}') of experiment '{exp_target}' for group "
                f"'{group}', got {len(data_t_target)}"
            )
        if len(data_t_control) != 1:
            raise ValueError(
                f"Expected exactly 1 dataset for temperature (variable name "
                f"'{var_t}') of experiment '{exp_control}' for group "
                f"'{group}', got {len(data_t_control)}"
            )
        data_e_target = data_e_target[0]
        data_t_target = data_t_target[0]
        data_t_control = data_t_control[0]

        cube_e_target = data_e_target["cube"]
        cube_t_target = data_t_target["cube"]
        cube_t_control = data_t_control["cube"]
        if cube_e_target.shape != cube_t_target.shape:
            raise ValueError(
                f"Expected identical shapes for all data of group '{group}', "
                f"got {cube_e_target.shape} for {data_e_target['filename']} "
                f"and {cube_t_target.shape} for {data_t_target['filename']}"
            )
        if cube_e_target.shape != cube_t_control.shape:
            raise ValueError(
                f"Expected identical shapes for all data of group '{group}', "
                f"got {cube_e_target.shape} for {data_e_target['filename']} "
                f"and {cube_t_control.shape} for {data_t_control['filename']}"
            )

        # Calculate temperature anomaly of target experiment relative to linear
        # regression of control experiment
        time_points = cube_t_control.coord("time").points
        reg = linregress(time_points, cube_t_control.data)
        ref = (reg.slope * time_points + reg.intercept).astype(
            cube_t_target.dtype
        )
        cube_t_anomaly = cube_t_target.copy(cube_t_target.data - ref)

        grouped_anomaly_data[group] = (cube_e_target, cube_t_anomaly)

    return grouped_anomaly_data


def _get_plot_kwargs(all_plot_kwargs: dict, group: str) -> dict:
    """Get keyword arguments for plot functions."""
    all_plot_kwargs = deepcopy(all_plot_kwargs)

    # First get default kwargs, then overwrite them with group-specific
    # ones
    plot_kwargs = all_plot_kwargs.get("default", {})
    plot_kwargs.update(all_plot_kwargs.get(group, {}))
    plot_kwargs.setdefault("label", group)

    return deepcopy(plot_kwargs)


def _load_and_preprocess_data(cfg: dict) -> list[dict]:
    """Load and preprocess data."""
    input_data = list(cfg["input_data"].values())

    for dataset in input_data:
        filename = dataset["filename"]
        logger.info("Loading %s", filename)
        cubes = iris.load(filename)
        if len(cubes) == 1:
            cube = cubes[0]
        else:
            var_name = dataset["short_name"]
            try:
                cube = cubes.extract_cube(
                    iris.NameConstraint(var_name=var_name)
                )
            except ConstraintMismatchError as exc:
                var_names = [c.var_name for c in cubes]
                raise ValueError(
                    f"Cannot load data: multiple variables ({var_names}) "
                    f"are available in file {filename}, but not the "
                    f"requested '{var_name}'"
                ) from exc

        if cube.ndim != 1:
            raise ValueError(
                f"Expcected 1D data for {filename}, got {cube.ndim}D data"
            )
        if not cube.coords("time", dim_coords=True):
            raise ValueError(
                f"Dimensional coordinate 'time' not found for {filename}"
            )
        if cube.shape[0] < cfg["calc_tcre_period"][1]:
            raise ValueError(
                f"The index given by calc_tcre_period needs to be smaller or "
                f"equal to the array size of {filename}, got "
                f"{cfg['calc_tcre_period'][1]} and {cube.shape[0]}, "
                f"respectively"
            )

        ih.unify_time_coord(cube)

        dataset["cube"] = cube

    return input_data


def _plot(
    cfg: dict,
    grouped_anomaly_data: dict[str, tuple[Cube, Cube]],
) -> str:
    """Plot temperature anomaly vs. CO2 emissions."""
    logger.info("Plotting temperature anomaly vs. CO2 emissions")
    fig = plt.figure(**cfg["figure_kwargs"])
    axes = fig.add_subplot()

    # Plot each group individually
    for group, (cube_e, cube_t) in grouped_anomaly_data.items():
        logger.info("Plotting group '%s'", group)
        plot_kwargs = _get_plot_kwargs(cfg["plot_kwargs"], group)
        plot_kwargs["axes"] = axes
        iris.plot.plot(cube_e, cube_t, **plot_kwargs)

    # Plot appearance
    if cfg["gridline_kwargs"] is not False:
        axes.grid(**cfg["gridline_kwargs"])
    if cfg["legend_kwargs"] is not False:
        axes.legend(**cfg["legend_kwargs"])
    _process_pyplot_kwargs(**cfg["pyplot_kwargs"])

    # Save plot
    plot_path = get_plot_filename("tcre", cfg)
    fig.savefig(plot_path, **cfg["savefig_kwargs"])
    logger.info("Wrote %s", plot_path)

    return plot_path


def _process_pyplot_kwargs(**pyplot_kwargs) -> None:
    """Process functions for :mod:`matplotlib.pyplot`."""
    for func, arg in pyplot_kwargs.items():
        if arg is None:
            getattr(plt, func)()
        elif isinstance(arg, dict):
            getattr(plt, func)(**arg)
        else:
            getattr(plt, func)(arg)


def main(cfg: dict) -> None:
    """Run diagnostic."""
    sns.set_theme(**cfg["seaborn_settings"])

    # Load and group data
    input_data = _load_and_preprocess_data(cfg)
    grouped_anomaly_data = _get_grouped_anomaly_data(cfg, input_data)

    # Plot data
    plot_path = _plot(cfg, grouped_anomaly_data)
    provenance_record = {
        "authors": ["schlund_manuel"],
        "ancestors": [d["filename"] for d in input_data],
        "caption": (
            "Surface air temperature anomaly versus cumulative CO2 emissions "
            "for multiple datasets."
        ),
        "plot_types": ["line"],
        "references": ["sanderson24gmd"],
        "realms": ["atmos"],
        "themes": ["carbon", "bgphys"],
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(plot_path, provenance_record)

    # Calculate TCRE
    netcdf_path = _calculate_tcre(cfg, grouped_anomaly_data)
    provenance_record = {
        "authors": ["schlund_manuel"],
        "ancestors": [d["filename"] for d in input_data],
        "caption": (
            "Transient Climate Response to Cumulative CO2 Emissions (TCRE) "
            "for multiple datasets."
        ),
        "references": ["sanderson24gmd"],
        "realms": ["atmos"],
        "themes": ["carbon", "bgphys"],
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)


if __name__ == "__main__":
    with run_diagnostic() as config:
        config = _get_default_cfg(config)
        with mpl.rc_context(config["matplotlib_rc_params"]):
            main(config)
