#!/usr/bin/env python
"""Calculate time of maximum precipitation.

Description
-----------
This diagnostics calculates the time of maximum precipitation and plots it. The
input data needs to be subdaily and of shape (time, latitude, longitude). The
time dimension should be in local solar time (see
:func:`~esmvalcore.preprocessor.local_solar_time`).

Author
------
Manuel Schlund (DLR, Germany)

Configuration options in recipe
-------------------------------
caption: str, optional
    Figure caption used for provenance tracking. By default, uses "Global map
    of the hour of the daily maximum precipitation for {dataset}.".
cbar_label: str, optional (default: 'Hour of maximum precipitation')
    Colorbar label.
cbar_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.colorbar`. By
    default, uses ``{orientation: 'horizontal'}``.
figure_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.figure`. By
    default, uses ``constrained_layout: true, ticks: [0, 3, 6, 9, 12, 15, 18,
    21]``.
matplotlib_rc_params: dict, optional
    Optional :class:`matplotlib.RcParams` used to customize matplotlib plots.
    Options given here will be passed to :func:`matplotlib.rc_context` and used
    for all plots produced with this diagnostic.
plot_kwargs: dict, optional
    Optional keyword arguments for :func:`iris.plot.pcolormesh`. By default,
    uses ``cmap: twilight``.
projection: str, optional (default: None)
    Projection used for the plot. Needs to be a valid projection class of
    :mod:`cartopy.crs`. Keyword arguments can be specified using the option
    ``projection_kwargs``.
projection_kwargs: dict, optional
    Optional keyword arguments for the projection given by ``projection``. For
    map plots, the default keyword arguments ``{central_longitude: 10}`` are
    used.
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
threshold_factor: float, optional (default: 2.0)
    Threshold factor to mask grid points with low precipitation. Only points
    are shown, where
    max(pr) - min(pr) > threshold_factor / sqrt(N_days / 4 - 1) * max(std(pr))
    (Mooers et al., 2021, https://doi.org/10.1029/2020MS002385).
    At a given grid cell, max(pr)/min(pr) are the maximum/minimum precipitation
    across the entire time range, N_days the number of total days and
    max(std(pr)) the maximum of the standard deviation of precipitation across
    a single day.

"""

import logging
import warnings
from copy import deepcopy
from pathlib import Path
from typing import Any

import cartopy.crs as ccrs
import iris
import iris.plot
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from esmvalcore.preprocessor import climate_statistics
from iris.cube import Cube
from iris.warnings import IrisVagueMetadataWarning

from esmvaltool.diag_scripts.shared import io, run_diagnostic
from esmvaltool.diag_scripts.shared._base import (
    get_diagnostic_filename,
    get_plot_filename,
)

logger = logging.getLogger(Path(__file__).stem)


def _calculate_hour_of_max_precipitation(
    cube: Cube, cfg: dict[str, Any]
) -> Cube:
    """Calculate hour of maximum daily precipitation."""
    diurnal_cycle = climate_statistics(cube, operator="mean", period="hourly")
    time_dim = diurnal_cycle.coord_dims("hour")[0]
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            category=IrisVagueMetadataWarning,
            module="iris",
        )
        hour_of_max_precipitation = diurnal_cycle.collapsed(
            "hour", iris.analysis.MAX
        )
    hour_of_max_precipitation.data = diurnal_cycle.coord("hour").points[
        np.argmax(diurnal_cycle.core_data(), axis=time_dim)
    ]
    return hour_of_max_precipitation


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

    cfg.setdefault(
        "caption",
        "Global map of the hour of the daily maximum precipitation for {dataset}.",
    )
    cfg.setdefault("cbar_label", "Hour of maximum precipitation")
    cfg.setdefault(
        "cbar_kwargs",
        {"orientation": "horizontal", "ticks": [0, 3, 6, 9, 12, 15, 18, 21]},
    )
    cfg.setdefault("figure_kwargs", {"constrained_layout": True})
    cfg.setdefault("matplotlib_rc_params", {})
    cfg.setdefault("plot_kwargs", {"cmap": "twilight"})
    cfg.setdefault("projection", "Robinson")
    cfg.setdefault("projection_kwargs", {"central_longitude": 10})
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
    cfg.setdefault("threshold_factor", 2.0)

    return cfg


def _get_projection(cfg: dict[str, Any]) -> Any:
    """Get plot projection."""
    projection = cfg["projection"]
    projection_kwargs = cfg["projection_kwargs"]

    if not hasattr(ccrs, projection):
        msg = f"Got invalid projection '{projection}', expected class of cartopy.crs"
        raise AttributeError(msg)

    return getattr(ccrs, projection)(**projection_kwargs)


def _plot(cube: Cube, cfg: dict[str, Any], plot_path: Path) -> Cube:
    """Plot map."""
    fig = plt.figure(**cfg["figure_kwargs"])
    axes = fig.add_subplot(projection=_get_projection(cfg))
    plot_kwargs = cfg["plot_kwargs"]
    plot_kwargs["axes"] = axes
    map_plot = iris.plot.pcolormesh(cube, **plot_kwargs)

    _process_pyplot_kwargs(**cfg["pyplot_kwargs"])
    cbar = plt.colorbar(map_plot, ax=axes, **cfg["cbar_kwargs"])
    cbar.set_label(cfg["cbar_label"])

    fig.savefig(plot_path, **cfg["savefig_kwargs"])
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
    cfg = _get_default_cfg(cfg)
    sns.set_theme(**cfg["seaborn_settings"])

    for dataset in cfg["input_data"].values():
        filename = dataset["filename"]
        logger.info("Loading %s", filename)
        cube = iris.load_cube(filename)
        cube = _calculate_hour_of_max_precipitation(cube, cfg)

        logger.info("Plotting map for %s", dataset["alias"])
        plot_path = Path(get_plot_filename(Path(filename).stem, cfg))
        _plot(cube, cfg, plot_path)
        logger.info("Created %s", plot_path)

    # Load and group data
    # input_data = _load_and_preprocess_data(cfg)
    # grouped_anomaly_data = _get_grouped_anomaly_data(cfg, input_data)

    # # Plot data
    # with mpl.rc_context(cfg["matplotlib_rc_params"]):
    #     plot_path = _plot(cfg, grouped_anomaly_data)
    # provenance_record = {
    #     "authors": ["schlund_manuel"],
    #     "ancestors": [d["filename"] for d in input_data],
    #     "plot_types": ["line"],
    #     "references": ["sanderson24gmd"],
    #     "realms": ["atmos"],
    #     "themes": ["carbon", "bgphys"],
    # }
    # provenance_record["caption"] = cfg["caption"]
    # with ProvenanceLogger(cfg) as provenance_logger:
    #     provenance_logger.log(plot_path, provenance_record)

    # # Calculate TCRE
    # netcdf_path = _calculate_tcre(cfg, grouped_anomaly_data)
    # provenance_record = {
    #     "authors": ["schlund_manuel"],
    #     "ancestors": [d["filename"] for d in input_data],
    #     "caption": (
    #         "Transient Climate Response to Cumulative CO2 Emissions (TCRE) "
    #         "for multiple datasets."
    #     ),
    #     "references": ["sanderson24gmd"],
    #     "realms": ["atmos"],
    #     "themes": ["carbon", "bgphys"],
    # }
    # with ProvenanceLogger(cfg) as provenance_logger:
    #     provenance_logger.log(netcdf_path, provenance_record)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
