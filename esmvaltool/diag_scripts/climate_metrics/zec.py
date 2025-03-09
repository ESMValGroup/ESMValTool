#!/usr/bin/env python
"""Diagnostic script to calculate and plot ZEC.

Description
-----------
Calculate and plot ZEC (Zero Emission Commitment) Temperature.
Requires input data from a base simulation (e.g. 1pctCO2, flat-10)
and a dedicated ZEC simulation (e.g. esm-flat10-zec, esm-1pct-brch-1000PgC).

Configuration options in recipe
-------------------------------
zec_year : list, optional (default: [50])
    Calculate ZEC for the 20-year average centered around year x,
    multiple values are possible.
experiments: dict, optional (default: {
    'reference': ['esm-flat10', '1pctCO2'],
    'simulation': ['esm-flat10-zec', 'esm-1pct-brch-1000PgC',
    'esm-1pct-brch-750PgC', 'esm-1pct-brch-2000PgC']
    })
    When using non-default experiments
    to calculate ZEC, the experiment setting is required with values for
    ``reference`` and ``simulation``.
"""

import logging
import os
from copy import deepcopy

import iris
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    io,
    run_diagnostic,
    select_metadata,
)

logger = logging.getLogger(os.path.basename(__file__))


def _get_default_cfg(cfg):
    """Get default options for configuration dictionary."""
    cfg = deepcopy(cfg)

    cfg.setdefault("zec_year", [50])
    cfg.setdefault(
        "experiments",
        {
            "reference": ["esm-flat10", "1pctCO2"],
            "simulation": [
                "esm-flat10-zec",
                "esm-1pct-brch-1000PgC",
                "esm-1pct-brch-750PgC",
                "esm-1pct-brch-2000PgC",
            ],
        },
    )
    return cfg


def grouped_data(cfg):
    """Group input data in preparation of ZEC calculation."""
    input_data = cfg["input_data"].values()
    group_data = group_metadata(input_data, "exp")
    data_grouped = {}
    for exp in group_data:
        if exp in cfg["experiments"]["simulation"]:
            data_grouped["zecmip_data"] = select_metadata(
                input_data, short_name="tas", exp=exp
            )
        elif exp in cfg["experiments"]["reference"]:
            data_grouped["base"] = select_metadata(
                input_data, short_name="tas", exp=exp
            )
        else:
            raise ValueError(
                f"{exp} is not a valid experiment for calculating ZEC, "
                f"please check the configuration value of 'experiments'. "
                f"Current accepted experiments are {cfg['experiments']}"
            )
    if "zecmip_data" not in data_grouped or "base" not in data_grouped:
        raise ValueError(
            f"Data does not include experiments valid for ZEC computation, "
            f"please check the configuration value of 'experiments'. "
            f"Current accepted experiments are {cfg['experiments']}, "
            f"received experiments are {list(group_data.keys())}"
        )
    return data_grouped


def calculate_zec(cfg):
    """Calculate ZEC for each model."""
    zec = {}
    data_grouped = grouped_data(cfg)
    zecmip_data = data_grouped["zecmip_data"]
    base = data_grouped["base"]
    for data in zecmip_data:
        # Account for ensembles by using alias, remove exp name
        name = data["alias"].replace("_" + data["exp"], "")
        tas = iris.load_cube(data["filename"])
        # Match correct anomaly base data, no ensemble key for means
        if "_r" in data["alias"]:
            match_base = select_metadata(
                base, dataset=data["dataset"], ensemble=data["ensemble"]
            )
        else:
            match_base = select_metadata(base, dataset=data["dataset"])
        tas_base = iris.load_cube(match_base[0]["filename"])
        zec_model = deepcopy(tas)
        # Fix time to start at 0
        fix_time(zec_model)
        zec[name] = zec_model - tas_base.data
    save_zec(zec, cfg)
    return zec


def save_zec(zec, cfg):
    """Write ZEC data to file."""
    netcdf_path = get_diagnostic_filename("zec", cfg)
    var_attrs = {
        "short_name": "zec",
        "long_name": "Zero Emissions Commitment (ZEC)",
        "units": "K",
    }
    io.save_1d_data(zec, netcdf_path, "time", var_attrs)
    write_provenance(cfg, netcdf_path)


def fix_time(cube):
    """Fix time to start at year 0."""
    time_coord = iris.coords.DimCoord(
        np.arange(cube.coord("time").shape[0]),
        var_name="time",
        standard_name="time",
        long_name="time",
        units="years",
    )
    cube.remove_coord("time")
    cube.add_dim_coord(time_coord, 0)


def plot_zec_timeseries(zec, cfg):
    """Plot all ZEC timeseries."""
    fig, axes = plt.subplots(figsize=(10, 6))
    for model in zec:
        iris.plot.plot(zec[model], axes=axes, label=model.replace("_", " "))
        axes.axhline(color="lightgrey", linestyle="--")
    axes.set_title("ZEC")
    axes.set_xlabel("Time [yr]")
    axes.set_ylabel("ZEC [K]")
    axes.set_xlim([0, 100])
    axes.legend(
        bbox_to_anchor=(0.50, -0.25),
        loc="center",
        ncol=3,
        handlelength=3.5,
        fontsize=12,
    )

    # Save plot
    plot_path = get_plot_filename("zec_timeseries_all_models", cfg)
    plt.tight_layout()
    fig.savefig(plot_path)
    plt.close()
    logger.info("Wrote %s", plot_path)
    prov_dict = {"caption": "ZEC timeseries", "plot_type": "times"}
    write_provenance(cfg, plot_path, prov_dict)


def calc_zec_x(zec, x_i):
    """Calculate ZEC mean of 20 years centered on year x."""
    zec_x = {}
    for model in zec:
        # ZEC_X is the 20-year anomaly centered at year x
        zec_x[model] = np.mean(zec[model].data[x_i - 10 : x_i + 9])
    return zec_x


def write_zec_x(zec_x, x_i, cfg):
    """Save ZEC_x data with provenance."""
    netcdf_path = get_diagnostic_filename(f"zec_{x_i}", cfg)
    var_attrs = {
        "short_name": "zec",
        "long_name": "Zero Emissions Commitment (ZEC)",
        "units": "K",
    }
    io.save_scalar_data(zec_x, netcdf_path, var_attrs)
    write_provenance(cfg, netcdf_path)


def plot_zec_x_bar(zec_x, x_i, cfg):
    """Plot a barplot of all ZEC_x in ascending order."""
    # Sort by value
    data = dict(sorted(zec_x.items(), key=lambda x: x[1]))
    labels = [x.replace("_", " \n ") for x in data]
    # Make plot
    fig, axes = plt.subplots(figsize=(10, 6))
    axes.bar(labels, list(data.values()))
    axes.axhline(color="lightgrey", linestyle="--")
    axes.set_ylabel(f"ZEC$_{{{x_i}}}$")
    # Save plot
    plot_path = get_plot_filename(f"zec_{x_i}_barplot", cfg)
    plt.tight_layout()
    fig.savefig(plot_path)
    plt.close()
    logger.info("Wrote %s", plot_path)
    prov_dict = {"caption": f"Barplot of ZEC_{x_i}", "plot_type": "bar"}
    write_provenance(cfg, plot_path, prov_dict)


def write_provenance(cfg, path, *args):
    """Write provenance record.

    Optionally pass additional dictionary entries as args.
    """
    input_data = cfg["input_data"].values()
    provenance_record = {
        "authors": ["gier_bettina"],
        "references": ["macdougall20"],
        "ancestors": [d["filename"] for d in input_data],
    }
    for add_dict in args:
        for key in add_dict:
            provenance_record[key] = add_dict[key]
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(path, provenance_record)


def main(cfg):
    """Execute diagnostic."""
    # Calculate ZEC
    zec = calculate_zec(cfg)

    # Plot ZEC Timeseries
    plot_zec_timeseries(zec, cfg)

    # Calculate ZEC at the years given in the recipe with default 50
    x_val = (
        cfg["zec_year"]
        if isinstance(cfg["zec_year"], list)
        else [cfg["zec_year"]]
    )
    for x_i in x_val:
        zec_x = calc_zec_x(zec, x_i)
        write_zec_x(zec_x, x_i, cfg)
        plot_zec_x_bar(zec_x, x_i, cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        config = _get_default_cfg(config)
        main(config)
