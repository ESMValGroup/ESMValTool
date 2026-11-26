# Copyright 2026 ESMValTool contributors.
"""Fixed-depth historical hydrographic differences and piControl drift.

Use annual means at common depths and on a common regular grid. The
historical-minus-WOA metric is a climatological reference difference;
only the piControl slope is called model drift. All spatial comparisons
use the intersection of valid model and reference cells.
"""

import logging
import sys
from pathlib import Path

import gsw
import iris
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    run_diagnostic,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cosima_common as cc

LOGGER = logging.getLogger(__name__)
REGIONS = (
    ("Global", -90.0, 90.0),
    ("Southern extratropics", -90.0, -30.0),
    ("Tropics", -30.0, 30.0),
    ("Northern extratropics", 30.0, 90.0),
)
REFERENCE_PERIOD = (1981, 2010)
SPATIAL_NDIM = 3
MIN_TREND_SAMPLES = 2
MIN_COMPARISON_MODELS = 2
NORTH_POLE_LATITUDE = 90.0


def annual_fields(cube):
    """Return annual fields as year, depth, latitude, longitude."""
    data = cc.masked_data(cube)
    depth_name = cc.depth_coord_name(cube)
    if depth_name is None:
        msg = "input has no depth coordinate"
        raise ValueError(msg)
    depth = cube.coord(depth_name).points.astype(float)
    depth_axis = cube.coord_dims(depth_name)[0]
    latitude_axis = cube.coord_dims("latitude")[0]
    longitude_axis = cube.coord_dims("longitude")[-1]
    if cube.coords("time") and cube.coord_dims("time"):
        time_axis = cube.coord_dims("time")[0]
        order = (time_axis, depth_axis, latitude_axis, longitude_axis)
        years = np.array(
            [
                date.year
                for date in cube.coord("time").units.num2date(
                    cube.coord("time").points
                )
            ],
            dtype=float,
        )
    else:
        order = (depth_axis, latitude_axis, longitude_axis)
        years = np.array([0.0])
    if len(set(order)) != len(order) or data.ndim != len(order):
        msg = "expected time/depth/latitude/longitude input"
        raise ValueError(msg)
    data = np.transpose(data, order)
    if len(order) == SPATIAL_NDIM:
        data = data[None, ...]
    if years.size != np.unique(years).size:
        msg = "input must have one annual sample per year"
        raise ValueError(msg)
    return data, depth, years


def compatible(model_cube, reference_cube):
    """Reject silently misaligned depth levels or regridding."""
    for name in ("depth", "latitude", "longitude"):
        model_name = (
            cc.depth_coord_name(model_cube) if name == "depth" else name
        )
        ref_name = (
            cc.depth_coord_name(reference_cube) if name == "depth" else name
        )
        model_points = model_cube.coord(model_name).points
        ref_points = reference_cube.coord(ref_name).points
        if model_points.shape != ref_points.shape or not np.allclose(
            model_points, ref_points, atol=1e-5
        ):
            msg = f"model and reference {name} grids differ"
            raise ValueError(msg)


def woa_potential_temperature(temperature_cube, salinity_cube):
    """Convert WOA in-situ temperature to potential temperature at 0 dbar.

    WOA `t_an` is in-situ temperature even when ESMValTool's WOA CMORizer
    labels the field `thetao`. Use its practical salinity, location, and
    pressure to perform the TEOS-10 conversion before comparison.
    """
    compatible(temperature_cube, salinity_cube)
    temperature_cube = temperature_cube.copy()
    salinity_cube = salinity_cube.copy()
    temperature_cube.convert_units("degC")
    salinity_cube.convert_units("1e-3")
    temperature, depth, _ = annual_fields(temperature_cube)
    salinity, _, _ = annual_fields(salinity_cube)
    if temperature.shape != salinity.shape:
        msg = "WOA temperature and salinity shapes differ"
        raise ValueError(msg)
    latitude = cc.lat_2d(temperature_cube)
    longitude = cc.lon_2d(temperature_cube)
    pressure = gsw.p_from_z(-depth[:, None, None], latitude[None, :, :])
    valid = (
        ~np.ma.getmaskarray(temperature)
        & ~np.ma.getmaskarray(salinity)
        & np.isfinite(np.ma.filled(temperature, np.nan))
        & np.isfinite(np.ma.filled(salinity, np.nan))
    )
    practical = np.ma.filled(salinity, np.nan)
    in_situ = np.ma.filled(temperature, np.nan)
    absolute = gsw.SA_from_SP(
        practical,
        pressure[None, ...],
        longitude[None, None, ...],
        latitude[None, None, ...],
    )
    potential = gsw.pt0_from_t(absolute, in_situ, pressure[None, ...])
    return np.ma.masked_where(~valid | ~np.isfinite(potential), potential)


def check_reference_period(years, period=REFERENCE_PERIOD):
    """Require every annual historical field in the WOA normal interval."""
    expected = np.arange(period[0], period[1] + 1)
    if years.shape != expected.shape or not np.array_equal(years, expected):
        msg = f"historical years must match WOA normal {period}"
        raise ValueError(msg)


def region_weights(cube):
    """Return area weights for the global and latitude-band regions."""
    area = np.ma.filled(cc.cell_area(cube), 0.0)
    latitude = cc.lat_2d(cube)
    weights = []
    for _, south, north in REGIONS:
        inside = (latitude >= south) & (
            latitude < north
            if north < NORTH_POLE_LATITUDE
            else latitude <= north
        )
        weights.append(np.where(inside, area, 0.0))
    return np.asarray(weights)


def _weighted_mean(data, weights, valid):
    weight = np.where(valid, weights, 0.0)
    denominator = weight.sum()
    if denominator <= 0:
        return np.nan
    return float(np.sum(np.where(valid, data, 0.0) * weight) / denominator)


def _reference_metrics(historical, reference, weights):
    """Area-weighted mean difference, RMSE and paired-cell coverage."""
    n_region, n_depth = weights.shape[0], historical.shape[1]
    bias = np.full((n_region, n_depth), np.nan)
    rmse = bias.copy()
    coverage = bias.copy()
    model_clim = np.ma.mean(historical, axis=0)
    ref_clim = np.ma.mean(reference, axis=0)
    for region in range(n_region):
        for level in range(n_depth):
            delta = model_clim[level] - ref_clim[level]
            valid = ~np.ma.getmaskarray(delta) & np.isfinite(
                np.ma.filled(delta, np.nan)
            )
            reference_valid = ~np.ma.getmaskarray(
                ref_clim[level]
            ) & np.isfinite(np.ma.filled(ref_clim[level], np.nan))
            total = np.where(reference_valid, weights[region], 0.0).sum()
            if total <= 0:
                continue
            paired = np.where(valid, weights[region], 0.0).sum()
            coverage[region, level] = paired / total
            values = np.ma.filled(delta, 0.0)
            bias[region, level] = _weighted_mean(
                values, weights[region], valid
            )
            squared = _weighted_mean(values**2, weights[region], valid)
            rmse[region, level] = (
                np.sqrt(squared) if np.isfinite(squared) else np.nan
            )
    return bias, rmse, coverage


def _control_metrics(control, years, weights):
    """Common-mask area means and least-squares slopes per century."""
    n_time, n_depth = control.shape[:2]
    n_region = weights.shape[0]
    means = np.full((n_time, n_region, n_depth), np.nan)
    coverage = np.full((n_region, n_depth), np.nan)
    slope = coverage.copy()
    if n_time < MIN_TREND_SAMPLES:
        msg = "piControl drift needs at least two annual means"
        raise ValueError(msg)
    for region in range(n_region):
        for level in range(n_depth):
            field = control[:, level]
            common = np.all(
                ~np.ma.getmaskarray(field)
                & np.isfinite(np.ma.filled(field, np.nan)),
                axis=0,
            )
            first_valid = ~np.ma.getmaskarray(field[0]) & np.isfinite(
                np.ma.filled(field[0], np.nan)
            )
            total = np.where(first_valid, weights[region], 0.0).sum()
            if total <= 0:
                continue
            coverage[region, level] = (
                np.where(common, weights[region], 0.0).sum() / total
            )
            for time in range(n_time):
                means[time, region, level] = _weighted_mean(
                    np.ma.filled(field[time], 0.0), weights[region], common
                )
            if np.isfinite(means[:, region, level]).all():
                slope[region, level] = (
                    np.polyfit(years - years[0], means[:, region, level], 1)[0]
                    * 100.0
                )
    anomaly = means - means[0:1]
    return means, anomaly, slope, coverage


def _global_mean_metrics(cube):
    """Mean, first-year anomaly and slope of archived CMIP6 `thetaoga`."""
    time = cube.coord("time")
    values = cc.masked_data(cube)
    if values.ndim != 1 or cube.coord_dims(time) != (0,):
        msg = "thetaoga must be an annual one-dimensional series"
        raise ValueError(msg)
    years = np.array(
        [date.year for date in time.units.num2date(time.points)], dtype=float
    )
    if len(years) < MIN_TREND_SAMPLES or len(np.unique(years)) != len(years):
        msg = "thetaoga needs at least two unique annual means"
        raise ValueError(msg)
    if np.ma.getmaskarray(values).any() or not np.isfinite(values).all():
        msg = "thetaoga contains missing annual means"
        raise ValueError(msg)
    anomaly = values - values[0]
    slope = float(np.polyfit(years - years[0], values, 1)[0] * 100.0)
    return years, values, anomaly, slope


def _cube(data, name, units, depth, years=None):
    dims = []
    if years is not None:
        dims.append((years, "year", "1"))
    dims.extend(
        ((np.arange(len(REGIONS)), "region", "1"), (depth, "depth", "m"))
    )
    cube = cc.make_cube(data, dims, name, units, name.replace("_", " "))
    cube.attributes["region_names"] = "; ".join(
        region[0] for region in REGIONS
    )
    return cube


def _plot_bias(dataset, variable, depth, bias, rmse, units, cfg):
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    for region, (label, _, _) in enumerate(REGIONS):
        axes[0].plot(bias[region], depth, marker="o", label=label)
        axes[1].plot(rmse[region], depth, marker="o", label=label)
    axes[0].axvline(0.0, color="black", linewidth=0.7)
    axes[0].set_xlabel(f"Model minus WOA ({units})")
    axes[1].set_xlabel(f"RMSE against WOA ({units})")
    axes[0].set_ylabel("Depth (m)")
    axes[0].invert_yaxis()
    axes[0].legend(fontsize=8)
    fig.suptitle(f"{dataset} {variable}: historical climatology vs WOA18")
    path = get_plot_filename(
        f"hydrography_{dataset}_{variable}_reference", cfg
    )
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    return path


def _plot_multimodel_reference(models, variable, metrics, depth, units, cfg):
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)
    for model in models:
        axes[0].plot(
            metrics[model]["reference_bias"][0], depth, marker="o", label=model
        )
        axes[1].plot(
            metrics[model]["reference_rmse"][0], depth, marker="o", label=model
        )
    axes[2].plot(
        metrics[models[0]]["reference_coverage"][0],
        depth,
        marker="o",
        color="black",
    )
    axes[0].axvline(0.0, color="black", linewidth=0.7)
    axes[0].set(xlabel=f"Global model minus WOA ({units})", ylabel="Depth (m)")
    axes[1].set_xlabel(f"Global spatial RMSE ({units})")
    axes[2].set(xlabel="Shared WOA-valid area fraction", xlim=(0, 1))
    axes[0].invert_yaxis()
    axes[0].legend(fontsize=8)
    fig.suptitle(f"Historical climatology: {variable} model comparison")
    plot_path = get_plot_filename(
        f"hydrography_multimodel_{variable}_reference", cfg
    )
    fig.savefig(plot_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    return plot_path


def _reference_comparison(models, variable, results, cfg):
    """Save metrics on the WOA/model-common mask for several models."""
    depth = results[models[0]]["depth"]
    units = results[models[0]]["units"]
    for model in models[1:]:
        if not np.allclose(results[model]["depth"], depth):
            msg = "models have different comparison depth levels"
            raise ValueError(msg)
        if results[model]["units"] != units:
            msg = "models have different comparison units"
            raise ValueError(msg)
    reference_valid = results[models[0]]["reference_valid"]
    weights = results[models[0]]["weights"]
    common = np.logical_and.reduce(
        [
            ~np.ma.getmaskarray(results[model]["difference"])
            & np.isfinite(np.ma.filled(results[model]["difference"], np.nan))
            for model in models
        ]
    )
    common &= reference_valid
    metrics = {}
    for model in models:
        difference = np.ma.filled(results[model]["difference"], 0.0)
        bias = np.full((len(REGIONS), len(depth)), np.nan)
        rmse = bias.copy()
        coverage = bias.copy()
        for region in range(len(REGIONS)):
            for level in range(len(depth)):
                possible = np.where(
                    reference_valid[level], weights[region], 0.0
                ).sum()
                used = np.where(common[level], weights[region], 0.0).sum()
                if possible <= 0:
                    continue
                coverage[region, level] = used / possible
                if used > 0:
                    bias[region, level] = (
                        np.sum(
                            difference[level]
                            * weights[region]
                            * common[level],
                        )
                        / used
                    )
                    rmse[region, level] = np.sqrt(
                        np.sum(
                            difference[level] ** 2
                            * weights[region]
                            * common[level]
                        )
                        / used
                    )
        metrics[model] = {
            "reference_bias": bias,
            "reference_rmse": rmse,
            "reference_coverage": coverage,
        }
    coords = (
        (np.arange(len(models)), "model_index", "1"),
        (np.arange(len(REGIONS)), "region", "1"),
        (depth, "depth", "m"),
    )
    cubes = iris.cube.CubeList()
    for name in ("reference_bias", "reference_rmse", "reference_coverage"):
        values = np.stack([metrics[model][name] for model in models])
        cube = cc.make_cube(
            values,
            coords,
            name,
            "1" if name.endswith("coverage") else units,
            name.replace("_", " "),
        )
        cube.attributes["model_names"] = "; ".join(models)
        cube.attributes["region_names"] = "; ".join(
            region[0] for region in REGIONS
        )
        cube.attributes["comparison_mask"] = (
            "intersection of WOA and every model at each depth"
        )
        cubes.append(cube)
    path = get_diagnostic_filename(
        f"hydrography_multimodel_{variable}_reference", cfg
    )
    iris.save(cubes, path)

    plot_path = _plot_multimodel_reference(
        models, variable, metrics, depth, units, cfg
    )
    return path, plot_path


def _plot_drift(dataset, variable, depth, years, anomaly, slope, cfg):
    fig, axes = plt.subplots(
        2, 2, figsize=(12, 7), sharex=True, constrained_layout=True
    )
    vmax = np.nanmax(np.abs(anomaly))
    vmax = vmax if np.isfinite(vmax) and vmax > 0 else 1.0
    display_units = "°C" if variable == "thetao" else "psu"
    year_edges = np.r_[
        years[0] - 0.5,
        (years[:-1] + years[1:]) / 2,
        years[-1] + 0.5,
    ]
    depth_edges = np.arange(len(depth) + 1) - 0.5
    for region, (label, _, _) in enumerate(REGIONS):
        ax = axes.flat[region]
        image = ax.pcolormesh(
            year_edges,
            depth_edges,
            anomaly[:, region, :].T,
            shading="flat",
            cmap="RdBu_r",
            vmin=-vmax,
            vmax=vmax,
        )
        ax.set_title(
            f"{label}\n10 m slope: {slope[region, 0]:+.2f} "
            f"{display_units}/century",
            fontsize=10,
        )
        if region >= axes.shape[1]:
            ax.set_xlabel("piControl year")
        ax.set_ylabel("Sampled depth (m)")
        ax.set_yticks(np.arange(len(depth)), [f"{level:g}" for level in depth])
        ax.set_ylim(len(depth) - 0.5, -0.5)
    fig.colorbar(
        image,
        ax=axes.ravel().tolist(),
        label=f"Anomaly ({display_units})",
        shrink=0.85,
        pad=0.03,
    )
    fig.suptitle(f"{dataset} {variable}: unforced control drift")
    path = get_plot_filename(f"hydrography_{dataset}_{variable}_drift", cfg)
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    return path


def _plot_global_drift(dataset, years, anomaly, slope, cfg):
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(years, anomaly, color="tab:red", linewidth=1.8)
    ax.axhline(0.0, color="black", linewidth=0.7)
    ax.set(
        xlabel="piControl year",
        ylabel="Global mean anomaly (°C)",
        title=f"{dataset} whole-ocean thetaoga drift: {slope:+.3f} °C/century",
    )
    path = get_plot_filename(f"hydrography_{dataset}_global_mean_drift", cfg)
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    return path


def _reference_products(
    dataset, variable, historical, reference, metadata, reference_name, cfg
):
    """Compute reference metrics, cubes, plots, and multimodel inputs."""
    cubes = []
    ancestors = []
    plot_paths = []
    model_cube = cc.load_cube(historical["filename"], variable)
    ref_cube = cc.load_cube(reference["filename"], variable)
    expected_version = cfg.get("reference_version")
    if expected_version and str(reference.get("version")) != str(
        expected_version
    ):
        msg = "WOA reference version does not match recipe"
        raise ValueError(msg)
    canonical_units = "degC" if variable == "thetao" else "1e-3"
    model_cube.convert_units(canonical_units)
    compatible(model_cube, ref_cube)
    model, depth, model_years = annual_fields(model_cube)
    if cfg.get("reference_period"):
        check_reference_period(model_years, tuple(cfg["reference_period"]))
    if variable == "thetao":
        salinity_ref = next(
            (
                item
                for item in metadata
                if item["dataset"] == reference_name
                and item["short_name"] == "so"
            ),
            None,
        )
        if salinity_ref is None:
            msg = "WOA salinity is needed to convert in-situ temperature"
            raise ValueError(msg)
        salinity_cube = cc.load_cube(salinity_ref["filename"], "so")
        obs = woa_potential_temperature(ref_cube, salinity_cube)
        ancestors.append(salinity_ref["filename"])
    else:
        ref_cube.convert_units(canonical_units)
        obs, _, _ = annual_fields(ref_cube)
    bias, rmse, coverage = _reference_metrics(
        model, obs, region_weights(model_cube)
    )
    comparison = {
        "depth": depth,
        "units": str(model_cube.units),
        "difference": (np.ma.mean(model, axis=0) - np.ma.mean(obs, axis=0)),
        "reference_valid": (
            ~np.ma.getmaskarray(np.ma.mean(obs, axis=0))
            & np.isfinite(np.ma.filled(np.ma.mean(obs, axis=0), np.nan))
        ),
        "weights": region_weights(ref_cube),
        "ancestors": [
            historical["filename"],
            reference["filename"],
        ]
        + ([salinity_ref["filename"]] if variable == "thetao" else []),
    }
    for name, values in (
        ("reference_bias", bias),
        ("reference_rmse", rmse),
        ("reference_coverage", coverage),
    ):
        units = "1" if name.endswith("coverage") else str(model_cube.units)
        cubes.append(_cube(values, name, units, depth))
    plot_paths.append(
        _plot_bias(
            dataset,
            variable,
            depth,
            bias,
            rmse,
            str(model_cube.units),
            cfg,
        )
    )
    ancestors.extend((historical["filename"], reference["filename"]))
    return depth, cubes, ancestors, plot_paths, comparison


def _control_products(dataset, variable, control, depth, cfg):
    """Compute piControl metrics and compare depth levels with historical data."""
    cubes = []
    ancestors = []
    plot_paths = []
    control_cube = cc.load_cube(control["filename"], variable)
    control_data, control_depth, years = annual_fields(control_cube)
    if depth is not None and not np.allclose(depth, control_depth):
        msg = "historical and control depth levels differ"
        raise ValueError(msg)
    depth = control_depth
    means, anomaly, slope, coverage = _control_metrics(
        control_data, years, region_weights(control_cube)
    )
    units = str(control_cube.units)
    cubes.extend(
        (
            _cube(means, "control_mean", units, depth, years),
            _cube(anomaly, "control_anomaly", units, depth, years),
            _cube(
                slope,
                "control_trend_per_century",
                f"{units}/(100 yr)",
                depth,
            ),
            _cube(coverage, "control_coverage", "1", depth),
        )
    )
    plot_paths.append(
        _plot_drift(
            dataset,
            variable,
            depth,
            years,
            anomaly,
            slope,
            cfg,
        )
    )
    ancestors.append(control["filename"])
    return depth, cubes, ancestors, plot_paths


def _write_model_outputs(dataset, variable, cubes, ancestors, plot_paths, cfg):
    """Save per-model NetCDF, plots, and provenance."""
    path = get_diagnostic_filename(f"hydrography_{dataset}_{variable}", cfg)
    iris.save(iris.cube.CubeList(cubes), path)
    record = {
        "caption": (
            f"{dataset} {variable} fixed-depth hydrographic "
            "reference metrics and piControl drift"
        ),
        "statistics": ["mean", "rmsd", "trend"],
        "domains": ["global"],
        "plot_types": ["vert", "times"],
        "authors": ["beucher_romain"],
        "ancestors": ancestors,
    }
    with ProvenanceLogger(cfg) as provenance:
        provenance.log(path, record)
        for plot_path in plot_paths:
            provenance.log(plot_path, record)
    LOGGER.info("Wrote %s", path)


def _write_global_drift(dataset, global_entry, cfg):
    """Save the archived whole-ocean piControl series and trend."""
    global_cube = cc.load_cube(global_entry["filename"], "thetaoga")
    global_cube.convert_units("degC")
    years, mean, anomaly, slope = _global_mean_metrics(global_cube)
    global_cubes = iris.cube.CubeList(
        [
            cc.make_cube(
                mean,
                [(years, "year", "1")],
                "global_mean",
                "degC",
                "whole-ocean mean potential temperature",
            ),
            cc.make_cube(
                anomaly,
                [(years, "year", "1")],
                "global_anomaly",
                "degC",
                "whole-ocean first-year anomaly",
            ),
            cc.make_cube(
                np.array(slope),
                [],
                "global_trend_per_century",
                "degC/(100 yr)",
                "whole-ocean piControl drift",
            ),
        ]
    )
    global_path = get_diagnostic_filename(
        f"hydrography_{dataset}_global_mean_drift", cfg
    )
    iris.save(global_cubes, global_path)
    plot_path = _plot_global_drift(dataset, years, anomaly, slope, cfg)
    record = {
        "caption": (
            f"{dataset} archived whole-ocean thetaoga piControl mean and drift"
        ),
        "statistics": ["mean", "trend"],
        "domains": ["global"],
        "plot_types": ["times"],
        "authors": ["beucher_romain"],
        "ancestors": [global_entry["filename"]],
    }
    with ProvenanceLogger(cfg) as provenance:
        provenance.log(global_path, record)
        provenance.log(plot_path, record)
    LOGGER.info("Wrote %s", global_path)


def _write_multimodel_comparisons(comparisons, cfg):
    """Save reference metrics using one wet mask shared by all models."""
    for variable, results in comparisons.items():
        if len(results) < MIN_COMPARISON_MODELS:
            continue
        models = sorted(results)
        path, plot_path = _reference_comparison(models, variable, results, cfg)
        record = {
            "caption": (
                f"Multi-model {variable} climatological comparison "
                "against WOA; model names follow model_index order"
            ),
            "statistics": ["mean", "rmsd"],
            "domains": ["global"],
            "plot_types": ["vert"],
            "authors": ["beucher_romain"],
            "ancestors": list(
                dict.fromkeys(
                    filename
                    for model in models
                    for filename in results[model]["ancestors"]
                )
            ),
        }
        with ProvenanceLogger(cfg) as provenance:
            provenance.log(path, record)
            provenance.log(plot_path, record)


def _process_variable(
    dataset, variable, entries, metadata, reference_name, cfg, comparisons
):
    """Combine historical reference and piControl products for one model."""
    historical = entries.get((dataset, variable, "historical"))
    control = entries.get((dataset, variable, "piControl"))
    reference = next(
        (
            item
            for item in metadata
            if item["dataset"] == reference_name
            and item["short_name"] == variable
        ),
        None,
    )
    if not any((historical, control)):
        return
    cubes = []
    ancestors = []
    plot_paths = []
    depth = None
    if historical and reference:
        depth, cubes, ancestors, plot_paths, comparison = _reference_products(
            dataset,
            variable,
            historical,
            reference,
            metadata,
            reference_name,
            cfg,
        )
        comparisons[variable][dataset] = comparison
    elif historical:
        LOGGER.warning(
            "%s %s: no %s reference", dataset, variable, reference_name
        )
    if control:
        depth, control_cubes, control_ancestors, control_plots = (
            _control_products(dataset, variable, control, depth, cfg)
        )
        cubes.extend(control_cubes)
        ancestors.extend(control_ancestors)
        plot_paths.extend(control_plots)
    if cubes:
        _write_model_outputs(
            dataset, variable, cubes, ancestors, plot_paths, cfg
        )


def main(cfg):
    """Run the diagnostic for each model and hydrographic variable."""
    reference_name = cfg.get("reference_dataset", "WOA")
    metadata = list(cfg["input_data"].values())
    entries = {}
    for item in metadata:
        key = (item["dataset"], item["short_name"], item.get("exp"))
        if key in entries:
            msg = f"duplicate preprocessed input: {key}"
            raise ValueError(msg)
        entries[key] = item
    comparisons = {"thetao": {}, "so": {}}
    datasets = sorted(
        {
            item["dataset"]
            for item in metadata
            if item.get("project") == "CMIP6"
        }
    )
    for dataset in datasets:
        for variable in ("thetao", "so"):
            _process_variable(
                dataset,
                variable,
                entries,
                metadata,
                reference_name,
                cfg,
                comparisons,
            )
        global_entry = entries.get((dataset, "thetaoga", "piControl"))
        if global_entry:
            _write_global_drift(dataset, global_entry, cfg)
    _write_multimodel_comparisons(comparisons, cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
