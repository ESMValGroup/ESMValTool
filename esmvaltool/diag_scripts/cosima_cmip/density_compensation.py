"""Explain hydrographic density bias as temperature and salinity effects.

The input fields are potential temperature and practical salinity on a
common fixed-depth grid. WOA in-situ temperature is converted by the
hydrographic benchmark helper before these functions are called.
"""

import os
import sys

import gsw
import iris
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    run_diagnostic,
)

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import cosima_common as cc
import hydrographic_benchmark as hb


def density_components(
    model_t, model_s, reference_t, reference_s, depth, latitude, longitude
):
    """Return exact T/S contributions to model-minus-reference density.

    Four matched three-dimensional fields have shape (depth, y, x).
    TEOS-10 calculates in-situ density at the pressure of each depth.
    The two contributions are the symmetric average over the two possible
    replacement orders, so their sum equals the total density difference
    even with a nonlinear equation of state.
    """
    fields = [
        np.ma.asarray(field, dtype=float)
        for field in (model_t, model_s, reference_t, reference_s)
    ]
    if len({field.shape for field in fields}) != 1 or fields[0].ndim != 3:
        raise ValueError("T and S fields must share (depth, y, x) shape")
    depth = np.asarray(depth, dtype=float)
    latitude = np.asarray(latitude, dtype=float)
    longitude = np.asarray(longitude, dtype=float)
    if (
        depth.size != fields[0].shape[0]
        or latitude.shape != fields[0].shape[1:]
        or longitude.shape != latitude.shape
    ):
        raise ValueError(
            "depth and geographic coordinates do not match fields"
        )
    valid = np.logical_and.reduce(
        [
            ~np.ma.getmaskarray(field)
            & np.isfinite(np.ma.filled(field, np.nan))
            for field in fields
        ]
    )
    pressure = gsw.p_from_z(-depth[:, None, None], latitude[None, :, :])

    def density(potential, practical):
        absolute = gsw.SA_from_SP(
            np.ma.filled(practical, np.nan),
            pressure,
            longitude[None, :, :],
            latitude[None, :, :],
        )
        conservative = gsw.CT_from_pt(
            absolute, np.ma.filled(potential, np.nan)
        )
        return gsw.rho(absolute, conservative, pressure)

    model = density(fields[0], fields[1])
    reference = density(fields[2], fields[3])
    warm_only = density(fields[0], fields[3])
    salty_only = density(fields[2], fields[1])
    total = model - reference
    thermal = 0.5 * ((warm_only - reference) + (model - salty_only))
    haline = 0.5 * ((salty_only - reference) + (model - warm_only))
    mask = (
        ~valid
        | ~np.isfinite(total)
        | ~np.isfinite(thermal)
        | ~np.isfinite(haline)
    )
    return tuple(
        np.ma.masked_where(mask, value) for value in (thermal, haline, total)
    )


def density_metrics(thermal, haline, total, reference_valid, weights):
    """Return regional means, RMSE, coverage and cancellation fraction.

    Cancellation is the area-and-effect-weighted fraction of the two
    absolute density contributions that offsets locally. Zero means
    reinforcement; one means perfect cancellation. Cells with tiny
    contributions have correspondingly tiny influence on this ratio.
    """
    n_region, n_depth = weights.shape[0], total.shape[0]
    result = {
        name: np.full((n_region, n_depth), np.nan)
        for name in (
            "thermal",
            "haline",
            "total",
            "rmse",
            "coverage",
            "cancellation",
        )
    }
    for region in range(n_region):
        for level in range(n_depth):
            weight = np.asarray(weights[region], dtype=float)
            paired = ~np.ma.getmaskarray(total[level]) & np.isfinite(
                np.ma.filled(total[level], np.nan)
            )
            ref = np.asarray(reference_valid[level], dtype=bool)
            possible = np.where(ref, weight, 0.0).sum()
            area = np.where(paired, weight, 0.0)
            used = area.sum()
            if possible <= 0:
                continue
            result["coverage"][region, level] = used / possible
            if used <= 0:
                continue
            for name, field in (
                ("thermal", thermal),
                ("haline", haline),
                ("total", total),
            ):
                result[name][region, level] = (
                    np.where(paired, np.ma.filled(field[level], 0.0), 0.0)
                    * area
                ).sum() / used
            values = np.ma.filled(total[level], 0.0)
            result["rmse"][region, level] = np.sqrt(
                (area * values**2).sum() / used
            )
            t = np.ma.filled(thermal[level], 0.0)
            s = np.ma.filled(haline[level], 0.0)
            magnitude = (area * (np.abs(t) + np.abs(s))).sum()
            if magnitude > 0:
                cancelled = (
                    area * (np.abs(t) + np.abs(s) - np.abs(values))
                ).sum()
                result["cancellation"][region, level] = np.clip(
                    cancelled / magnitude, 0.0, 1.0
                )
    return result


def _output_cube(data, name, units, depth):
    dims = ((np.arange(len(hb.REGIONS)), "region", "1"), (depth, "depth", "m"))
    cube = cc.make_cube(data, dims, name, units, name.replace("_", " "))
    cube.attributes["region_names"] = "; ".join(
        region[0] for region in hb.REGIONS
    )
    return cube


def _plot(dataset, depth, metrics, cfg):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    for region, (label, _, _) in enumerate(hb.REGIONS):
        axes[0].plot(metrics["total"][region], depth, marker="o", label=label)
        axes[1].plot(
            metrics["cancellation"][region], depth, marker="o", label=label
        )
    axes[0].axvline(0.0, color="black", linewidth=0.7)
    axes[0].set_xlabel("Model minus WOA density (kg m$^{-3}$)")
    axes[1].set_xlabel("Density cancellation fraction (0–1)")
    axes[1].set_xlim(0, 1)
    axes[0].set_ylabel("Depth (m)")
    axes[0].invert_yaxis()
    axes[0].legend(fontsize=8)
    fig.suptitle(f"{dataset}: hydrographic density bias")
    path = get_plot_filename(f"density_compensation_{dataset}", cfg)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def _multimodel_comparison(models, results, cfg):
    """Save density metrics on cells shared by WOA and every model."""
    depth = results[models[0]]["depth"]
    for model in models[1:]:
        if not np.allclose(results[model]["depth"], depth):
            raise ValueError("models have different comparison depth levels")
    common = np.logical_and.reduce(
        [
            ~np.ma.getmaskarray(results[model]["total"])
            & np.isfinite(np.ma.filled(results[model]["total"], np.nan))
            for model in models
        ]
    )
    reference_valid = results[models[0]]["reference_valid"]
    common &= reference_valid
    metrics_by_model = {}
    for model in models:
        fields = [
            np.ma.masked_where(~common, results[model][name])
            for name in ("thermal", "haline", "total")
        ]
        metrics_by_model[model] = density_metrics(
            *fields, reference_valid, results[models[0]]["weights"]
        )
    coords = (
        (np.arange(len(models)), "model_index", "1"),
        (np.arange(len(hb.REGIONS)), "region", "1"),
        (depth, "depth", "m"),
    )
    cubes = iris.cube.CubeList()
    for name in (
        "thermal",
        "haline",
        "total",
        "rmse",
        "coverage",
        "cancellation",
    ):
        values = np.stack([metrics_by_model[model][name] for model in models])
        cube = cc.make_cube(
            values,
            coords,
            name,
            "1" if name in ("coverage", "cancellation") else "kg m-3",
            name.replace("_", " "),
        )
        cube.attributes["model_names"] = "; ".join(models)
        cube.attributes["region_names"] = "; ".join(
            region[0] for region in hb.REGIONS
        )
        cube.attributes["comparison_mask"] = (
            "intersection of WOA and every model at each depth"
        )
        cubes.append(cube)
    path = get_diagnostic_filename("density_compensation_multimodel", cfg)
    iris.save(cubes, path)

    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)
    for model in models:
        metrics = metrics_by_model[model]
        axes[0].plot(metrics["total"][0], depth, marker="o", label=model)
        axes[1].plot(
            metrics["cancellation"][0], depth, marker="o", label=model
        )
    axes[0].axvline(0.0, color="black", linewidth=0.7)
    axes[0].set(
        xlabel="Global model minus WOA density (kg m$^{-3}$)",
        ylabel="Depth (m)",
    )
    axes[1].set(xlabel="Global local cancellation fraction", xlim=(0, 1))
    axes[2].plot(
        metrics_by_model[models[0]]["coverage"][0],
        depth,
        marker="o",
        color="black",
    )
    axes[2].set(xlabel="Shared WOA-valid area fraction", xlim=(0, 1))
    axes[0].invert_yaxis()
    axes[0].legend(fontsize=8)
    fig.suptitle("Historical climatology: model density comparison")
    plot_path = get_plot_filename("density_compensation_multimodel", cfg)
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path, plot_path


def main(cfg):
    """Compute density attribution for each complete model/WOA quartet."""
    metadata = list(cfg["input_data"].values())
    reference_name = cfg.get("reference_dataset", "WOA")
    references = {
        item["short_name"]: item
        for item in metadata
        if item["dataset"] == reference_name
    }
    if set(references) != {"thetao", "so"}:
        raise ValueError("WOA temperature and salinity are both required")
    if any(
        item.get("version") != cfg.get("reference_version")
        for item in references.values()
    ):
        raise ValueError("WOA version does not match recipe")
    ref_t_cube = cc.load_cube(references["thetao"]["filename"], "thetao")
    ref_s_cube = cc.load_cube(references["so"]["filename"], "so")
    ref_t = hb._woa_potential_temperature(ref_t_cube, ref_s_cube)
    ref_s_cube.convert_units("1e-3")
    ref_s, _, _ = hb._annual_fields(ref_s_cube)
    ref_t = np.ma.mean(ref_t, axis=0)
    ref_s = np.ma.mean(ref_s, axis=0)
    reference_valid = (
        ~np.ma.getmaskarray(ref_t)
        & ~np.ma.getmaskarray(ref_s)
        & np.isfinite(np.ma.filled(ref_t, np.nan))
        & np.isfinite(np.ma.filled(ref_s, np.nan))
    )

    model_inputs = [
        item for item in metadata if item.get("project") == "CMIP6"
    ]
    if any(item.get("exp") != "historical" for item in model_inputs):
        raise ValueError("density comparison requires historical model inputs")
    models = sorted({item["dataset"] for item in model_inputs})
    comparisons = {}
    for dataset in models:
        entries = {
            item["short_name"]: item
            for item in model_inputs
            if item["dataset"] == dataset
        }
        if sum(item["dataset"] == dataset for item in model_inputs) != len(
            entries
        ):
            raise ValueError(f"{dataset} has duplicate model variables")
        if set(entries) != {"thetao", "so"}:
            raise ValueError(f"{dataset} needs historical thetao and so")
        temperature = cc.load_cube(entries["thetao"]["filename"], "thetao")
        salinity = cc.load_cube(entries["so"]["filename"], "so")
        temperature.convert_units("degC")
        salinity.convert_units("1e-3")
        hb._compatible(temperature, salinity)
        hb._compatible(temperature, ref_t_cube)
        model_t, depth, years = hb._annual_fields(temperature)
        model_s, _, salinity_years = hb._annual_fields(salinity)
        hb._check_reference_period(years, tuple(cfg["reference_period"]))
        if not np.array_equal(years, salinity_years):
            raise ValueError("temperature and salinity years differ")
        thermal, haline, total = density_components(
            np.ma.mean(model_t, axis=0),
            np.ma.mean(model_s, axis=0),
            ref_t,
            ref_s,
            depth,
            cc.lat_2d(temperature),
            cc.lon_2d(temperature),
        )
        metrics = density_metrics(
            thermal,
            haline,
            total,
            reference_valid,
            hb._region_weights(temperature),
        )
        comparisons[dataset] = {
            "depth": depth,
            "thermal": thermal,
            "haline": haline,
            "total": total,
            "reference_valid": reference_valid,
            "weights": hb._region_weights(ref_t_cube),
            "ancestors": [
                entries["thetao"]["filename"],
                entries["so"]["filename"],
            ],
        }
        cubes = iris.cube.CubeList(
            [
                _output_cube(
                    value,
                    name,
                    "1" if name in ("coverage", "cancellation") else "kg m-3",
                    depth,
                )
                for name, value in metrics.items()
            ]
        )
        path = get_diagnostic_filename(f"density_compensation_{dataset}", cfg)
        iris.save(cubes, path)
        plot_path = _plot(dataset, depth, metrics, cfg)
        record = {
            "caption": (
                f"{dataset} density bias attributed to temperature "
                "and salinity differences against WOA"
            ),
            "statistics": ["mean", "rms"],
            "domains": ["global"],
            "plot_types": ["vert"],
            "authors": ["beucher_romain"],
            "ancestors": [
                entries["thetao"]["filename"],
                entries["so"]["filename"],
                references["thetao"]["filename"],
                references["so"]["filename"],
            ],
        }
        with ProvenanceLogger(cfg) as provenance:
            provenance.log(path, record)
            provenance.log(plot_path, record)

    if len(comparisons) > 1:
        path, plot_path = _multimodel_comparison(models, comparisons, cfg)
        record = {
            "caption": (
                "Multi-model density comparison against WOA; "
                "model names follow model_index order"
            ),
            "statistics": ["mean", "rms"],
            "domains": ["global"],
            "plot_types": ["vert"],
            "authors": ["beucher_romain"],
            "ancestors": list(
                dict.fromkeys(
                    filename
                    for model in models
                    for filename in comparisons[model]["ancestors"]
                ),
            )
            + [references["thetao"]["filename"], references["so"]["filename"]],
        }
        with ProvenanceLogger(cfg) as provenance:
            provenance.log(path, record)
            provenance.log(plot_path, record)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
