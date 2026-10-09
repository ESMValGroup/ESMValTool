"""ML-based downscaling — extra diagnostics.

This diagnostic implements the additional figures that are not already
covered by ``ml_downscaling_evaluation.py``. It is
dispatched on ``cfg['analysis_type']``; the same script handles all of:

  * ``climatology``                — time-mean spatial maps for each variable
                                      and method, with ensemble std panel.
  * ``percentile_map``             — 99th-percentile spatial maps (pr and tas).
  * ``cdd``                        — consecutive dry days for precipitation.
  * ``rx1day``                     — annual maximum daily precipitation.
  * ``conditional_rank_histogram`` — rank histograms restricted to time steps
                                      where domain-mean pr exceeds its 99th
                                      percentile.
  * ``return_level``               — return-level plot for 1-, 2-, and 3-day
                                      precipitation accumulations at one
                                      reference grid cell per region.
  * ``temporal_spectrum``          — temporal power spectrum (per region,
                                      area-averaged), pr and tas.
  * ``cross_variable_correlation`` — pointwise Pearson correlation between
                                      pairs of variables (tas–huss, huss–pr,
                                      uas–vas).
  * ``case_study``                 — synoptic snapshot panel: LR / ref /
                                      AFM / PC-AFM ensemble mean / ensemble
                                      std for all output variables.
  * ``quantile_mae_table``         — absolute Quantile-MAE values for the
                                      baseline reference (AFM-baseline-2M),
                                      to accompany the existing relative
                                      heatmaps.

The diagnostic re-uses the data-loading utilities of
``ml_downscaling_evaluation.py`` where possible.
"""

import logging
import os
from pathlib import Path

import iris
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvaltool.diag_scripts.shared._base import ProvenanceLogger

# Re-use helpers from the main evaluation diagnostic where convenient
try:
    from esmvaltool.diag_scripts.regional_downscaling.ml_downscaling_evaluation import (  # noqa: E501
        load_ensemble_data,
        update_variable_units,
        UNITS,
        avg_pool_2d_weighted,
        calculate_wbgt,
        calculate_wind_speed,
    )
except ImportError:  # pragma: no cover - defensive
    load_ensemble_data = None
    UNITS = {
        "pr": "mm h-1", "huss": "g kg-1", "tas": "degC",
        "ps": "hPa", "uas": "m s-1", "vas": "m s-1",
    }

logger = logging.getLogger(os.path.basename(__file__))

# Colorblind-friendly diverging palette consistent with the heatmap update
DIVERGING_CMAP = "RdBu_r"
SEQUENTIAL_CMAP = "viridis"

# Reference sites for return-level plots: one per region (lon, lat, label)
REFERENCE_SITES = {
    "CE": (8.7, 49.5, "Frankfurt"),
    "IBE": (-3.7, 40.4, "Madrid"),
    "SCA": (18.1, 59.3, "Stockholm"),
}


# ---------------------------------------------------------------------------
# Provenance & helpers
# ---------------------------------------------------------------------------


def _provenance(cfg, plot_file, caption, statistics=None, plot_types=None):
    """Write an ESMValTool provenance record for ``plot_file``."""
    record = {
        "caption": caption,
        "statistics": statistics or ["other"],
        "domains": ["reg"],
        "plot_types": plot_types or ["map"],
        "authors": ["debeire_kevin"],
        "references": [],
        "plot_file": plot_file,
        "ancestors": [
            d["filename"] for d in cfg.get("input_data", {}).values()
        ],
    }
    with ProvenanceLogger(cfg) as plog:
        plog.log(plot_file, record)


def _save_figure(fig, cfg, name, caption, statistics=None, plot_types=None):
    """Save figure to ``cfg['plot_dir']`` (png+pdf) and write provenance."""
    plot_file = os.path.join(cfg["plot_dir"], f"{name}.png")
    fig.savefig(plot_file, dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(plot_file.replace(".png", ".pdf"),
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    _provenance(cfg, plot_file, caption,
                statistics=statistics, plot_types=plot_types)
    logger.info("Saved %s", plot_file)
    return plot_file


def _extract_extent(truth_cube):
    """Return (lon_min, lon_max, lat_min, lat_max) or a default."""
    try:
        lat = truth_cube.coord("latitude").points
        lon = truth_cube.coord("longitude").points
        if lat.ndim == 2:
            lat = np.mean(lat, axis=1)
        if lon.ndim == 2:
            lon = np.mean(lon, axis=0)
        return [float(lon.min()), float(lon.max()),
                float(lat.min()), float(lat.max())]
    except iris.exceptions.CoordinateNotFoundError:
        return [-180.0, 180.0, -90.0, 90.0]


def _load_method_truth(grouped_data, var_name, reference_name, ml_methods):
    """Load reference cube + per-method ensemble arrays for a variable.

    Returns (truth_cube, truth_array, preds_methods, method_names) where
    ``preds_methods`` is a list of (n_ens, T, X, Y) arrays.
    """
    var_datasets = grouped_data[var_name]
    ref_datasets = [d for d in var_datasets if d["dataset"] == reference_name]
    if not ref_datasets:
        return None, None, [], []
    truth_cube = iris.load_cube(ref_datasets[0]["filename"])
    if truth_cube.ndim == 4:
        truth_cube = truth_cube[:, :, :, 0]
        truth_cube.transpose([0, 2, 1])
    elif truth_cube.ndim == 3:
        truth_cube.transpose([0, 2, 1])
    truth_array = truth_cube.data

    preds_methods, method_names = [], []
    for method in ml_methods:
        method_datasets = [d for d in var_datasets if d["dataset"] == method]
        if not method_datasets:
            continue
        if load_ensemble_data is None:
            logger.warning("load_ensemble_data unavailable; skipping %s",
                           method)
            continue
        pred_ens = load_ensemble_data(method_datasets)
        preds_methods.append(pred_ens)
        method_names.append(method)
    return truth_cube, truth_array, preds_methods, method_names


# ---------------------------------------------------------------------------
# Scalar-summary helpers used to annotate map subplots
# ---------------------------------------------------------------------------


def _field_metrics(pred, ref):
    """Domain-mean, bias, RMSE and spatial Pearson r of ``pred`` vs ``ref``.

    NaN-aware. Returns ``(mean, bias, rmse, corr)`` over the pixels where both
    fields are finite. If no overlap, returns four NaNs.
    """
    p = np.asarray(pred, dtype=float).ravel()
    r = np.asarray(ref, dtype=float).ravel()
    mask = np.isfinite(p) & np.isfinite(r)
    if not np.any(mask):
        return np.nan, np.nan, np.nan, np.nan
    p, r = p[mask], r[mask]
    mean = float(np.mean(p))
    bias = float(np.mean(p - r))
    rmse = float(np.sqrt(np.mean((p - r) ** 2)))
    if p.size > 1 and np.std(p) > 0 and np.std(r) > 0:
        corr = float(np.corrcoef(p, r)[0, 1])
    else:
        corr = np.nan
    return mean, bias, rmse, corr


def _g(x, p=2):
    """Compact float formatter; ``NaN`` falls back to a short string."""
    return "NaN" if not np.isfinite(x) else f"{x:.{p}g}"


def _ref_metric_label(field, precision=2):
    """One-line label for a reference panel: just the domain-mean."""
    val = float(np.nanmean(field)) if np.any(np.isfinite(field)) else np.nan
    return rf"$\mu$={_g(val, precision)}"


def _pred_metric_label(pred, ref, precision=2, include_rmse=True):
    """One-line label for a method panel: mean, bias, [RMSE], spatial r."""
    mean, bias, rmse, corr = _field_metrics(pred, ref)
    parts = [rf"$\mu$={_g(mean, precision)}",
             rf"bias={_g(bias, precision)}"]
    if include_rmse:
        parts.append(rf"RMSE={_g(rmse, precision)}")
    parts.append(rf"r={_g(corr, 2)}")
    return "  ".join(parts)


# ---------------------------------------------------------------------------
# 1. Climatology maps + ensemble-std panel
# ---------------------------------------------------------------------------


def plot_climatology(truth, preds_methods, method_names, var_name, extent, cfg):
    """Time-mean (climatology) maps + ensemble standard-deviation maps.

    Layout: n rows = 1 (reference) + len(methods); n cols = 2
    (mean, ensemble std). The reference has no ensemble std; that cell is
    left blank.
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    units = UNITS.get(var_name, "")
    region = cfg.get("region_name", "")
    n_rows = 1 + len(method_names)
    fig = plt.figure(figsize=(8.0, 3.0 * n_rows))
    projection = ccrs.PlateCarree()

    # Time means
    truth_mean = np.mean(truth, axis=0)
    method_means = [np.mean(np.mean(p, axis=0), axis=0)
                    for p in preds_methods]
    method_ens_std = [np.std(np.mean(p, axis=1), axis=0)
                      for p in preds_methods]

    vmin = float(np.nanmin([truth_mean.min(),
                            *[m.min() for m in method_means]]))
    vmax = float(np.nanmax([truth_mean.max(),
                            *[m.max() for m in method_means]]))

    for row in range(n_rows):
        # Mean column
        ax = fig.add_subplot(n_rows, 2, 2 * row + 1, projection=projection)
        if row == 0:
            data = truth_mean
            label = "Reference"
            metric_str = _ref_metric_label(truth_mean)
        else:
            data = method_means[row - 1]
            label = method_names[row - 1]
            metric_str = _pred_metric_label(data, truth_mean)
        im = ax.imshow(data, cmap=SEQUENTIAL_CMAP, vmin=vmin, vmax=vmax,
                       origin="lower", extent=extent, transform=projection,
                       aspect="auto")
        ax.coastlines(resolution="50m", linewidth=0.8, color="black")
        ax.add_feature(cfeature.BORDERS, linewidth=0.5, alpha=0.5)
        ax.text(-0.10, 0.5, label, transform=ax.transAxes, rotation=90,
                ha="right", va="center", fontsize=12, fontweight="bold")
        if row == 0:
            ax.set_title(f"Time mean ({units})\n{metric_str}",
                         fontsize=11, fontweight="bold")
        else:
            ax.set_title(metric_str, fontsize=10)
        if row == n_rows - 1:
            fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.05,
                         fraction=0.06)

        # Ensemble std column
        if row == 0:
            # Deterministic reference → no ensemble std. Use a plain axes
            # (no cartopy projection) so we get a clean blank panel
            # instead of a global map.
            ax = fig.add_subplot(n_rows, 2, 2 * row + 2)
            ax.text(0.5, 0.5, "N/A\n(reference is deterministic)",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=11, color="#555555")
            ax.set_axis_off()
            ax.set_title(f"Ensemble std ({units})",
                         fontsize=12, fontweight="bold")
        else:
            ax = fig.add_subplot(n_rows, 2, 2 * row + 2, projection=projection)
            data = method_ens_std[row - 1]
            im = ax.imshow(data, cmap="cividis", origin="lower",
                           extent=extent, transform=projection, aspect="auto")
            ax.coastlines(resolution="50m", linewidth=0.8, color="black")
            ax.add_feature(cfeature.BORDERS, linewidth=0.5, alpha=0.5)
            if row == n_rows - 1:
                fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.05,
                             fraction=0.06)

    fig.suptitle(f"{var_name} climatology — {region}", fontsize=14,
                 fontweight="bold", y=1.02)
    fig.subplots_adjust(wspace=0.05, hspace=0.25)
    caption = (
        f"Time-mean climatology and ensemble standard deviation of {var_name} "
        f"over the {region} domain. Reference is deterministic, so its "
        f"ensemble std panel is blank."
    )
    _save_figure(fig, cfg, f"climatology_{var_name}", caption,
                 statistics=["mean", "stddev"], plot_types=["map"])


# ---------------------------------------------------------------------------
# 2. 99th-percentile spatial maps  (pr & tas)
# ---------------------------------------------------------------------------


def plot_percentile_map(truth, preds_methods, method_names, var_name, extent,
                        cfg, percentile=99.0):
    """Spatial maps of the qth percentile (default 99) per pixel."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    units = UNITS.get(var_name, "")
    region = cfg.get("region_name", "")
    n_cols = 1 + len(method_names)
    fig = plt.figure(figsize=(4.2 * n_cols, 3.6))
    projection = ccrs.PlateCarree()

    truth_q = np.nanpercentile(truth, percentile, axis=0)
    preds_q = []
    for p in preds_methods:
        ens_mean = np.mean(p, axis=0)
        preds_q.append(np.nanpercentile(ens_mean, percentile, axis=0))

    vmin = float(np.nanmin([truth_q.min(), *[q.min() for q in preds_q]]))
    vmax = float(np.nanmax([truth_q.max(), *[q.max() for q in preds_q]]))
    cmap = "magma_r" if var_name == "pr" else "inferno"

    for col in range(n_cols):
        ax = fig.add_subplot(1, n_cols, col + 1, projection=projection)
        if col == 0:
            data, title = truth_q, "Reference"
            metric_str = _ref_metric_label(truth_q)
        else:
            data = preds_q[col - 1]
            title = method_names[col - 1]
            metric_str = _pred_metric_label(data, truth_q)
        im = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax,
                       origin="lower", extent=extent, transform=projection,
                       aspect="auto")
        ax.coastlines(resolution="50m", linewidth=0.8, color="black")
        ax.add_feature(cfeature.BORDERS, linewidth=0.5, alpha=0.5)
        ax.set_title(f"{title}\n{metric_str}", fontsize=10, fontweight="bold")

    cbar = fig.colorbar(im, ax=fig.get_axes(), orientation="horizontal",
                        fraction=0.04, pad=0.08, aspect=40)
    cbar.set_label(f"{var_name} {percentile:.0f}th percentile ({units})",
                   fontsize=10)
    fig.suptitle(f"{var_name} P{percentile:.0f} — {region}",
                 fontsize=13, fontweight="bold", y=1.08)
    caption = (
        f"Per-pixel {percentile:.0f}th percentile of {var_name} over the "
        f"test period for the {region} domain. For ML methods the percentile "
        f"is computed on the ensemble-mean field."
    )
    _save_figure(fig, cfg, f"percentile{int(percentile)}_{var_name}", caption,
                 statistics=["perc"], plot_types=["map"])


# ---------------------------------------------------------------------------
# 3. CDD (consecutive dry days)
# ---------------------------------------------------------------------------


def _daily_aggregate(arr, hours_per_step=3):
    """Aggregate sub-daily timeseries to daily totals along axis 0."""
    steps_per_day = max(1, 24 // hours_per_step)
    n = (arr.shape[0] // steps_per_day) * steps_per_day
    if n == 0:
        return arr.copy()
    arr = arr[:n]
    new_shape = (n // steps_per_day, steps_per_day) + arr.shape[1:]
    return arr.reshape(new_shape).sum(axis=1) * hours_per_step


def _cdd_field(daily_pr, threshold=1.0):
    """Maximum consecutive dry days per pixel (daily_pr in mm/day)."""
    T = daily_pr.shape[0]
    dry = daily_pr < threshold
    # Pixelwise scan
    out = np.zeros(daily_pr.shape[1:], dtype=np.float32)
    cur = np.zeros_like(out)
    for t in range(T):
        cur = np.where(dry[t], cur + 1, 0)
        out = np.maximum(out, cur)
    return out


def plot_cdd(truth, preds_methods, method_names, var_name, extent, cfg):
    """Spatial map of maximum consecutive dry days (CDD) for pr."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    if var_name != "pr":
        logger.warning("CDD is computed for pr only, got %s", var_name)
        return

    # Aggregate to daily totals (pr in mm/h → mm/day)
    truth_daily = _daily_aggregate(truth)
    truth_cdd = _cdd_field(truth_daily)
    preds_cdd = []
    for p in preds_methods:
        ens_mean = np.mean(p, axis=0)
        preds_cdd.append(_cdd_field(_daily_aggregate(ens_mean)))

    n_cols = 1 + len(method_names)
    fig = plt.figure(figsize=(4.2 * n_cols, 3.6))
    projection = ccrs.PlateCarree()
    vmax = float(np.nanmax([truth_cdd.max(), *[m.max() for m in preds_cdd]]))

    for col in range(n_cols):
        ax = fig.add_subplot(1, n_cols, col + 1, projection=projection)
        if col == 0:
            data, title = truth_cdd, "Reference"
            metric_str = _ref_metric_label(truth_cdd)
        else:
            data = preds_cdd[col - 1]
            title = method_names[col - 1]
            metric_str = _pred_metric_label(data, truth_cdd)
        im = ax.imshow(data, cmap="YlOrBr", vmin=0, vmax=vmax,
                       origin="lower", extent=extent, transform=projection,
                       aspect="auto")
        ax.coastlines(resolution="50m", linewidth=0.8, color="black")
        ax.add_feature(cfeature.BORDERS, linewidth=0.5, alpha=0.5)
        ax.set_title(f"{title}\n{metric_str}", fontsize=10, fontweight="bold")

    cbar = fig.colorbar(im, ax=fig.get_axes(), orientation="horizontal",
                        fraction=0.04, pad=0.08, aspect=40)
    cbar.set_label("Max consecutive dry days (threshold 1 mm)", fontsize=10)
    fig.suptitle(f"CDD — {cfg.get('region_name', '')}",
                 fontsize=13, fontweight="bold", y=1.08)
    caption = (
        "Maximum number of consecutive dry days (daily precipitation < 1 mm) "
        "per pixel over the test period, computed on the ensemble mean for "
        "ML methods."
    )
    _save_figure(fig, cfg, "cdd_pr", caption,
                 statistics=["other"], plot_types=["map"])


# ---------------------------------------------------------------------------
# 4. Rx1day (annual maximum daily pr)
# ---------------------------------------------------------------------------


def plot_rx1day(truth, preds_methods, method_names, var_name, extent, cfg):
    """Spatial map of the annual maximum daily precipitation."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    if var_name != "pr":
        logger.warning("Rx1day is computed for pr only, got %s", var_name)
        return

    truth_daily = _daily_aggregate(truth)
    truth_max = np.max(truth_daily, axis=0)
    preds_max = []
    for p in preds_methods:
        ens_mean = np.mean(p, axis=0)
        preds_max.append(np.max(_daily_aggregate(ens_mean), axis=0))

    n_cols = 1 + len(method_names)
    fig = plt.figure(figsize=(4.2 * n_cols, 3.6))
    projection = ccrs.PlateCarree()
    vmax = float(np.nanmax([truth_max.max(), *[m.max() for m in preds_max]]))

    for col in range(n_cols):
        ax = fig.add_subplot(1, n_cols, col + 1, projection=projection)
        if col == 0:
            data, title = truth_max, "Reference"
            metric_str = _ref_metric_label(truth_max)
        else:
            data = preds_max[col - 1]
            title = method_names[col - 1]
            metric_str = _pred_metric_label(data, truth_max)
        im = ax.imshow(data, cmap="magma_r", vmin=0, vmax=vmax,
                       origin="lower", extent=extent, transform=projection,
                       aspect="auto")
        ax.coastlines(resolution="50m", linewidth=0.8, color="black")
        ax.add_feature(cfeature.BORDERS, linewidth=0.5, alpha=0.5)
        ax.set_title(f"{title}\n{metric_str}", fontsize=10, fontweight="bold")

    cbar = fig.colorbar(im, ax=fig.get_axes(), orientation="horizontal",
                        fraction=0.04, pad=0.08, aspect=40)
    cbar.set_label("Annual maximum daily precipitation (mm)", fontsize=10)
    fig.suptitle(f"Rx1day — {cfg.get('region_name', '')}",
                 fontsize=13, fontweight="bold", y=1.08)
    caption = (
        "Per-pixel annual maximum daily precipitation total over the test "
        "year, computed on the ensemble mean for ML methods."
    )
    _save_figure(fig, cfg, "rx1day_pr", caption,
                 statistics=["other"], plot_types=["map"])


# ---------------------------------------------------------------------------
# 5. Conditional rank histogram (>99th percentile pr days)
# ---------------------------------------------------------------------------


def plot_conditional_rank_histogram(truth, preds_methods, method_names,
                                    var_name, cfg):
    """Rank histograms restricted to the time steps where the domain-mean
    precipitation exceeds its 99th percentile.
    """
    if var_name != "pr":
        logger.warning("Conditional rank hist computed for pr only, got %s",
                       var_name)
        return

    region = cfg.get("region_name", "")
    domain_mean = np.nanmean(truth, axis=(1, 2))
    threshold = np.nanpercentile(domain_mean, 99.0)
    extreme_mask = domain_mean >= threshold
    logger.info("Conditional rank hist: %d / %d time steps are extreme",
                int(extreme_mask.sum()), domain_mean.size)
    if extreme_mask.sum() < 5:
        logger.warning("Too few extreme time steps; skipping plot.")
        return

    n_methods = len(method_names)
    fig, axes = plt.subplots(1, max(1, n_methods),
                             figsize=(4.5 * max(1, n_methods), 4.0),
                             squeeze=False)
    axes = axes[0]
    for ax, method, pred_ens in zip(axes, method_names, preds_methods):
        n_ens = pred_ens.shape[0]
        # Restrict to extreme time steps
        truth_e = truth[extreme_mask]
        pred_e = pred_ens[:, extreme_mask]
        # Rank of truth within ensemble per pixel and time
        ranks = (pred_e < truth_e[np.newaxis]).sum(axis=0)
        ranks = ranks.flatten()
        ranks = ranks[np.isfinite(ranks)]
        hist, edges = np.histogram(ranks, bins=np.arange(n_ens + 2))
        hist = hist / hist.sum() if hist.sum() > 0 else hist
        ax.bar(edges[:-1], hist, width=1.0, color="#2166ac",
               edgecolor="#053061")
        ax.axhline(1.0 / (n_ens + 1), color="black", linestyle="--",
                   linewidth=1.2)
        mcb = np.sum((hist - 1.0 / (n_ens + 1)) ** 2) / (n_ens + 1)
        ax.set_title(f"{method}\nExtreme MCB = {mcb:.3f}",
                     fontsize=11, fontweight="bold")
        ax.set_xlabel("Rank")
        ax.set_ylabel("Relative frequency")
    fig.suptitle(f"Rank histogram | pr > P99 — {region}",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    caption = (
        f"Rank histograms computed only on time steps where the "
        f"domain-mean precipitation exceeds its 99th percentile. "
        f"MCB shown per method; perfect calibration corresponds to a flat "
        f"histogram at 1/(n_ens+1) = {1.0 / 13.0:.3f}."
    )
    _save_figure(fig, cfg, "conditional_rank_histogram_pr", caption,
                 statistics=["other"], plot_types=["other"])


# ---------------------------------------------------------------------------
# 6. Return-level plots
# ---------------------------------------------------------------------------


def plot_return_level(truth, preds_methods, method_names, var_name, extent,
                      cfg):
    """Return-level plot for 1-, 2-, 3-day pr accumulations at one site."""
    if var_name != "pr":
        logger.warning("Return-level computed for pr only, got %s", var_name)
        return

    region_key = cfg.get("region_key", "CE")
    lon_site, lat_site, site_label = REFERENCE_SITES.get(
        region_key, (None, None, "site")
    )

    # Nearest-pixel selection from the spatial extent
    if lon_site is None:
        ix, iy = truth.shape[2] // 2, truth.shape[1] // 2
    else:
        lon_min, lon_max, lat_min, lat_max = extent
        ix = int(np.clip(
            round((lon_site - lon_min) / (lon_max - lon_min) *
                  (truth.shape[2] - 1)), 0, truth.shape[2] - 1
        ))
        iy = int(np.clip(
            round((lat_site - lat_min) / (lat_max - lat_min) *
                  (truth.shape[1] - 1)), 0, truth.shape[1] - 1
        ))

    accumulations = [1, 2, 3]
    fig, axes = plt.subplots(1, len(accumulations),
                             figsize=(4.5 * len(accumulations), 4.0),
                             sharey=True)
    for j, (ax, n_days) in enumerate(zip(axes, accumulations)):
        truth_daily = _daily_aggregate(truth)[:, iy, ix]
        # Sliding-window accumulation
        if n_days > 1:
            kernel = np.ones(n_days)
            truth_acc = np.convolve(truth_daily, kernel, mode="valid")
        else:
            truth_acc = truth_daily
        ranks = np.arange(1, truth_acc.size + 1)
        sorted_truth = np.sort(truth_acc)[::-1]
        emp_rl = sorted_truth
        return_periods = (ranks / (sorted_truth.size + 1.0)) ** -1 / 365.0
        ax.plot(return_periods, emp_rl, "k-", lw=2.0, label="Reference")
        for method, pred_ens in zip(method_names, preds_methods):
            ens_mean = np.mean(pred_ens, axis=0)
            p_daily = _daily_aggregate(ens_mean)[:, iy, ix]
            if n_days > 1:
                p_acc = np.convolve(p_daily, np.ones(n_days), mode="valid")
            else:
                p_acc = p_daily
            sorted_p = np.sort(p_acc)[::-1]
            ax.plot(return_periods[:sorted_p.size], sorted_p, "--", lw=2.0,
                    label=method)
        ax.set_xscale("log")
        ax.set_xlabel("Return period (years)")
        if j == 0:
            ax.set_ylabel("Precipitation accumulation (mm)")
        ax.set_title(f"{n_days}-day", fontsize=11, fontweight="bold")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)

    fig.suptitle(f"Return-level @ {site_label} — {cfg.get('region_name', '')}",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    caption = (
        f"Empirical return levels for {accumulations}-day precipitation "
        f"accumulation at {site_label} ({lat_site}°N, {lon_site}°E). "
        f"Sorted decreasing values are plotted against empirical return "
        f"period. Only ~1 year of test data is available, so return periods "
        f"beyond a few months should be read as indicative."
    )
    _save_figure(fig, cfg, "return_level_pr", caption,
                 statistics=["other"], plot_types=["other"])


# ---------------------------------------------------------------------------
# 7. Temporal power spectrum
# ---------------------------------------------------------------------------


def plot_temporal_spectrum(truth, preds_methods, method_names, var_name, cfg):
    """Domain-averaged temporal power spectrum (frequency = 1/step)."""
    region = cfg.get("region_name", "")
    dt_hours = float(cfg.get("dt_hours", 3.0))
    truth_ts = np.nanmean(truth, axis=(1, 2))  # (T,)
    truth_ts = truth_ts - truth_ts.mean()
    n = truth_ts.size
    freqs = np.fft.rfftfreq(n, d=dt_hours)  # cycles/hour
    truth_psd = np.abs(np.fft.rfft(truth_ts)) ** 2

    fig, ax = plt.subplots(figsize=(8, 5))
    # Reference plotted first (lower zorder, lighter) so the method curves
    # render on top and remain readable.
    ax.loglog(freqs[1:], truth_psd[1:], color="black", lw=1.4, alpha=0.55,
              label="Reference", zorder=2)
    for method, pred_ens in zip(method_names, preds_methods):
        ens_mean = np.mean(pred_ens, axis=0)
        ts = np.nanmean(ens_mean, axis=(1, 2))
        ts = ts - ts.mean()
        psd = np.abs(np.fft.rfft(ts)) ** 2
        ax.loglog(freqs[1:], psd[1:], "--", lw=1.1, alpha=0.95, label=method,
                  zorder=3)
    ax.set_xlabel("Frequency (1/hour)")
    ax.set_ylabel("Power spectral density (arb. units)")
    ax.set_title(f"Temporal PSD — {var_name} — {region}",
                 fontsize=13, fontweight="bold")
    ax.legend(loc="lower left", fontsize=9)
    ax.grid(True, which="both", alpha=0.3)
    # Annotate diurnal & weekly bands at the bottom so the labels do not
    # collide with the legend or the title.
    ymin, ymax = ax.get_ylim()
    label_y = ymin * (ymax / ymin) ** 0.04  # ~4% above the bottom in log scale
    for f, label in [(1 / 24.0, "diurnal"), (1 / (24.0 * 7), "weekly")]:
        ax.axvline(f, color="#777777", linestyle=":", linewidth=1.0,
                   zorder=1)
        ax.text(f, label_y, f" {label}", rotation=90, va="bottom",
                ha="left", color="#222222", fontsize=9)
    plt.tight_layout()
    caption = (
        f"Temporal power spectral density of domain-averaged {var_name} "
        f"over the {region} domain. Diurnal and weekly bands marked for "
        f"reference. ML methods that independently downscale each time step "
        f"are expected to underestimate the low-frequency power."
    )
    _save_figure(fig, cfg, f"temporal_spectrum_{var_name}", caption,
                 statistics=["other"], plot_types=["other"])


# ---------------------------------------------------------------------------
# 8. Cross-variable correlation maps
# ---------------------------------------------------------------------------


def _pointwise_correlation(a, b):
    """Pixelwise Pearson r between time series in two (T, X, Y) arrays."""
    a = a - np.nanmean(a, axis=0, keepdims=True)
    b = b - np.nanmean(b, axis=0, keepdims=True)
    num = np.nansum(a * b, axis=0)
    den = np.sqrt(np.nansum(a ** 2, axis=0) * np.nansum(b ** 2, axis=0))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan)


def plot_cross_variable_correlation(grouped_data, reference_name, ml_methods,
                                    cfg, extent):
    """Pointwise Pearson r between (tas, huss), (huss, pr), (uas, vas)."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    pairs = [("tas", "huss"), ("huss", "pr"), ("uas", "vas")]
    needed = sorted({v for pair in pairs for v in pair})
    if not all(v in grouped_data for v in needed):
        logger.warning("Missing variables for cross-corr: need %s, have %s",
                       needed, list(grouped_data.keys()))
        return

    # Load reference + method arrays for each needed variable
    truth_by_var, preds_by_var, method_names = {}, {}, None
    for v in needed:
        out = _load_method_truth(grouped_data, v, reference_name, ml_methods)
        truth_cube, truth_arr, preds, names = out
        truth_by_var[v] = truth_arr
        preds_by_var[v] = preds
        if method_names is None:
            method_names = names
            extent = _extract_extent(truth_cube)

    projection = ccrs.PlateCarree()
    n_rows = len(pairs)
    n_cols = 1 + len(method_names)
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(3.6 * n_cols, 3.0 * n_rows),
                             subplot_kw={"projection": projection},
                             squeeze=False)

    for r, (v1, v2) in enumerate(pairs):
        ref_corr = _pointwise_correlation(truth_by_var[v1], truth_by_var[v2])
        ax = axes[r, 0]
        im = ax.imshow(ref_corr, cmap="RdBu_r", vmin=-1.0, vmax=1.0,
                       origin="lower", extent=extent, transform=projection,
                       aspect="auto")
        ax.coastlines(resolution="50m", linewidth=0.8, color="black")
        ax.add_feature(cfeature.BORDERS, linewidth=0.5, alpha=0.5)
        ref_label = _ref_metric_label(ref_corr)
        ax.set_title(f"Reference — {v1} vs {v2}\n{ref_label}",
                     fontsize=10, fontweight="bold")
        for c, method in enumerate(method_names, start=1):
            ax = axes[r, c]
            p1 = np.mean(preds_by_var[v1][c - 1], axis=0)
            p2 = np.mean(preds_by_var[v2][c - 1], axis=0)
            corr = _pointwise_correlation(p1, p2)
            im = ax.imshow(corr, cmap="RdBu_r", vmin=-1.0, vmax=1.0,
                           origin="lower", extent=extent,
                           transform=projection, aspect="auto")
            ax.coastlines(resolution="50m", linewidth=0.8, color="black")
            ax.add_feature(cfeature.BORDERS, linewidth=0.5, alpha=0.5)
            method_label = _pred_metric_label(corr, ref_corr,
                                              include_rmse=False)
            header = f"{method} — {v1} vs {v2}" if r == 0 else \
                     f"{v1} vs {v2}"
            ax.set_title(f"{header}\n{method_label}",
                         fontsize=10, fontweight="bold")

    cbar = fig.colorbar(im, ax=axes, orientation="horizontal",
                        fraction=0.04, pad=0.06, aspect=40)
    cbar.set_label("Pointwise Pearson r", fontsize=10)
    fig.suptitle(f"Cross-variable correlation — {cfg.get('region_name', '')}",
                 fontsize=13, fontweight="bold")
    caption = (
        "Pointwise Pearson correlation between pairs of variables over the "
        "test period: (tas, huss), (huss, pr), (uas, vas). Reference and "
        "ensemble-mean ML predictions are shown side-by-side."
    )
    _save_figure(fig, cfg, "cross_variable_correlation", caption,
                 statistics=["corr"], plot_types=["map"])


# ---------------------------------------------------------------------------
# 9. Case-study snapshot panel
# ---------------------------------------------------------------------------


def plot_case_study(grouped_data, reference_name, ml_methods, cfg):
    """Single random time step, all variables, side-by-side."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    rng = np.random.default_rng(cfg.get("case_study_seed", 42))
    variables = ["pr", "huss", "tas", "ps", "uas", "vas"]
    variables = [v for v in variables if v in grouped_data]
    if not variables:
        return

    truth_by_var, preds_by_var, method_names = {}, {}, None
    extent = None
    ref_time_cube = None
    for v in variables:
        out = _load_method_truth(grouped_data, v, reference_name, ml_methods)
        truth_cube, truth_arr, preds, names = out
        if truth_arr is None:
            continue
        truth_by_var[v] = truth_arr
        preds_by_var[v] = preds
        if method_names is None:
            method_names = names
            extent = _extract_extent(truth_cube)
            ref_time_cube = truth_cube

    if extent is None:
        return

    # Pick a time step where the domain-mean pr is in the top decile
    if "pr" in truth_by_var:
        dmean = np.nanmean(truth_by_var["pr"], axis=(1, 2))
        candidates = np.where(dmean >= np.nanpercentile(dmean, 90.0))[0]
        t_idx = int(rng.choice(candidates))
    else:
        T = next(iter(truth_by_var.values())).shape[0]
        t_idx = int(rng.integers(T))

    # Resolve t_idx to an absolute timestamp from the reference cube.
    timestamp_label = f"t={t_idx}"
    try:
        tcoord = ref_time_cube.coord("time")
        dt = tcoord.units.num2date(tcoord.points[t_idx])
        timestamp_label = dt.strftime("%Y-%m-%d %H:%M UTC")
    except Exception:
        logger.warning("Could not resolve case-study timestamp; "
                       "falling back to integer index.")

    projection = ccrs.PlateCarree()
    n_cols = 3 + len(method_names)  # LR / Ref / [methods mean, ens std]
    n_rows = len(variables)
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(3.0 * n_cols, 2.6 * n_rows),
                             subplot_kw={"projection": projection},
                             squeeze=False)

    for r, v in enumerate(variables):
        truth_field = truth_by_var[v][t_idx]
        cmap = ("magma_r" if v == "pr"
                else "viridis" if v in {"huss", "tas", "ps"}
                else "RdBu_r")
        vmin = float(np.nanmin(truth_field))
        vmax = float(np.nanmax(truth_field))
        units = UNITS.get(v, "")

        im_data = None
        im_std = None

        # Column 0: a coarsened "LR" view (5x box average for visualisation)
        ax = axes[r, 0]
        try:
            pooled, _ = avg_pool_2d_weighted(
                truth_field[np.newaxis], np.linspace(extent[2], extent[3],
                                                     truth_field.shape[0])
            )
            lr_view = pooled[0]
        except Exception:
            lr_view = truth_field
        im_data = ax.imshow(lr_view, cmap=cmap, vmin=vmin, vmax=vmax,
                            origin="lower", extent=extent,
                            transform=projection, aspect="auto")
        ax.coastlines(resolution="50m", linewidth=0.6, color="black")
        if r == 0:
            ax.set_title("LR view", fontsize=10, fontweight="bold")
        ax.text(-0.12, 0.5, f"{v}\n({units})", transform=ax.transAxes,
                rotation=90, ha="right", va="center",
                fontsize=10, fontweight="bold")

        ax = axes[r, 1]
        im_data = ax.imshow(truth_field, cmap=cmap, vmin=vmin, vmax=vmax,
                            origin="lower", extent=extent,
                            transform=projection, aspect="auto")
        ax.coastlines(resolution="50m", linewidth=0.6, color="black")
        if r == 0:
            ax.set_title("Reference", fontsize=10, fontweight="bold")

        for c, method in enumerate(method_names):
            pred = preds_by_var[v][c]  # (n_ens, T, X, Y)
            ens_mean = np.mean(pred[:, t_idx], axis=0)
            ax = axes[r, 2 + c]
            im_data = ax.imshow(ens_mean, cmap=cmap, vmin=vmin, vmax=vmax,
                                origin="lower", extent=extent,
                                transform=projection, aspect="auto")
            ax.coastlines(resolution="50m", linewidth=0.6, color="black")
            if r == 0:
                ax.set_title(method, fontsize=10, fontweight="bold")
        # Last column: ensemble std of the most expressive method
        if method_names:
            std_field = np.std(preds_by_var[v][-1][:, t_idx], axis=0)
            ax = axes[r, -1]
            im_std = ax.imshow(std_field, cmap="cividis", origin="lower",
                               extent=extent, transform=projection,
                               aspect="auto")
            ax.coastlines(resolution="50m", linewidth=0.6, color="black")
            if r == 0:
                ax.set_title(f"{method_names[-1]} ens std",
                             fontsize=10, fontweight="bold")

        # Row-shared colorbars: one for the data cmap (cols 0..n-2),
        # one for the ensemble-std cmap (col n-1). Placed to the right.
        data_axes = list(axes[r, :n_cols - 1])
        cbar = fig.colorbar(im_data, ax=data_axes, location="right",
                            pad=0.005, fraction=0.018, shrink=0.85)
        cbar.set_label(units, fontsize=8)
        cbar.ax.tick_params(labelsize=7)
        if im_std is not None:
            cbar_std = fig.colorbar(im_std, ax=axes[r, -1], location="right",
                                    pad=0.02, fraction=0.05, shrink=0.85)
            cbar_std.set_label(f"std ({units})", fontsize=8)
            cbar_std.ax.tick_params(labelsize=7)

    fig.suptitle(f"Case study — {timestamp_label} — "
                 f"{cfg.get('region_name', '')}",
                 fontsize=14, fontweight="bold")
    caption = (
        f"Single test-set time step ({timestamp_label}, index {t_idx}) "
        f"showing LR view, reference, and per-method ensemble mean for all "
        f"output variables. The last column shows the ensemble standard "
        f"deviation of the last listed method."
    )
    _save_figure(fig, cfg, "case_study_panel", caption,
                 statistics=["other"], plot_types=["map"])


# ---------------------------------------------------------------------------
# 10. Absolute Quantile-MAE table
# ---------------------------------------------------------------------------


def write_absolute_qmae_table(cfg):
    """Aggregate quantile_mae_table.csv files (ancestors) into a single CSV
    of absolute values, to accompany the relative-performance heatmaps.
    """
    rows = []
    for ancestor_dir in cfg.get("input_files", []):
        for csv_path in Path(ancestor_dir).rglob("quantile_mae_table.csv"):
            df = pd.read_csv(csv_path)
            df["source"] = str(csv_path.parent.name)
            rows.append(df)
    if not rows:
        logger.warning("No quantile_mae_table.csv found in ancestors.")
        return
    full = pd.concat(rows, ignore_index=True)
    out_csv = os.path.join(cfg["work_dir"], "quantile_mae_absolute.csv")
    full.to_csv(out_csv, index=False)
    logger.info("Saved absolute Quantile-MAE table: %s", out_csv)


# ---------------------------------------------------------------------------
# 11. Per-pixel Quantile-MAE spatial map
# ---------------------------------------------------------------------------


def _per_pixel_quantile_mae(truth, pred_ens, quantiles):
    """Compute per-pixel quantile MAE between truth and ensemble-mean."""
    ens_mean = np.mean(pred_ens, axis=0)            # (T, X, Y)
    q_truth = np.nanquantile(truth, quantiles, axis=0)    # (Q, X, Y)
    q_pred = np.nanquantile(ens_mean, quantiles, axis=0)  # (Q, X, Y)
    return np.nanmean(np.abs(q_truth - q_pred), axis=0)   # (X, Y)


def plot_spatial_quantile_mae(truth, preds_methods, method_names, var_name,
                              extent, cfg):
    """Per-pixel Quantile-MAE map (replaces Fig 10 of the submitted MS).

    Computed over the upper-tail quantiles for precipitation
    (default 0.95--0.99) or full-distribution quantiles for other variables.
    Plots one map per method + a ratio map (method2 / method1) at the end.
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    if var_name == "pr":
        quantiles = np.linspace(0.95, 0.99, 5)
        title_q = "Q[0.95,0.99]"
    else:
        quantiles = np.linspace(0.01, 0.99, 99)
        title_q = "Q[0.01,0.99]"

    units = UNITS.get(var_name, "")
    region = cfg.get("region_name", "")
    qmae_maps = [
        _per_pixel_quantile_mae(truth, p, quantiles) for p in preds_methods
    ]
    vmax = float(np.nanmax([m.max() for m in qmae_maps]))

    n_cols = len(method_names) + (1 if len(method_names) >= 2 else 0)
    fig = plt.figure(figsize=(4.0 * n_cols, 3.6))
    projection = ccrs.PlateCarree()

    last_im = None
    for col, (method, qmap) in enumerate(zip(method_names, qmae_maps)):
        ax = fig.add_subplot(1, n_cols, col + 1, projection=projection)
        last_im = ax.imshow(qmap, cmap="magma_r", vmin=0.0, vmax=vmax,
                            origin="lower", extent=extent,
                            transform=projection, aspect="auto")
        ax.coastlines(resolution="50m", linewidth=0.8, color="black")
        ax.add_feature(cfeature.BORDERS, linewidth=0.5, alpha=0.5)
        qmae_label = _ref_metric_label(qmap)
        ax.set_title(f"{method}\nQuantile MAE ({units})  {qmae_label}",
                     fontsize=11, fontweight="bold")

    fig.colorbar(last_im, ax=fig.get_axes()[:len(method_names)],
                 orientation="horizontal", fraction=0.04, pad=0.10, aspect=40,
                 label=f"Quantile MAE — {title_q} ({units})")

    # Ratio map (method2 / method1) if there are at least two methods.
    if len(method_names) >= 2:
        from scipy.ndimage import gaussian_filter

        method1, method2 = method_names[0], method_names[1]
        ratio = qmae_maps[1] / np.where(qmae_maps[0] > 0,
                                        qmae_maps[0], np.nan)
        # NaN-aware Gaussian smoothing: weighted average of finite pixels.
        # Pixel-scale Q-MAE differences are dominated by sampling noise;
        # smoothing reveals coherent regional bias between the two methods.
        sigma = float(cfg.get("ratio_smoothing_sigma", 3.0))
        if sigma > 0:
            mask = np.isfinite(ratio).astype(float)
            filled = np.where(mask > 0, ratio, 0.0)
            num = gaussian_filter(filled, sigma=sigma, mode="nearest")
            den = gaussian_filter(mask, sigma=sigma, mode="nearest")
            with np.errstate(invalid="ignore", divide="ignore"):
                ratio_smooth = np.where(den > 0, num / den, np.nan)
            ratio_to_plot = np.where(mask > 0, ratio_smooth, np.nan)
            smoothing_note = rf"smoothed $\sigma$={sigma:g} px"
        else:
            ratio_to_plot = ratio
            smoothing_note = "unsmoothed"

        ax = fig.add_subplot(1, n_cols, n_cols, projection=projection)
        ratio_norm = mcolors.TwoSlopeNorm(vmin=0.6, vcenter=1.0, vmax=1.4)
        im = ax.imshow(ratio_to_plot, cmap=DIVERGING_CMAP, norm=ratio_norm,
                       origin="lower", extent=extent, transform=projection,
                       aspect="auto")
        ax.coastlines(resolution="50m", linewidth=0.8, color="black")
        ax.add_feature(cfeature.BORDERS, linewidth=0.5, alpha=0.5)
        # Scalar summary uses the unsmoothed ratio so it is independent of
        # the cosmetic smoothing choice.
        finite_ratio = ratio[np.isfinite(ratio)]
        if finite_ratio.size:
            med = float(np.median(finite_ratio))
            frac_better = float(np.mean(finite_ratio < 1.0))
            ratio_label = (rf"median={_g(med, 2)}  "
                           rf"area<1={frac_better * 100:.0f}%")
        else:
            ratio_label = "median=NaN  area<1=NaN"
        ax.set_title(f"{method2} / {method1}  ratio ({smoothing_note})\n"
                     f"{ratio_label}",
                     fontsize=11, fontweight="bold")
        fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.04,
                     pad=0.10, aspect=20, label="ratio")

    fig.suptitle(f"Spatial Quantile-MAE — {var_name} — {region}",
                 fontsize=13, fontweight="bold", y=1.08)
    plt.tight_layout()
    caption = (
        f"Per-pixel quantile mean-absolute error of {var_name} over the "
        f"{title_q} range, computed from the ensemble-mean field. "
        f"Maps are shown for each method; if two or more methods are "
        f"provided, the final panel shows the ratio of the second to the "
        f"first."
    )
    _save_figure(fig, cfg, f"spatial_quantile_mae_{var_name}", caption,
                 statistics=["other"], plot_types=["map"])


# ---------------------------------------------------------------------------
# 12. Combined Quantile-MAE summary table (replaces Figs 5/8/10)
# ---------------------------------------------------------------------------


def write_combined_quantile_mae_table(cfg):
    """Combine the per-region absolute Quantile-MAE values into one table.

    Walks the ancestor directories for ``quantile_mae_table.csv`` files
    (produced by ``ml_downscaling_evaluation.py`` with analysis_type =
    "quantile_MAE"). Each file is expected to contain one column per
    method plus a ``variable`` and ``metric`` column. The region for each
    file is inferred from the parent directory name (matching strings
    "CE", "IBE", "SCA" or the suffixes "Iberia"/"Scandinavia").

    Output: a wide-format CSV with columns
    ``[region, composite, method, abs_qmae, rel_to_<reference>]`` and a
    pretty-printed LaTeX table.
    """
    rel_ref = cfg.get("relative_reference", "AFM-baseline")
    rows = []

    def _infer_region(path):
        s = str(path).lower()
        if "scandin" in s or "sca" in s:
            return "SCA"
        if "iberia" in s or "ibe" in s:
            return "IBE"
        return "CE"

    for ancestor_dir in cfg.get("input_files", []):
        for csv_path in Path(ancestor_dir).rglob("quantile_mae_table.csv"):
            df = pd.read_csv(csv_path)
            region = _infer_region(csv_path)
            method_cols = [c for c in df.columns
                           if c not in ("variable", "metric", "source")]
            for _, r in df.iterrows():
                composite = r["variable"]
                for m in method_cols:
                    try:
                        val = float(r[m])
                    except (TypeError, ValueError):
                        val = float("nan")
                    rows.append({
                        "region": region,
                        "composite": composite,
                        "method": m,
                        "abs_qmae": val,
                    })
    if not rows:
        logger.warning("No quantile_mae_table.csv found in ancestors.")
        return
    long_df = pd.DataFrame(rows)

    # Relative to reference method.
    ref_lookup = long_df[long_df["method"] == rel_ref].set_index(
        ["region", "composite"])["abs_qmae"].to_dict()
    long_df["rel_to_" + rel_ref] = long_df.apply(
        lambda r: r["abs_qmae"] / ref_lookup.get(
            (r["region"], r["composite"]), float("nan")),
        axis=1,
    )

    out_csv = os.path.join(cfg["work_dir"], "quantile_mae_combined.csv")
    long_df.to_csv(out_csv, index=False)
    logger.info("Saved combined Q-MAE table: %s", out_csv)

    # Pretty wide-format table: region × composite rows, method columns.
    wide_abs = long_df.pivot_table(
        index=["region", "composite"], columns="method", values="abs_qmae"
    )
    wide_rel = long_df.pivot_table(
        index=["region", "composite"], columns="method",
        values="rel_to_" + rel_ref,
    )
    out_abs = os.path.join(cfg["work_dir"], "quantile_mae_absolute_wide.csv")
    out_rel = os.path.join(cfg["work_dir"], "quantile_mae_relative_wide.csv")
    wide_abs.to_csv(out_abs)
    wide_rel.to_csv(out_rel)
    logger.info("Saved wide tables: %s, %s", out_abs, out_rel)

    # LaTeX table for direct insertion into the manuscript.
    out_tex = os.path.join(cfg["work_dir"], "quantile_mae_combined.tex")
    with open(out_tex, "w") as fh:
        fh.write("% Auto-generated by ml_downscaling_extras.py\n")
        fh.write("\\begin{tabular}{ll" +
                 "r" * len(wide_abs.columns) + "}\n\\hline\n")
        fh.write("Region & Composite & " +
                 " & ".join(wide_abs.columns.astype(str)) + " \\\\\n\\hline\n")
        for (region, comp), row in wide_abs.iterrows():
            fh.write(f"{region} & {comp} & " +
                     " & ".join(f"{v:.3f}" for v in row.values) +
                     " \\\\\n")
        fh.write("\\hline\n\\end{tabular}\n")
    logger.info("Saved LaTeX table: %s", out_tex)


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------


DISPATCH = {
    "climatology": "per_variable",
    "percentile_map": "per_variable",
    "cdd": "per_variable",
    "rx1day": "per_variable",
    "conditional_rank_histogram": "per_variable",
    "return_level": "per_variable",
    "temporal_spectrum": "per_variable",
    "spatial_quantile_mae": "per_variable",
    "cross_variable_correlation": "all_variables",
    "case_study": "all_variables",
    "quantile_mae_table": "no_data",
    "quantile_mae_summary_table": "no_data",
}


def main(cfg):
    logger.setLevel(cfg.get("log_level", "INFO").upper())
    analysis_type = cfg.get("analysis_type", "climatology")
    ml_methods = cfg.get("ml_methods", [])
    reference_name = cfg.get("reference", "HIGHRES-REF")
    logger.info("ml_downscaling_extras: analysis_type=%s methods=%s",
                analysis_type, ml_methods)

    if DISPATCH.get(analysis_type) == "no_data":
        if analysis_type == "quantile_mae_table":
            write_absolute_qmae_table(cfg)
        elif analysis_type == "quantile_mae_summary_table":
            write_combined_quantile_mae_table(cfg)
        return

    input_data = cfg["input_data"].values()
    grouped_data = group_metadata(input_data, "short_name", sort="dataset")

    if DISPATCH.get(analysis_type) == "all_variables":
        # Single-call dispatch
        # Extent extracted inside the function from any one variable
        extent = None
        if analysis_type == "cross_variable_correlation":
            plot_cross_variable_correlation(grouped_data, reference_name,
                                            ml_methods, cfg, extent)
        elif analysis_type == "case_study":
            plot_case_study(grouped_data, reference_name, ml_methods, cfg)
        return

    # Per-variable dispatch
    for var_name in grouped_data:
        if analysis_type in {"cdd", "rx1day", "return_level",
                             "conditional_rank_histogram"} \
                and var_name != "pr":
            continue
        if analysis_type == "percentile_map" \
                and var_name not in {"pr", "tas"}:
            continue
        if analysis_type == "temporal_spectrum" \
                and var_name not in {"pr", "tas"}:
            continue
        truth_cube, truth_array, preds_methods, method_names = \
            _load_method_truth(grouped_data, var_name, reference_name,
                               ml_methods)
        if truth_cube is None or not preds_methods:
            continue
        ref_datasets = [d for d in grouped_data[var_name]
                        if d["dataset"] == reference_name]
        update_variable_units(ref_datasets)
        extent = _extract_extent(truth_cube)

        if analysis_type == "climatology":
            plot_climatology(truth_array, preds_methods, method_names,
                             var_name, extent, cfg)
        elif analysis_type == "percentile_map":
            plot_percentile_map(truth_array, preds_methods, method_names,
                                var_name, extent, cfg,
                                percentile=cfg.get("percentile", 99.0))
        elif analysis_type == "cdd":
            plot_cdd(truth_array, preds_methods, method_names, var_name,
                     extent, cfg)
        elif analysis_type == "rx1day":
            plot_rx1day(truth_array, preds_methods, method_names, var_name,
                        extent, cfg)
        elif analysis_type == "conditional_rank_histogram":
            plot_conditional_rank_histogram(truth_array, preds_methods,
                                            method_names, var_name, cfg)
        elif analysis_type == "return_level":
            plot_return_level(truth_array, preds_methods, method_names,
                              var_name, extent, cfg)
        elif analysis_type == "temporal_spectrum":
            plot_temporal_spectrum(truth_array, preds_methods, method_names,
                                   var_name, cfg)
        elif analysis_type == "spatial_quantile_mae":
            plot_spatial_quantile_mae(truth_array, preds_methods, method_names,
                                      var_name, extent, cfg)
        else:
            logger.warning("Unknown analysis_type: %s", analysis_type)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
