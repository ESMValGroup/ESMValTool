from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr


def _style():
    plt.rcParams.update(
        {
            "figure.dpi": 110,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "font.size": 10,
        }
    )


def _color_limits(values: np.ndarray, *, diverging: bool) -> dict:
    finite = np.asarray(values, dtype=np.float64)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return {}
    if diverging:
        vmax = float(np.nanpercentile(np.abs(finite), 98)) or 1.0
        return {"vmin": -vmax, "vmax": vmax}
    lo, hi = np.nanpercentile(finite, [2, 98])
    if not np.isfinite(lo) or hi <= lo:
        return {}
    return {"vmin": float(lo), "vmax": float(hi)}


def plot_cube_maps(
    ds: xr.Dataset,
    time_index: int = 0,
    cmap: str = "RdYlBu_r",
    title: str | None = None,
    reduce: str = "none",
):
    """One spatial map per variable.

    ``reduce='mean'`` averages over time (needed for SHAP cubes, which are
    filled only on test samples). Colour limits use the 2–98 percentiles so
    small CMOR units remain visible.
    """
    _style()
    names = list(ds.data_vars)
    n = len(names)
    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 3.6), squeeze=False)
    diverging = cmap.lower() in {"rdbu_r", "rdbu", "coolwarm", "seismic"}
    stamp = ""
    if "time" in ds.dims and reduce != "mean":
        stamp = str(pd.to_datetime(ds.time.values[time_index]))
    elif reduce == "mean":
        stamp = "time mean"
    for ax, name in zip(axes[0], names):
        da = ds[name]
        if "time" in da.dims:
            da = da.mean("time") if reduce == "mean" else da.isel(time=time_index)
        im = ax.pcolormesh(
            ds.lon,
            ds.lat,
            da,
            shading="auto",
            cmap=cmap,
            **_color_limits(da.values, diverging=diverging),
        )
        ax.set_xlabel("lon")
        ax.set_ylabel("lat")
        ax.set_title(name)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(title or f"spatial maps ({stamp})")
    fig.tight_layout()
    return fig


def plot_spatial_mean_timeseries(ds: xr.Dataset, title: str = "domain-mean time series"):
    _style()
    fig, ax = plt.subplots(figsize=(8, 3.5))
    for name, da in ds.data_vars.items():
        series = da.mean(("lat", "lon"))
        ax.plot(ds.time, series, label=name)
    ax.set_xlabel("time")
    ax.set_ylabel("value")
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    return fig


def plot_histograms(table_like: pd.DataFrame, columns: Sequence[str], title: str = "feature distributions"):
    _style()
    cols = list(columns)
    n = len(cols)
    fig, axes = plt.subplots(1, n, figsize=(3.2 * n, 3.2), squeeze=False)
    for ax, col in zip(axes[0], cols):
        ax.hist(table_like[col].dropna(), bins=25, color="steelblue", edgecolor="white")
        ax.set_title(col)
        ax.set_xlabel(col)
        ax.set_ylabel("count")
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_correlation(x: np.ndarray, names: Sequence[str], title: str = "feature correlation"):
    _style()
    corr = pd.DataFrame(x, columns=list(names)).corr()
    fig, ax = plt.subplots(figsize=(5.5, 4.8))
    im = ax.imshow(corr, cmap="coolwarm", vmin=-1, vmax=1)
    ax.set_xticks(range(len(names)))
    ax.set_yticks(range(len(names)))
    ax.set_xticklabels(names, rotation=45, ha="right")
    ax.set_yticklabels(names)
    for i in range(len(names)):
        for j in range(len(names)):
            ax.text(j, i, f"{corr.iloc[i, j]:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.046)
    ax.set_title(title)
    fig.tight_layout()
    return fig


def plot_split_coverage(train_index: pd.MultiIndex, test_index: pd.MultiIndex, title: str = "train / test coverage"):
    _style()
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
    for ax, index, label, color in (
        (axes[0], train_index, "train", "tab:blue"),
        (axes[1], test_index, "test", "tab:orange"),
    ):
        ax.scatter(index.get_level_values("lon"), index.get_level_values("lat"), s=8, alpha=0.35, c=color)
        ax.set_xlabel("lon")
        ax.set_ylabel("lat")
        ax.set_title(f"{label} n={len(index)}")
        ax.set_aspect("equal", adjustable="box")
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_split_time(train_index: pd.MultiIndex, test_index: pd.MultiIndex):
    _style()
    fig, ax = plt.subplots(figsize=(8, 2.8))
    tr = pd.to_datetime(np.unique(train_index.get_level_values("time")))
    te = pd.to_datetime(np.unique(test_index.get_level_values("time")))
    ax.scatter(tr, np.zeros(len(tr)), marker="|", s=400, label="train")
    ax.scatter(te, np.ones(len(te)), marker="|", s=400, label="test")
    ax.set_yticks([0, 1], ["train", "test"])
    ax.set_xlabel("time")
    ax.set_title("time split")
    ax.legend(loc="upper left")
    fig.tight_layout()
    return fig


def plot_sampling_map(all_index: pd.MultiIndex, sampled_index: pd.MultiIndex, title: str = "training sample subset"):
    _style()
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    ax.scatter(all_index.get_level_values("lon"), all_index.get_level_values("lat"), s=6, alpha=0.15, c="gray", label="all train")
    ax.scatter(sampled_index.get_level_values("lon"), sampled_index.get_level_values("lat"), s=12, alpha=0.6, c="crimson", label="sampled")
    ax.set_xlabel("lon")
    ax.set_ylabel("lat")
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    return fig


def plot_pred_scatter(y_true, y_pred, metrics: dict | None = None):
    _style()
    yt = np.asarray(y_true).ravel()
    yp = np.asarray(y_pred).ravel()
    fig, ax = plt.subplots(figsize=(4.8, 4.8))
    ax.scatter(yt, yp, s=10, alpha=0.4)
    lims = [np.nanmin([yt, yp]), np.nanmax([yt, yp])]
    ax.plot(lims, lims, "k--", lw=1)
    ax.set_xlabel("observed Y")
    ax.set_ylabel("predicted Y")
    title = "predicted vs observed"
    if metrics:
        title += f"\nR$^2$={metrics['r2']:.3f}  RMSE={metrics['rmse']:.3f}  ubRMSE={metrics['ubrmse']:.3f}"
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    return fig


def plot_error_map(obs: xr.DataArray, pred: xr.DataArray, time_index: int | None = None):
    _style()
    if time_index is None:
        obs_m = obs.mean("time")
        pred_m = pred.mean("time")
        subtitle = "time mean"
    else:
        obs_m = obs.isel(time=time_index)
        pred_m = pred.isel(time=time_index)
        subtitle = str(pd.to_datetime(obs.time.values[time_index]))
    err = pred_m - obs_m
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    for ax, da, title, cmap in (
        (axes[0], obs_m, "observed", "viridis"),
        (axes[1], pred_m, "predicted", "viridis"),
        (axes[2], err, "pred − obs", "RdBu_r"),
    ):
        kw = dict(
            shading="auto",
            cmap=cmap,
            **_color_limits(da.values, diverging="obs" in title and title.startswith("pred")),
        )
        im = ax.pcolormesh(da.lon, da.lat, da, **kw)
        ax.set_title(title)
        ax.set_xlabel("lon")
        ax.set_ylabel("lat")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(f"prediction maps ({subtitle})")
    fig.tight_layout()
    return fig


def plot_importance_bars(frame: pd.DataFrame, value_col: str = "importance", title: str = "feature importance"):
    _style()
    df = frame.sort_values(value_col)
    fig, ax = plt.subplots(figsize=(6, max(2.5, 0.35 * len(df) + 1)))
    ax.barh(df["feature"], df[value_col], color="steelblue", xerr=df["std"] if "std" in df.columns else None)
    ax.set_xlabel(value_col)
    ax.set_title(title)
    fig.tight_layout()
    return fig


def plot_shap_summary(shap_values: np.ndarray, x: np.ndarray, names: Sequence[str], title: str = "SHAP summary"):
    _style()
    names = list(names)
    order = np.argsort(np.abs(shap_values).mean(axis=0))
    fig, ax = plt.subplots(figsize=(6.5, max(2.8, 0.4 * len(names) + 1)))
    rng = np.random.default_rng(0)
    for row, j in enumerate(order):
        vals = shap_values[:, j]
        feat = x[:, j]
        colors = (feat - np.nanmin(feat)) / (np.nanmax(feat) - np.nanmin(feat) + 1e-12)
        jitter = rng.normal(0, 0.08, size=len(vals))
        sc = ax.scatter(vals, np.full(len(vals), row) + jitter, c=colors, cmap="coolwarm", s=8, alpha=0.7, vmin=0, vmax=1)
    ax.set_yticks(range(len(order)), [names[j] for j in order])
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("SHAP value")
    ax.set_title(title)
    fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.02, label="feature value (low→high)")
    fig.tight_layout()
    return fig


def plot_shap_maps(
    shap_cube: xr.Dataset,
    time_index: int | None = None,
    title: str = "SHAP restored to cube",
):
    """Map SHAP that has been written back onto ``(time, lat, lon)``.

    Default is the time mean: a chronological split leaves early times empty.
    """
    if time_index is None:
        return plot_cube_maps(shap_cube, reduce="mean", cmap="RdBu_r", title=title)
    return plot_cube_maps(
        shap_cube,
        time_index=min(time_index, shap_cube.sizes["time"] - 1),
        cmap="RdBu_r",
        title=title,
    )


def plot_comparison_bars(
    frame: pd.DataFrame,
    value_col: str,
    title: str,
    xlabel: str | None = None,
):
    """Bar chart of a comparison table (distance or slope difference)."""
    _style()
    df = frame.copy()
    if "feature" not in df.columns:
        df = df.reset_index().rename(columns={"index": "feature"})
    if "interacts_with" in df.columns:
        df["feature"] = df["feature"].astype(str) + " × " + df["interacts_with"].astype(str)
    df = df.dropna(subset=[value_col]).sort_values(value_col)
    fig, ax = plt.subplots(figsize=(6.2, max(2.6, 0.4 * len(df) + 1)))
    colors = ["tab:red" if v < 0 else "steelblue" for v in df[value_col]]
    ax.barh(df["feature"].astype(str), df[value_col], color=colors)
    ax.axvline(0, color="k", lw=0.7)
    ax.set_xlabel(xlabel or value_col)
    ax.set_title(title)
    fig.tight_layout()
    return fig


def plot_ale(ale: dict, title: str = "Accumulated Local Effects"):
    _style()
    names = list(ale)
    n = len(names)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 3.3), squeeze=False)
    for ax, name in zip(axes[0], names):
        df = ale[name]
        ax.plot(df["x"], df["ale"], color="tab:green")
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xlabel(name)
        ax.set_ylabel("ALE")
        ax.set_title(name)
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_ice_pdp(ice_pdp: dict, title: str = "ICE (thin) + PDP (thick)"):
    _style()
    names = list(ice_pdp)
    n = len(names)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 3.4), squeeze=False)
    for ax, name in zip(axes[0], names):
        pack = ice_pdp[name]
        grid, ice, pdp = pack["grid"], pack["ice"], pack["pdp"]
        for row in ice:
            ax.plot(grid, row, color="0.75", lw=0.7, alpha=0.7)
        ax.plot(grid, pdp, color="tab:red", lw=2.2, label="PDP")
        ax.set_xlabel(name)
        ax.set_ylabel("prediction")
        ax.set_title(name)
        ax.legend()
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def _subsample_rows(n: int, max_points: int, seed: int = 0) -> np.ndarray:
    if n <= max_points:
        return np.arange(n)
    return np.random.default_rng(seed).choice(n, size=max_points, replace=False)


def plot_shap_dependence(
    x,
    shap_values,
    names: Sequence[str],
    title: str = "SHAP dependence",
    n_bins: int = 20,
    max_points: int = 2500,
):
    """One panel per driver: x-axis is the variable, y-axis is its SHAP value."""
    from rf_shap.explain.evaluate import binned_dependence

    _style()
    features = np.asarray(x, dtype=np.float64)
    values = np.asarray(shap_values, dtype=np.float64)
    names = list(names)
    n = len(names)
    fig, axes = plt.subplots(1, n, figsize=(3.6 * n, 3.5), squeeze=False)
    keep = _subsample_rows(len(features), max_points)
    for ax, j, name in zip(axes[0], range(n), names):
        ax.scatter(features[keep, j], values[keep, j], s=6, alpha=0.25, c="0.45", linewidths=0)
        curve = binned_dependence(features[:, j], values[:, j], n_bins=n_bins)
        if len(curve):
            ax.plot(curve["x"], curve["y"], color="tab:red", lw=2.0, label="binned mean")
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xlabel(name)
        ax.set_ylabel(f"SHAP of {name}")
        ax.set_title(name)
    axes[0, 0].legend(loc="best", fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_shap_dependence_overlay(
    x_obs,
    shap_obs,
    names_obs: Sequence[str],
    x_cmip,
    shap_cmip,
    names_cmip: Sequence[str],
    title: str = "SHAP dependence: observations vs CMIP",
    obs_label: str = "observations",
    cmip_label: str = "CMIP",
    n_bins: int = 20,
):
    """Binned SHAP dependence curves for shared drivers."""
    from rf_shap.explain.evaluate import binned_dependence

    _style()
    names = [name for name in names_obs if name in set(names_cmip)]
    if not names:
        raise ValueError("No shared features for SHAP dependence overlay")
    obs_x = np.asarray(x_obs, dtype=np.float64)
    obs_s = np.asarray(shap_obs, dtype=np.float64)
    cmip_x = np.asarray(x_cmip, dtype=np.float64)
    cmip_s = np.asarray(shap_cmip, dtype=np.float64)
    obs_idx = {name: i for i, name in enumerate(names_obs)}
    cmip_idx = {name: i for i, name in enumerate(names_cmip)}
    fig, axes = plt.subplots(1, len(names), figsize=(3.6 * len(names), 3.5), squeeze=False)
    for ax, name in zip(axes[0], names):
        i, j = obs_idx[name], cmip_idx[name]
        left = binned_dependence(obs_x[:, i], obs_s[:, i], n_bins=n_bins)
        right = binned_dependence(cmip_x[:, j], cmip_s[:, j], n_bins=n_bins)
        if len(left):
            ax.plot(left["x"], left["y"], color="tab:blue", lw=2.2, label=obs_label)
        if len(right):
            ax.plot(right["x"], right["y"], color="tab:orange", lw=2.2, label=cmip_label)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xlabel(name)
        ax.set_ylabel(f"SHAP of {name}")
        ax.set_title(name)
    axes[0, 0].legend(loc="best", fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    return fig


_VERIFY_COLORS = {
    "true": "0.15",
    "m1": "tab:blue",
    "m2": "tab:red",
    "m3": "tab:green",
}
_VERIFY_LABELS = {
    "true": "true",
    "m1": "m1 linear",
    "m2": "m2 sign flip",
    "m3": "m3 weak curve",
}


def plot_shap_metric_verification(result: dict, title: str = "Fig. 2  SHAP metrics against a nonlinear truth"):
    """Three-panel figure: x–y maps, SHAP dependence, m1–m3 versus truth."""
    _style()
    metrics = result["metrics"]
    curves = result["curves"]
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.0))

    ax = axes[0]
    scatter = curves["true:scatter"]
    keep = _subsample_rows(len(scatter), 800, seed=0)
    ax.scatter(scatter["x"].to_numpy()[keep], scatter["y"].to_numpy()[keep], s=6, alpha=0.15, c="0.65", linewidths=0)
    for case in ("true", "m1", "m2", "m3"):
        curve = curves[f"{case}:y"]
        ax.plot(curve["x"], curve["y"], color=_VERIFY_COLORS[case], lw=2.2, label=_VERIFY_LABELS[case])
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$y$")
    ax.set_title("(a) $x$–$y$ relationship")
    ax.legend(fontsize=7.5, loc="best")

    ax = axes[1]
    for case in ("true", "m1", "m2", "m3"):
        curve = curves[f"{case}:x1"]
        ax.plot(curve["x"], curve["y"], color=_VERIFY_COLORS[case], lw=2.2, label=_VERIFY_LABELS[case])
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("SHAP of $x_1$")
    ax.set_title("(b) SHAP dependence")
    ax.legend(fontsize=7.5, loc="best")

    ax = axes[2]
    cases = ["m1", "m2", "m3"]
    x1_raw = metrics[metrics["feature"] == "x1"]
    xpos = np.arange(len(cases))
    width = 0.36

    def _median_band(case: str, column: str, absolute: bool = False) -> tuple[float, float, float]:
        values = x1_raw.loc[x1_raw["case"] == case, column].to_numpy(dtype=float)
        if absolute:
            values = np.abs(values)
        med = float(np.median(values))
        lo = med - float(np.quantile(values, 0.1))
        hi = float(np.quantile(values, 0.9)) - med
        return med, max(lo, 0.0), max(hi, 0.0)

    dist = [_median_band(c, "dependence_distance") for c in cases]
    slope = [_median_band(c, "slope_difference", absolute=True) for c in cases]
    ax.bar(
        xpos - width / 2,
        [d[0] for d in dist],
        width,
        yerr=np.vstack([[d[1] for d in dist], [d[2] for d in dist]]),
        color="0.35",
        label="curve distance",
        capsize=3,
    )
    ax.bar(
        xpos + width / 2,
        [s[0] for s in slope],
        width,
        yerr=np.vstack([[s[1] for s in slope], [s[2] for s in slope]]),
        color="tab:red",
        alpha=0.85,
        label="|slope difference|",
        capsize=3,
    )
    ax.set_xticks(xpos)
    ax.set_xticklabels([_VERIFY_LABELS[c] for c in cases])
    ax.set_ylabel("median over repeats")
    ax.set_title("(c) m1–m3 versus true SHAP")
    ax.legend(fontsize=8, loc="upper right")

    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    return fig


def plot_ale_metric_supplement(result: dict, title: str = "Supplement  ALE distances versus the nonlinear truth"):
    """ALE curves and distances for m1–m3 (not in Fig. 2)."""
    _style()
    summary = result["summary"]
    curves = result["curves"]
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.0))
    ax = axes[0]
    for case in ("true", "m1", "m2", "m3"):
        curve = curves[f"{case}:ale:x1"]
        ax.plot(curve["x"], curve["ale"], color=_VERIFY_COLORS[case], lw=2.0, label=_VERIFY_LABELS[case])
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("ALE of $x_1$")
    ax.set_title("(a) ALE curves")
    ax.legend(fontsize=7.5)

    ax = axes[1]
    cases = ["m1", "m2", "m3"]
    x1 = summary[(summary["feature"] == "x1")].set_index("case")
    xpos = np.arange(len(cases))
    dist = [float(x1.loc[c, "ale_distance"]) for c in cases]
    lo = [dist[i] - float(x1.loc[c, "ale_q10"]) for i, c in enumerate(cases)]
    hi = [float(x1.loc[c, "ale_q90"]) - dist[i] for i, c in enumerate(cases)]
    ax.bar(xpos, dist, yerr=np.vstack([lo, hi]), color="0.45", capsize=3)
    ax.set_xticks(xpos)
    ax.set_xticklabels([_VERIFY_LABELS[c] for c in cases])
    ax.set_ylabel("median ALE distance")
    ax.set_title("(b) ALE curve distance for $x_1$")
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    return fig


def plot_ale_overlay(
    ale_obs: dict,
    ale_cmip: dict,
    title: str = "ALE: observations vs CMIP",
    obs_label: str = "observations",
    cmip_label: str = "CMIP",
):
    _style()
    names = [name for name in ale_obs if name in ale_cmip]
    if not names:
        raise ValueError("No shared features for ALE overlay")
    fig, axes = plt.subplots(1, len(names), figsize=(3.8 * len(names), 3.4), squeeze=False)
    for ax, name in zip(axes[0], names):
        ax.plot(ale_obs[name]["x"], ale_obs[name]["ale"], color="tab:blue", lw=2.0, label=obs_label)
        ax.plot(ale_cmip[name]["x"], ale_cmip[name]["ale"], color="tab:orange", lw=2.0, label=cmip_label)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xlabel(name)
        ax.set_ylabel("ALE")
        ax.set_title(name)
    axes[0, 0].legend(loc="best", fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_pdp_overlay(
    ice_obs: dict,
    ice_cmip: dict,
    title: str = "PDP: observations vs CMIP",
    obs_label: str = "observations",
    cmip_label: str = "CMIP",
):
    _style()
    names = [name for name in ice_obs if name in ice_cmip]
    if not names:
        raise ValueError("No shared features for PDP overlay")
    fig, axes = plt.subplots(1, len(names), figsize=(3.8 * len(names), 3.4), squeeze=False)
    for ax, name in zip(axes[0], names):
        ax.plot(ice_obs[name]["grid"], ice_obs[name]["pdp"], color="tab:blue", lw=2.2, label=obs_label)
        ax.plot(ice_cmip[name]["grid"], ice_cmip[name]["pdp"], color="tab:orange", lw=2.2, label=cmip_label)
        ax.set_xlabel(name)
        ax.set_ylabel("prediction")
        ax.set_title(name)
    axes[0, 0].legend(loc="best", fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_importance_overlay(
    frame_obs: pd.DataFrame,
    frame_cmip: pd.DataFrame,
    title: str = "Importance: observations vs CMIP",
    obs_label: str = "observations",
    cmip_label: str = "CMIP",
    value_col: str = "importance",
):
    _style()
    left = frame_obs.set_index("feature")[value_col]
    right = frame_cmip.set_index("feature")[value_col]
    names = [name for name in left.index if name in right.index]
    y = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(6.5, max(2.6, 0.4 * len(names) + 1)))
    ax.barh(y - 0.18, [left.loc[n] for n in names], height=0.35, color="tab:blue", label=obs_label)
    ax.barh(y + 0.18, [right.loc[n] for n in names], height=0.35, color="tab:orange", label=cmip_label)
    ax.set_yticks(y, names)
    ax.set_xlabel(value_col)
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    return fig


def plot_metrics_compare(
    metrics_obs: dict,
    metrics_cmip: dict,
    title: str = "Mapping skill: observations vs CMIP",
    obs_label: str = "observations",
    cmip_label: str = "CMIP",
):
    _style()
    keys = [k for k in ("r2", "rmse", "ubrmse", "bias") if k in metrics_obs and k in metrics_cmip]
    x = np.arange(len(keys))
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.bar(x - 0.18, [metrics_obs[k] for k in keys], width=0.35, color="tab:blue", label=obs_label)
    ax.bar(x + 0.18, [metrics_cmip[k] for k in keys], width=0.35, color="tab:orange", label=cmip_label)
    ax.set_xticks(x, keys)
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    return fig


def plot_cube_difference(
    obs_ds: xr.Dataset,
    cmip_ds: xr.Dataset,
    title: str = "CMIP − observations (time mean)",
    cmap: str = "RdBu_r",
):
    """Time-mean difference maps for shared variables."""
    _style()
    names = [name for name in obs_ds.data_vars if name in cmip_ds.data_vars]
    if not names:
        raise ValueError("No shared variables for cube difference")
    fig, axes = plt.subplots(1, len(names), figsize=(4.2 * len(names), 3.6), squeeze=False)
    for ax, name in zip(axes[0], names):
        left = obs_ds[name]
        right = cmip_ds[name]
        if "time" in left.dims:
            left = left.mean("time")
        if "time" in right.dims:
            right = right.mean("time")
        left, right = xr.align(left, right, join="inner")
        diff = right - left
        im = ax.pcolormesh(
            diff.lon,
            diff.lat,
            diff,
            shading="auto",
            cmap=cmap,
            **_color_limits(diff.values, diverging=True),
        )
        ax.set_title(name)
        ax.set_xlabel("lon")
        ax.set_ylabel("lat")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def write_three_column_html(
    path: str | Path,
    title: str,
    obs_label: str,
    cmip_label: str,
    rows: Sequence[dict],
    intro: str = "",
) -> Path:
    """Write a three-column gallery: observations | CMIP | comparison."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    root = path.parent

    def _src(value) -> str:
        if not value:
            return ""
        target = Path(value)
        try:
            return target.resolve().relative_to(root.resolve()).as_posix()
        except ValueError:
            return target.as_posix()

    blocks = []
    for row in rows:
        cells = []
        for key in ("obs", "cmip", "compare"):
            src = _src(row.get(key))
            if src:
                cells.append(
                    f'<div class="cell"><img src="{src}" alt="{row.get("title", "")}"></div>'
                )
            else:
                cells.append('<div class="cell empty">—</div>')
        caption = row.get("caption") or ""
        blocks.append(
            f'<section class="row-block">'
            f'<h2>{row.get("title", "")}</h2>'
            f'<div class="grid">{"".join(cells)}</div>'
            f'<p class="caption">{caption}</p>'
            f"</section>"
        )
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title}</title>
  <style>
    body {{ font-family: "Helvetica Neue", Helvetica, Arial, sans-serif; color: #222;
           margin: 0 auto; max-width: 1400px; padding: 1.5rem 1.2rem 4rem; background: #f4f6f8; }}
    h1 {{ font-size: 1.6rem; margin: 0 0 0.4rem; }}
    .intro {{ color: #555; margin-bottom: 1.2rem; }}
    .headers {{ display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 10px; position: sticky;
                top: 0; z-index: 2; }}
    .headers div {{ background: #1f4e79; color: #fff; text-align: center; padding: 0.55rem;
                    font-weight: 600; border-radius: 6px; }}
    .row-block {{ background: #fff; border-radius: 8px; padding: 0.8rem 0.9rem 1rem;
                  margin: 1rem 0; box-shadow: 0 1px 3px rgba(0,0,0,0.08); }}
    .row-block h2 {{ font-size: 1.05rem; margin: 0 0 0.6rem; color: #1f4e79; }}
    .grid {{ display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 10px; align-items: start; }}
    .cell img {{ width: 100%; height: auto; background: #fff; border: 1px solid #e2e6ea; border-radius: 4px; }}
    .cell.empty {{ min-height: 4rem; display: flex; align-items: center; justify-content: center;
                   color: #aaa; border: 1px dashed #ccc; border-radius: 4px; }}
    .caption {{ color: #555; font-size: 0.92rem; margin: 0.6rem 0 0; }}
  </style>
</head>
<body>
  <h1>{title}</h1>
  <p class="intro">{intro}</p>
  <div class="headers">
    <div>Observations<br><span style="font-weight:400;font-size:0.9rem">{obs_label}</span></div>
    <div>CMIP<br><span style="font-weight:400;font-size:0.9rem">{cmip_label}</span></div>
    <div>Observations vs CMIP</div>
  </div>
  {"".join(blocks)}
</body>
</html>
"""
    path.write_text(html, encoding="utf-8")
    return path
