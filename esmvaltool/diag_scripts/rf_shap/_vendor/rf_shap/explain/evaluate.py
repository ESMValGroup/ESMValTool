"""Compare two explanations: curve distance and dependence-SHAP sensitivity."""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd


def _as_1d(values) -> np.ndarray:
    return np.asarray(values, dtype=np.float64).reshape(-1)


def _sorted_curve(x, y) -> tuple[np.ndarray, np.ndarray]:
    x = _as_1d(x)
    y = _as_1d(y)
    if len(x) != len(y):
        raise ValueError("A curve needs the same number of x and y values")
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) == 0:
        return x, y
    order = np.argsort(x)
    x, y = x[order], y[order]
    unique_x, inverse = np.unique(x, return_inverse=True)
    summed = np.zeros(len(unique_x), dtype=np.float64)
    counts = np.zeros(len(unique_x), dtype=np.float64)
    np.add.at(summed, inverse, y)
    np.add.at(counts, inverse, 1.0)
    return unique_x, summed / counts


def _overlap_curves(x_ref, y_ref, x_other, y_other, n_grid: int = 40):
    """Interpolate both curves onto one grid inside their overlapping x range."""
    x_ref, y_ref = _sorted_curve(x_ref, y_ref)
    x_other, y_other = _sorted_curve(x_other, y_other)
    if len(x_ref) < 2 or len(x_other) < 2:
        return None
    lower = max(float(x_ref[0]), float(x_other[0]))
    upper = min(float(x_ref[-1]), float(x_other[-1]))
    if not np.isfinite(lower) or upper <= lower:
        return None
    grid = np.linspace(lower, upper, int(n_grid))
    reference = np.interp(grid, x_ref, y_ref)
    other = np.interp(grid, x_other, y_other)
    return x_ref, y_ref, x_other, y_other, lower, upper, reference, other


def curve_distance(x_ref, y_ref, x_other, y_other, *, n_grid: int = 40) -> float:
    """Root-mean-square difference of two curves on their overlapping x range.

    Both curves are interpolated onto one shared grid. Identical curves return 0.
    A constant vertical offset is kept: shifting a curve by 1.5 gives distance 1.5.
    """
    overlapped = _overlap_curves(x_ref, y_ref, x_other, y_other, n_grid=n_grid)
    if overlapped is None:
        return float("nan")
    _x_ref, _y_ref, _x_other, _y_other, _lower, _upper, reference, other = overlapped
    return float(np.sqrt(np.mean((reference - other) ** 2)))


def curve_fidelity(x_ref, y_ref, x_other, y_other, *, n_grid: int = 40) -> dict[str, float]:
    """Level, shape, interval-slope and direction agreement of two curves.

    Level distance is the RMSE on the shared overlap grid. Shape distance
    subtracts each curve's mean on that same grid before the RMSE. Interval
    slopes are first differences at the reference curve's nodes inside the
    overlap; the other curve is interpolated onto those nodes. Direction
    agreement is the fraction of those intervals with the same nonzero sign.
    A zero slope does not count as agreement.
    """
    empty = {
        "level_distance": float("nan"),
        "shape_distance": float("nan"),
        "slope_distance": float("nan"),
        "slope_scale": float("nan"),
        "direction_agreement": float("nan"),
    }
    overlapped = _overlap_curves(x_ref, y_ref, x_other, y_other, n_grid=n_grid)
    if overlapped is None:
        return empty
    x_left, y_left, x_right, _y_right, lower, upper, reference, other = overlapped
    level = float(np.sqrt(np.mean((reference - other) ** 2)))
    shape = float(np.sqrt(np.mean(((reference - reference.mean()) - (other - other.mean())) ** 2)))
    inside = (x_left >= lower) & (x_left <= upper)
    slopes_ref = []
    slopes_other = []
    for index in range(len(x_left) - 1):
        if not (inside[index] and inside[index + 1]):
            continue
        span = float(x_left[index + 1] - x_left[index])
        if span <= 0:
            continue
        slopes_ref.append(float(y_left[index + 1] - y_left[index]) / span)
        other_left = float(np.interp(x_left[index], x_right, _y_right))
        other_right = float(np.interp(x_left[index + 1], x_right, _y_right))
        slopes_other.append((other_right - other_left) / span)
    if not slopes_ref:
        empty["level_distance"] = level
        empty["shape_distance"] = shape
        return empty
    ref_slope = np.asarray(slopes_ref, dtype=np.float64)
    other_slope = np.asarray(slopes_other, dtype=np.float64)
    same = (np.sign(ref_slope) == np.sign(other_slope)) & (ref_slope != 0.0) & (other_slope != 0.0)
    scale = float(np.mean(np.abs(ref_slope)))
    return {
        "level_distance": level,
        "shape_distance": shape,
        "slope_distance": float(np.sqrt(np.mean((ref_slope - other_slope) ** 2))),
        "slope_scale": scale,
        "direction_agreement": float(np.mean(same)),
    }


def binned_dependence(x, y, n_bins: int = 20) -> pd.DataFrame:
    """Mean of ``y`` in quantile bins of ``x``.

    Columns are ``x`` (bin mean), ``y`` (bin mean) and ``n`` (bin count).
    This is the dependence curve used for distance and slope.
    """
    x = _as_1d(x)
    y = _as_1d(y)
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) == 0:
        return pd.DataFrame(columns=["x", "y", "n"])
    edges = np.unique(np.quantile(x, np.linspace(0.0, 1.0, int(n_bins) + 1)))
    if len(edges) < 2:
        return pd.DataFrame({"x": [float(x.mean())], "y": [float(y.mean())], "n": [int(len(x))]})
    bins = np.digitize(x, edges[1:-1], right=False)
    rows = []
    for bin_id in range(len(edges) - 1):
        chosen = bins == bin_id
        if not np.any(chosen):
            continue
        rows.append({"x": float(x[chosen].mean()), "y": float(y[chosen].mean()), "n": int(chosen.sum())})
    return pd.DataFrame(rows)


def _weighted_slope(curve: pd.DataFrame) -> float:
    if len(curve) < 2:
        return float("nan")
    x = curve["x"].to_numpy(dtype=np.float64)
    y = curve["y"].to_numpy(dtype=np.float64)
    weights = curve["n"].to_numpy(dtype=np.float64)
    if np.ptp(x) == 0 or weights.sum() <= 0:
        return float("nan")
    x_mean = np.average(x, weights=weights)
    y_mean = np.average(y, weights=weights)
    variance = np.average((x - x_mean) ** 2, weights=weights)
    if variance <= 0:
        return float("nan")
    covariance = np.average((x - x_mean) * (y - y_mean), weights=weights)
    return float(covariance / variance)


def _feature_names(names: Sequence[str] | None, n_features: int) -> list[str]:
    if names is None:
        return [f"f{i}" for i in range(n_features)]
    names = [str(name) for name in names]
    if len(names) != n_features:
        raise ValueError(f"Expected {n_features} feature names, got {len(names)}")
    return names


def _matrix(values, n_features: int | None = None) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 2:
        raise ValueError("SHAP values must have shape (n_samples, n_features)")
    if n_features is not None and array.shape[1] != n_features:
        raise ValueError("SHAP columns do not match x")
    return array


def dependence_slopes(x, shap_values, feature_names: Sequence[str] | None = None, n_bins: int = 20) -> pd.DataFrame:
    """Slope of each dependence-SHAP curve, a sensitivity of SHAP to that input.

    The slope is the weighted least-squares slope of bin-mean SHAP against
    bin-mean x. Units are SHAP per unit of the feature.
    """
    features = np.asarray(x, dtype=np.float64)
    values = _matrix(shap_values, features.shape[1])
    names = _feature_names(feature_names, features.shape[1])
    rows = []
    for index, name in enumerate(names):
        curve = binned_dependence(features[:, index], values[:, index], n_bins=n_bins)
        rows.append({"feature": name, "slope": _weighted_slope(curve)})
    return pd.DataFrame(rows)


def dependence_slope_difference(
    x_ref,
    shap_ref,
    x_other,
    shap_other,
    feature_names: Sequence[str] | None = None,
    n_bins: int = 20,
) -> pd.DataFrame:
    """Difference of dependence-SHAP slopes, ``slope_other - slope_ref``."""
    reference = dependence_slopes(x_ref, shap_ref, feature_names, n_bins=n_bins)
    other = dependence_slopes(x_other, shap_other, feature_names, n_bins=n_bins)
    merged = reference.merge(other, on="feature", suffixes=("_ref", "_other"))
    merged["difference"] = merged["slope_other"] - merged["slope_ref"]
    return merged


def dependence_fidelity(
    x_ref,
    shap_ref,
    x_other,
    shap_other,
    feature_names: Sequence[str] | None = None,
    n_bins: int = 20,
    n_grid: int = 40,
) -> pd.DataFrame:
    """Per-feature level, shape, interval-slope and direction fidelity."""
    left = np.asarray(x_ref, dtype=np.float64)
    right = np.asarray(x_other, dtype=np.float64)
    shap_left = _matrix(shap_ref, left.shape[1])
    shap_right = _matrix(shap_other, right.shape[1])
    if left.shape[1] != right.shape[1]:
        raise ValueError("The two dependence inputs have different feature counts")
    names = _feature_names(feature_names, left.shape[1])
    rows = []
    for index, name in enumerate(names):
        curve_left = binned_dependence(left[:, index], shap_left[:, index], n_bins=n_bins)
        curve_right = binned_dependence(right[:, index], shap_right[:, index], n_bins=n_bins)
        scores = curve_fidelity(
            curve_left["x"], curve_left["y"], curve_right["x"], curve_right["y"], n_grid=n_grid
        )
        scores["feature"] = name
        rows.append(scores)
    return pd.DataFrame(rows)


def dependence_curve_distance(
    x_ref,
    shap_ref,
    x_other,
    shap_other,
    feature_names: Sequence[str] | None = None,
    n_bins: int = 20,
    n_grid: int = 40,
) -> pd.DataFrame:
    """Distance between binned dependence-SHAP curves, one row per feature."""
    left = np.asarray(x_ref, dtype=np.float64)
    right = np.asarray(x_other, dtype=np.float64)
    shap_left = _matrix(shap_ref, left.shape[1])
    shap_right = _matrix(shap_other, right.shape[1])
    if left.shape[1] != right.shape[1]:
        raise ValueError("The two dependence inputs have different feature counts")
    names = _feature_names(feature_names, left.shape[1])
    rows = []
    for index, name in enumerate(names):
        curve_left = binned_dependence(left[:, index], shap_left[:, index], n_bins=n_bins)
        curve_right = binned_dependence(right[:, index], shap_right[:, index], n_bins=n_bins)
        rows.append({
            "feature": name,
            "distance": curve_distance(curve_left["x"], curve_left["y"], curve_right["x"], curve_right["y"], n_grid=n_grid),
        })
    return pd.DataFrame(rows)


def ale_curve_distance(ale_ref: dict, ale_other: dict, *, n_grid: int = 40) -> pd.DataFrame:
    """Distance between 1D ALE curves from :func:`rf_shap.explain.ale.ale_1d`."""
    rows = []
    for name in sorted(set(ale_ref) & set(ale_other)):
        left = ale_ref[name]
        right = ale_other[name]
        rows.append({
            "feature": str(name),
            "distance": curve_distance(left["x"], left["ale"], right["x"], right["ale"], n_grid=n_grid),
        })
    return pd.DataFrame(rows)


def _interaction_array(values, n_features: int) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 3 or array.shape[1] != n_features or array.shape[2] != n_features:
        raise ValueError(f"Interaction SHAP must have shape (n_samples, {n_features}, {n_features})")
    return array


def interaction_curve_distance(
    x_ref,
    interaction_ref,
    x_other,
    interaction_other,
    feature_names: Sequence[str] | None = None,
    n_bins: int = 15,
    n_grid: int = 40,
) -> pd.DataFrame:
    """Distance between interactive-SHAP curves.

    For each ordered pair ``(feature, interacts_with)`` with the two names
    different, the curve is mean interaction SHAP against ``feature``.
    The diagonal main effect is left to dependence SHAP.
    """
    left = np.asarray(x_ref, dtype=np.float64)
    right = np.asarray(x_other, dtype=np.float64)
    if left.shape[1] != right.shape[1]:
        raise ValueError("The two interaction inputs have different feature counts")
    n_features = left.shape[1]
    inter_left = _interaction_array(interaction_ref, n_features)
    inter_right = _interaction_array(interaction_other, n_features)
    names = _feature_names(feature_names, n_features)
    rows = []
    for i, feature in enumerate(names):
        for j, other_name in enumerate(names):
            if i == j:
                continue
            curve_left = binned_dependence(left[:, i], inter_left[:, i, j], n_bins=n_bins)
            curve_right = binned_dependence(right[:, i], inter_right[:, i, j], n_bins=n_bins)
            rows.append({
                "feature": feature,
                "interacts_with": other_name,
                "distance": curve_distance(
                    curve_left["x"], curve_left["y"], curve_right["x"], curve_right["y"], n_grid=n_grid
                ),
            })
    return pd.DataFrame(rows)


def interaction_shap(model, x: np.ndarray) -> np.ndarray:
    """Pairwise SHAP interaction values, shape ``(n_samples, n_features, n_features)``.

    Calls ``shap.TreeExplainer.shap_interaction_values``. Off-diagonal entries
    are the interactive SHAP used by :func:`interaction_curve_distance`.
    """
    import shap

    features = np.asarray(x, dtype=np.float64)
    explainer = shap.TreeExplainer(model)
    values = explainer.shap_interaction_values(features)
    array = np.asarray(values, dtype=np.float64)
    if array.ndim == 2:
        array = array[:, :, None]
    if array.ndim != 3:
        raise ValueError(f"Unexpected interaction shape {array.shape}")
    return array


def compare_explanations(
    *,
    ale_ref: dict | None = None,
    ale_other: dict | None = None,
    x_ref=None,
    shap_ref=None,
    x_other=None,
    shap_other=None,
    interaction_ref=None,
    interaction_other=None,
    feature_names: Sequence[str] | None = None,
    n_bins: int = 20,
) -> dict[str, pd.DataFrame]:
    """Distances and sensitivity differences between a reference and another explanation.

    Provide the pairs you have. Returned keys are among ``ale_distance``,
    ``dependence_distance``, ``interaction_distance`` and ``slope_difference``.
    """
    compared: dict[str, pd.DataFrame] = {}
    if ale_ref is not None or ale_other is not None:
        if ale_ref is None or ale_other is None:
            raise ValueError("ALE comparison needs both ale_ref and ale_other")
        compared["ale_distance"] = ale_curve_distance(ale_ref, ale_other)
    if shap_ref is not None or shap_other is not None:
        if any(value is None for value in (x_ref, shap_ref, x_other, shap_other)):
            raise ValueError("Dependence comparison needs x and SHAP for both explanations")
        compared["dependence_distance"] = dependence_curve_distance(
            x_ref, shap_ref, x_other, shap_other, feature_names, n_bins=n_bins
        )
        compared["slope_difference"] = dependence_slope_difference(
            x_ref, shap_ref, x_other, shap_other, feature_names, n_bins=n_bins
        )
    if interaction_ref is not None or interaction_other is not None:
        if any(value is None for value in (x_ref, interaction_ref, x_other, interaction_other)):
            raise ValueError("Interaction comparison needs x and interaction SHAP for both explanations")
        compared["interaction_distance"] = interaction_curve_distance(
            x_ref, interaction_ref, x_other, interaction_other, feature_names, n_bins=min(n_bins, 15)
        )
    if not compared:
        raise ValueError("Provide ALE curves, dependence SHAP, or interaction SHAP to compare")
    return compared
