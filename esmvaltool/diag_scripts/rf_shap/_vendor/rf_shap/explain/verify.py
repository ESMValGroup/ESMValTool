"""Controlled checks of SHAP comparison metrics (GMD Fig. 2).

The true map is nonlinear. Three alternatives share the same inputs
(corr(x1, x2) = 0.2):

* ``m1`` linear least-squares approximation
* ``m2`` sign flip of the nonlinear term
* ``m3`` same family as the truth but weaker curvature

Each map is fit with ExtraTrees + TreeSHAP. The main figure compares
``x`` versus ``y``, SHAP dependence, and metric distance to the truth.
ALE distances are stored for the supplement.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from rf_shap.cube.tabular import SampleTable
from rf_shap.explain.ale import ale_1d
from rf_shap.explain.evaluate import binned_dependence, compare_explanations
from rf_shap.explain.shap_explain import compute_shap
from rf_shap.models.registry import build_model

FEATURE_NAMES = ("x1", "x2")
COEFF_A = 1.5
COEFF_B = 0.8
WEAK_CURVATURE = 0.5 * COEFF_A
INPUT_CORR = 0.2

MODELS: dict[str, dict[str, Any]] = {
    "m1": {"label": "m1 linear", "kind": "linear_approx"},
    "m2": {"label": "m2 sign flip", "kind": "sign_flip"},
    "m3": {"label": "m3 weak curve", "kind": "weak_curve"},
}

# Kept so older imports of CASES still resolve to the three models.
CASES = MODELS
EXPECTED_SLOPE_SIGN_X1: dict[str, int] = {}


def sample_verification_inputs(
    n: int,
    *,
    corr: float = INPUT_CORR,
    loc: float = 0.0,
    scale: float = 1.0,
    seed: int = 0,
) -> np.ndarray:
    """Draw bivariate Gaussian (x1, x2) with the requested correlation."""
    rng = np.random.default_rng(seed)
    cov = np.array([[scale**2, corr * scale], [corr * scale, 1.0]], dtype=np.float64)
    mean = np.array([loc, 0.0], dtype=np.float64)
    return rng.multivariate_normal(mean, cov, size=int(n))


def verification_target(x: np.ndarray, kind: str = "true") -> np.ndarray:
    """Known maps. ``true`` is quadratic in x1; ``linear_approx`` is not used here."""
    x1 = np.asarray(x, dtype=np.float64)[:, 0]
    x2 = np.asarray(x, dtype=np.float64)[:, 1]
    if kind in {"true", "nonlinear"}:
        return COEFF_A * (x1**2) + COEFF_B * x2
    if kind == "sign_flip":
        return -COEFF_A * (x1**2) + COEFF_B * x2
    if kind == "weak_curve":
        return WEAK_CURVATURE * (x1**2) + COEFF_B * x2
    if kind == "linear":
        return COEFF_A * x1 + COEFF_B * x2
    raise ValueError(f"Unknown target kind {kind!r}")


def linear_approximation(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """OLS fit of y on (1, x1, x2), the linear approximation of the true map."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64).reshape(-1)
    design = np.column_stack([np.ones(len(y)), x[:, 0], x[:, 1]])
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    return design @ coef


def model_target(x: np.ndarray, kind: str, y_true: np.ndarray | None = None) -> np.ndarray:
    if kind == "true":
        return verification_target(x, "true")
    if kind == "linear_approx":
        if y_true is None:
            y_true = verification_target(x, "true")
        return linear_approximation(x, y_true)
    return verification_target(x, kind)


def _table(x: np.ndarray, y: np.ndarray) -> SampleTable:
    n = len(y)
    index = pd.MultiIndex.from_arrays(
        [np.arange(n), np.zeros(n), np.zeros(n)],
        names=("time", "lat", "lon"),
    )
    return SampleTable(x, y, list(FEATURE_NAMES), index)


def fit_and_explain(
    x: np.ndarray,
    y: np.ndarray,
    *,
    seed: int,
    noise: float = 0.15,
    n_estimators: int = 80,
    max_depth: int = 8,
    max_samples: int | None = 2000,
    n_jobs: int = -1,
) -> dict[str, Any]:
    """Fit ExtraTrees on a noisy copy of y = f(x) and return SHAP plus ALE."""
    rng = np.random.default_rng(seed + 17)
    y_obs = np.asarray(y, dtype=np.float64) + noise * rng.normal(size=len(y))
    model = build_model(
        "extra_trees",
        params={"n_estimators": n_estimators, "max_depth": max_depth, "n_jobs": n_jobs, "min_samples_leaf": 2},
        seed=seed,
    )
    model.fit(x, y_obs)
    table = _table(x, y_obs)
    shap = compute_shap(model, table, method="auto", max_samples=max_samples, seed=seed, n_jobs=n_jobs)
    positions = np.asarray(shap.sample_positions)
    x_used = np.asarray(x, dtype=np.float64)[positions]
    pred = np.asarray(model.predict(x_used), dtype=np.float64)
    expected = np.asarray(shap.expected_value, dtype=np.float64).reshape(-1)
    if expected.size == 1:
        expected = np.full(len(pred), float(expected[0]))
    recon = shap.values.sum(axis=1) + expected
    add_rmse = float(np.sqrt(np.mean((pred - recon) ** 2)))
    ale = ale_1d(model, x_used, FEATURE_NAMES, n_grid=20)
    return {
        "x": x_used,
        "y": y_obs[positions],
        "y_map": np.asarray(y, dtype=np.float64)[positions],
        "shap": shap.values,
        "ale": ale,
        "additivity_rmse": add_rmse,
        "pred_rmse": float(np.sqrt(np.mean((pred - y_obs[positions]) ** 2))),
    }


def _compare_pair(ref: dict[str, Any], other: dict[str, Any], n_bins: int = 20) -> dict[str, pd.DataFrame]:
    return compare_explanations(
        ale_ref=ref["ale"],
        ale_other=other["ale"],
        x_ref=ref["x"],
        shap_ref=ref["shap"],
        x_other=other["x"],
        shap_other=other["shap"],
        feature_names=list(FEATURE_NAMES),
        n_bins=n_bins,
    )


def _curve(x: np.ndarray, values: np.ndarray, feature: str = "x1", n_bins: int = 20) -> pd.DataFrame:
    j = FEATURE_NAMES.index(feature)
    curve = binned_dependence(x[:, j], values, n_bins=n_bins)
    return curve.copy()


def run_attribution_verification(
    *,
    n: int = 2500,
    repeats: int = 12,
    n_bins: int = 20,
    noise: float = 0.15,
    n_estimators: int = 80,
    max_depth: int = 8,
    max_samples: int | None = 2000,
    n_jobs: int = -1,
    seed0: int = 0,
) -> dict[str, Any]:
    """Fit the true nonlinear map and m1–m3 on shared inputs; compare SHAP to truth."""
    metric_rows: list[dict[str, Any]] = []
    add_rows: list[dict[str, Any]] = []
    curves: dict[str, pd.DataFrame] = {}
    display_seed = seed0

    for repeat in range(int(repeats)):
        seed = seed0 + repeat
        x = sample_verification_inputs(n, seed=seed, corr=INPUT_CORR)
        y_true = model_target(x, "true")
        truth = fit_and_explain(
            x,
            y_true,
            seed=seed,
            noise=noise,
            n_estimators=n_estimators,
            max_depth=max_depth,
            max_samples=max_samples,
            n_jobs=n_jobs,
        )
        add_rows.append({"seed": seed, "case": "true", "additivity_rmse": truth["additivity_rmse"]})
        if seed == display_seed:
            curves["true:y"] = _curve(x, y_true, n_bins=n_bins)
            curves["true:x1"] = _curve(truth["x"], truth["shap"][:, 0], n_bins=n_bins)
            curves["true:ale:x1"] = truth["ale"]["x1"].assign(feature="x1")
            curves["true:scatter"] = pd.DataFrame({"x": x[:, 0], "y": y_true})

        for case, spec in MODELS.items():
            y_model = model_target(x, str(spec["kind"]), y_true=y_true)
            other = fit_and_explain(
                x,
                y_model,
                seed=seed + 3,
                noise=noise,
                n_estimators=n_estimators,
                max_depth=max_depth,
                max_samples=max_samples,
                n_jobs=n_jobs,
            )
            add_rows.append({"seed": seed, "case": case, "additivity_rmse": other["additivity_rmse"]})
            compared = _compare_pair(truth, other, n_bins=n_bins)
            dep = compared["dependence_distance"].set_index("feature")
            slope = compared["slope_difference"].set_index("feature")
            ale = compared["ale_distance"].set_index("feature")
            for feature in FEATURE_NAMES:
                expected = EXPECTED_SLOPE_SIGN_X1.get(case) if feature == "x1" else None
                difference = float(slope.loc[feature, "difference"])
                metric_rows.append({
                    "seed": seed,
                    "case": case,
                    "label": spec["label"],
                    "feature": feature,
                    "dependence_distance": float(dep.loc[feature, "distance"]),
                    "ale_distance": float(ale.loc[feature, "distance"]),
                    "slope_ref": float(slope.loc[feature, "slope_ref"]),
                    "slope_other": float(slope.loc[feature, "slope_other"]),
                    "slope_difference": difference,
                    "expected_sign": expected,
                    "sign_hit": (
                        float(np.sign(difference) == expected) if expected is not None else np.nan
                    ),
                })
            if seed == display_seed:
                curves[f"{case}:y"] = _curve(x, y_model, n_bins=n_bins)
                curves[f"{case}:x1"] = _curve(other["x"], other["shap"][:, 0], n_bins=n_bins)
                curves[f"{case}:ale:x1"] = other["ale"]["x1"].assign(feature="x1")

    metrics = pd.DataFrame(metric_rows)
    additivity = pd.DataFrame(add_rows)
    summary = (
        metrics.groupby(["case", "label", "feature"], sort=False)
        .agg(
            dependence_distance=("dependence_distance", "median"),
            dependence_q10=("dependence_distance", lambda s: float(s.quantile(0.1))),
            dependence_q90=("dependence_distance", lambda s: float(s.quantile(0.9))),
            ale_distance=("ale_distance", "median"),
            ale_q10=("ale_distance", lambda s: float(s.quantile(0.1))),
            ale_q90=("ale_distance", lambda s: float(s.quantile(0.9))),
            slope_difference=("slope_difference", "median"),
            slope_q10=("slope_difference", lambda s: float(s.quantile(0.1))),
            slope_q90=("slope_difference", lambda s: float(s.quantile(0.9))),
            sign_hit=("sign_hit", "mean"),
        )
        .reset_index()
    )
    return {"metrics": metrics, "summary": summary, "curves": curves, "additivity": additivity}
