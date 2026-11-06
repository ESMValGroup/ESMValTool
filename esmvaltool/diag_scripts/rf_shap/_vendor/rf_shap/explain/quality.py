"""Quality metrics for SHAP, ALE and PDP explanations.

Properties scored without a ground-truth explanation:

* **Faithfulness** — does the attribution track the model's response when
  features are changed?
* **Locality / complexity** — nearby samples, sparse attributions, effective
  number of drivers
* **Robustness / stability** — noise, bootstrap, or retrain/seed changes
* **Contrastivity** — do high- and low-target samples get different SHAP
* **Agreement of global views** — ALE vs PDP vs dependence-SHAP distances
* **Interaction strength** — off-diagonal interaction mass vs main effects

References: Yeh et al. (2019); Bhatt et al. (2020); Ancona et al. (2018)
Sensitivity-n; OpenXAI PGI/PGU; Agarwal et al. (2022) RIS.
"""

from __future__ import annotations

from typing import Callable, Sequence

import numpy as np
import pandas as pd

from rf_shap.explain.ale import ale_1d
from rf_shap.explain.evaluate import binned_dependence, curve_distance
from rf_shap.explain.ice_pdp import ice_pdp


def _matrix(values) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 2:
        raise ValueError(
            "Expected a 2-D array of shape (n_samples, n_features)"
        )
    return array


def _feature_names(names: Sequence[str] | None, n_features: int) -> list[str]:
    if names is None:
        return [f"f{i}" for i in range(n_features)]
    names = [str(name) for name in names]
    if len(names) != n_features:
        raise ValueError(
            f"Expected {n_features} feature names, got {len(names)}"
        )
    return names


def _predict(model, x: np.ndarray) -> np.ndarray:
    return np.asarray(
        model.predict(np.asarray(x, dtype=np.float64)), dtype=np.float64
    ).reshape(-1)


def _standardize(x: np.ndarray) -> np.ndarray:
    scale = np.nanstd(x, axis=0)
    scale = np.where(scale > 0, scale, 1.0)
    return (x - np.nanmean(x, axis=0)) / scale


def shap_additivity(model, x, shap_values, expected_value=None) -> float:
    """RMSE of ``sum(SHAP) + E[f]`` against the model prediction.

    Additive TreeSHAP explanations should be near zero. A large value means
    the attributions do not reconstruct the fitted map.
    """
    features = _matrix(x)
    values = _matrix(shap_values)
    if features.shape != values.shape:
        raise ValueError("x and SHAP must share the same shape")
    pred = _predict(model, features)
    if expected_value is None:
        baseline = float(np.mean(pred))
    else:
        expected = np.asarray(expected_value, dtype=np.float64).reshape(-1)
        baseline = float(expected[0]) if expected.size == 1 else expected
        if np.ndim(baseline) == 0:
            baseline = float(baseline)
        else:
            baseline = np.asarray(baseline, dtype=np.float64)
            if baseline.shape[0] != len(pred):
                raise ValueError(
                    "expected_value length must match the number of samples"
                )
    recon = values.sum(axis=1) + baseline
    return float(np.sqrt(np.mean((pred - recon) ** 2)))


def shap_faithfulness_correlation(
    model,
    x,
    shap_values,
    *,
    n_perturb: int = 50,
    seed: int = 0,
) -> float:
    """Spearman correlation between attributions and single-feature impacts.

    For random (sample, feature) pairs the feature is replaced by its column
    mean. The prediction change is correlated with the absolute attribution.
    Higher is more faithful (Bhatt et al., 2020 / Alvarez-Melis & Jaakkola).
    """
    features = _matrix(x)
    values = _matrix(shap_values)
    if features.shape != values.shape:
        raise ValueError("x and SHAP must share the same shape")
    n, p = features.shape
    if n == 0 or p == 0:
        return float("nan")
    rng = np.random.default_rng(seed)
    means = np.nanmean(features, axis=0)
    pred = _predict(model, features)
    n_perturb = int(min(max(n_perturb, 1), n * p))
    impacts = np.empty(n_perturb, dtype=np.float64)
    attributions = np.empty(n_perturb, dtype=np.float64)
    for index in range(n_perturb):
        i = int(rng.integers(0, n))
        j = int(rng.integers(0, p))
        probe = features[i].copy()
        probe[j] = means[j]
        impacts[index] = abs(
            float(pred[i] - _predict(model, probe[None, :])[0])
        )
        attributions[index] = abs(float(values[i, j]))
    if np.allclose(impacts, impacts[0]) or np.allclose(
        attributions, attributions[0]
    ):
        return float("nan")
    return float(
        pd.Series(attributions).corr(pd.Series(impacts), method="spearman")
    )


def shap_prediction_gaps(
    model,
    x,
    shap_values,
    *,
    top_k: int | None = None,
    seed: int = 0,
) -> dict[str, float]:
    """PGI and PGU: prediction change after masking important / unimportant features.

    Important features (largest ``|SHAP|``) are replaced by column means for
    PGI; the complementary set is replaced for PGU. Faithful attributions give
    a large PGI and a small PGU.
    """
    features = _matrix(x)
    values = _matrix(shap_values)
    if features.shape != values.shape:
        raise ValueError("x and SHAP must share the same shape")
    n, p = features.shape
    if p < 2:
        return {
            "pgi": float("nan"),
            "pgu": float("nan"),
            "top_k": float("nan"),
        }
    k = int(top_k) if top_k is not None else max(1, p // 2)
    k = min(max(k, 1), p - 1)
    means = np.nanmean(features, axis=0)
    pred = _predict(model, features)
    order = np.argsort(-np.abs(values), axis=1)
    important = np.zeros_like(features, dtype=bool)
    rows = np.arange(n)[:, None]
    important[rows, order[:, :k]] = True
    x_important = features.copy()
    x_unimportant = features.copy()
    x_important[important] = np.take(means, np.where(important)[1])
    x_unimportant[~important] = np.take(means, np.where(~important)[1])
    pgi = float(np.mean(np.abs(pred - _predict(model, x_important))))
    pgu = float(np.mean(np.abs(pred - _predict(model, x_unimportant))))
    return {"pgi": pgi, "pgu": pgu, "top_k": float(k)}


def shap_infidelity(
    model,
    x,
    shap_values,
    *,
    n_perturb: int = 40,
    noise: float = 0.1,
    seed: int = 0,
) -> float:
    """Mean squared gap between ``phi · δ`` and ``f(x) - f(x+δ)`` (Yeh et al.).

    Lower is better. ``noise`` is a fraction of each feature's standard deviation.
    """
    features = _matrix(x)
    values = _matrix(shap_values)
    if features.shape != values.shape:
        raise ValueError("x and SHAP must share the same shape")
    n, p = features.shape
    rng = np.random.default_rng(seed)
    scale = np.nanstd(features, axis=0)
    scale = np.where(scale > 0, scale, 1.0)
    pred = _predict(model, features)
    scores = []
    for _ in range(int(n_perturb)):
        delta = rng.normal(0.0, 1.0, size=(n, p)) * (noise * scale)
        moved = features + delta
        change = pred - _predict(model, moved)
        approx = np.sum(values * delta, axis=1)
        scores.append(float(np.mean((approx - change) ** 2)))
    return float(np.mean(scores)) if scores else float("nan")


def explanation_locality(
    x,
    shap_values,
    *,
    k: int = 5,
    seed: int = 0,
) -> dict[str, float]:
    """Locality: neighbor SHAP cosine minus random-pair cosine.

    Features are standardised before nearest-neighbour search. A positive
    ``locality_gap`` means nearby samples share more similar attributions than
    random pairs.
    """
    features = _matrix(x)
    values = _matrix(shap_values)
    if features.shape[0] != values.shape[0]:
        raise ValueError("x and SHAP must have the same number of samples")
    n = features.shape[0]
    k = int(min(max(k, 1), max(n - 1, 1)))
    if n < 3:
        return {
            "neighbor_cosine": float("nan"),
            "random_cosine": float("nan"),
            "locality_gap": float("nan"),
        }
    z = _standardize(features)
    distance = np.sqrt(((z[:, None, :] - z[None, :, :]) ** 2).sum(axis=2))
    np.fill_diagonal(distance, np.inf)
    neighbors = np.argpartition(distance, k, axis=1)[:, :k]
    value_norm = values / np.clip(
        np.linalg.norm(values, axis=1, keepdims=True), 1e-12, None
    )
    similarity = value_norm @ value_norm.T
    neighbor_score = float(
        np.take_along_axis(similarity, neighbors, axis=1).mean()
    )
    rng = np.random.default_rng(seed)
    left = rng.integers(0, n, size=n * k)
    right = rng.integers(0, n, size=n * k)
    random_score = float(similarity[left, right].mean())
    return {
        "neighbor_cosine": neighbor_score,
        "random_cosine": random_score,
        "locality_gap": neighbor_score - random_score,
    }


def attribution_sparseness(shap_values) -> float:
    """Gini coefficient of mean ``|SHAP|`` across features (higher = sparser)."""
    values = np.mean(np.abs(_matrix(shap_values)), axis=0)
    values = np.sort(np.asarray(values, dtype=np.float64))
    total = values.sum()
    if total <= 0 or len(values) < 2:
        return float("nan")
    index = np.arange(1, len(values) + 1, dtype=np.float64)
    return float(
        (2.0 * np.sum(index * values) / (len(values) * total))
        - (len(values) + 1.0) / len(values)
    )


def shap_sensitivity_n(
    model,
    x,
    shap_values,
    *,
    n: int = 2,
    n_draws: int = 40,
    seed: int = 0,
) -> float:
    """Sensitivity-n: correlation of subset attribution sums with prediction gaps.

    Random subsets of size ``n`` are masked to column means. Spearman
    correlation is taken between ``sum(phi_S)`` and ``f(x) - f(x with S
    replaced)``. Higher is more complete (Ancona et al., 2018).
    """
    features = _matrix(x)
    values = _matrix(shap_values)
    if features.shape != values.shape:
        raise ValueError("x and SHAP must share the same shape")
    n_samples, n_features = features.shape
    subset_size = int(min(max(n, 1), n_features))
    if n_samples == 0 or n_features == 0:
        return float("nan")
    rng = np.random.default_rng(seed)
    means = np.nanmean(features, axis=0)
    pred = _predict(model, features)
    attribution_sums = np.empty(n_draws, dtype=np.float64)
    prediction_gaps = np.empty(n_draws, dtype=np.float64)
    for draw in range(int(n_draws)):
        i = int(rng.integers(0, n_samples))
        chosen = rng.choice(n_features, size=subset_size, replace=False)
        probe = features[i].copy()
        probe[chosen] = means[chosen]
        attribution_sums[draw] = float(values[i, chosen].sum())
        prediction_gaps[draw] = float(
            pred[i] - _predict(model, probe[None, :])[0]
        )
    if np.allclose(attribution_sums, attribution_sums[0]) or np.allclose(
        prediction_gaps, prediction_gaps[0]
    ):
        return float("nan")
    return float(
        pd.Series(attribution_sums).corr(
            pd.Series(prediction_gaps), method="spearman"
        )
    )


def effective_complexity(
    shap_values,
    *,
    mass: float = 0.9,
) -> dict[str, float]:
    """How many features carry most of the mean absolute attribution.

    ``effective_n`` is the smallest k such that the top-k features cover
    ``mass`` of total mean ``|SHAP|``. ``top_k_mass`` is the fraction carried
    by the first ``ceil(p/2)`` features. Higher sparseness / lower effective_n
    means a simpler local story.
    """
    values = np.mean(np.abs(_matrix(shap_values)), axis=0)
    total = float(values.sum())
    p = len(values)
    if total <= 0 or p == 0:
        return {
            "effective_n": float("nan"),
            "top_k_mass": float("nan"),
            "mass": float(mass),
        }
    order = np.sort(values)[::-1]
    cumulative = np.cumsum(order) / total
    threshold = float(np.clip(mass, 0.0, 1.0))
    effective_n = int(np.searchsorted(cumulative, threshold) + 1)
    half = max(1, (p + 1) // 2)
    return {
        "effective_n": float(effective_n),
        "top_k_mass": float(order[:half].sum() / total),
        "mass": threshold,
    }


def explanation_contrastivity(
    shap_values,
    y,
    *,
    low_q: float = 0.25,
    high_q: float = 0.75,
) -> dict[str, float]:
    """Distance between mean SHAP of low-y and high-y samples.

    Larger ``contrastivity`` means the explanation separates the tails of the
    target. Also returns the cosine between the two mean attribution vectors.
    """
    values = _matrix(shap_values)
    target = np.asarray(y, dtype=np.float64).reshape(-1)
    if len(target) != values.shape[0]:
        raise ValueError("y and SHAP must have the same number of samples")
    low = float(np.nanquantile(target, low_q))
    high = float(np.nanquantile(target, high_q))
    low_mask = target <= low
    high_mask = target >= high
    if low_mask.sum() < 2 or high_mask.sum() < 2:
        return {
            "contrastivity": float("nan"),
            "contrast_cosine": float("nan"),
            "n_low": float(low_mask.sum()),
            "n_high": float(high_mask.sum()),
        }
    mean_low = values[low_mask].mean(axis=0)
    mean_high = values[high_mask].mean(axis=0)
    denom = float(np.linalg.norm(mean_low) * np.linalg.norm(mean_high))
    cosine = (
        float(np.dot(mean_low, mean_high) / denom)
        if denom > 0
        else float("nan")
    )
    return {
        "contrastivity": float(np.linalg.norm(mean_high - mean_low)),
        "contrast_cosine": cosine,
        "n_low": float(low_mask.sum()),
        "n_high": float(high_mask.sum()),
    }


def seed_stability(
    explanation_list: Sequence[np.ndarray] | Sequence[dict],
    *,
    kind: str = "shap",
    n_grid: int = 40,
) -> dict[str, float] | pd.DataFrame:
    """Agreement across explanations from different seeds or retrains.

    For ``kind="shap"`` pass a sequence of ``(n_samples, n_features)`` arrays
    aligned on the same rows; returns mean pairwise cosine of flattened
    matrices and mean pairwise RMSE. For ``kind="ale"`` or ``kind="pdp"``
    pass dicts like ``ale_1d`` / ``ice_pdp`` and get a per-feature table of
    mean pairwise curve RMSE.
    """
    if len(explanation_list) < 2:
        raise ValueError("seed_stability needs at least two explanations")
    if kind == "shap":
        matrices = [_matrix(item) for item in explanation_list]
        shape0 = matrices[0].shape
        if any(item.shape != shape0 for item in matrices):
            raise ValueError("All SHAP arrays must share the same shape")
        cosines = []
        rmses = []
        for i in range(len(matrices)):
            for j in range(i + 1, len(matrices)):
                a = matrices[i].reshape(-1)
                b = matrices[j].reshape(-1)
                denom = float(np.linalg.norm(a) * np.linalg.norm(b))
                cosines.append(
                    float(np.dot(a, b) / denom) if denom > 0 else float("nan")
                )
                rmses.append(float(np.sqrt(np.mean((a - b) ** 2))))
        return {
            "mean_pairwise_cosine": float(np.nanmean(cosines)),
            "mean_pairwise_rmse": float(np.nanmean(rmses)),
            "n_runs": float(len(matrices)),
        }
    if kind not in {"ale", "pdp"}:
        raise ValueError("kind must be 'shap', 'ale', or 'pdp'")
    names = sorted(
        set.intersection(*[set(map(str, item)) for item in explanation_list])
    )
    rows = []
    for name in names:
        distances = []
        for i in range(len(explanation_list)):
            for j in range(i + 1, len(explanation_list)):
                left = explanation_list[i][name]
                right = explanation_list[j][name]
                if kind == "ale":
                    distances.append(
                        curve_distance(
                            left["x"],
                            left["ale"],
                            right["x"],
                            right["ale"],
                            n_grid=n_grid,
                        )
                    )
                else:
                    distances.append(
                        curve_distance(
                            left["grid"],
                            left["pdp"],
                            right["grid"],
                            right["pdp"],
                            n_grid=n_grid,
                        )
                    )
        rows.append(
            {
                "feature": name,
                "mean_pairwise_distance": float(np.nanmean(distances))
                if distances
                else float("nan"),
                "n_pairs": int(len(distances)),
            }
        )
    return pd.DataFrame(rows)


def global_method_agreement(
    *,
    ale_dict: dict | None = None,
    ice_pdp_dict: dict | None = None,
    x=None,
    shap_values=None,
    feature_names: Sequence[str] | None = None,
    n_bins: int = 20,
    n_grid: int = 40,
) -> pd.DataFrame:
    """Pairwise curve distances among ALE, PDP and dependence-SHAP.

    Provide at least two of the three views. Dependence curves are quantile
    bins of SHAP against each feature. Lower distance means the global
    methods tell a more similar story.
    """
    curves: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]] = {}
    if ale_dict is not None:
        curves["ale"] = {
            str(name): (
                np.asarray(frame["x"], dtype=np.float64),
                np.asarray(frame["ale"], dtype=np.float64),
            )
            for name, frame in ale_dict.items()
        }
    if ice_pdp_dict is not None:
        curves["pdp"] = {
            str(name): (
                np.asarray(payload["grid"], dtype=np.float64),
                np.asarray(payload["pdp"], dtype=np.float64),
            )
            for name, payload in ice_pdp_dict.items()
        }
    if shap_values is not None:
        if x is None:
            raise ValueError(
                "Dependence agreement needs both x and shap_values"
            )
        features = _matrix(x)
        values = _matrix(shap_values)
        names = _feature_names(feature_names, features.shape[1])
        dep = {}
        for index, name in enumerate(names):
            curve = binned_dependence(
                features[:, index], values[:, index], n_bins=n_bins
            )
            dep[name] = (
                curve["x"].to_numpy(dtype=np.float64),
                curve["y"].to_numpy(dtype=np.float64),
            )
        curves["dependence"] = dep
    methods = list(curves)
    if len(methods) < 2:
        raise ValueError(
            "Provide at least two of ale_dict, ice_pdp_dict, and shap_values"
        )
    rows = []
    for i, left_name in enumerate(methods):
        for right_name in methods[i + 1 :]:
            shared = sorted(set(curves[left_name]) & set(curves[right_name]))
            for feature in shared:
                x_left, y_left = curves[left_name][feature]
                x_right, y_right = curves[right_name][feature]
                # Center each curve so level offsets from different methods
                # do not dominate the agreement score.
                y_left = y_left - np.nanmean(y_left)
                y_right = y_right - np.nanmean(y_right)
                rows.append(
                    {
                        "feature": feature,
                        "pair": f"{left_name}_vs_{right_name}",
                        "distance": curve_distance(
                            x_left, y_left, x_right, y_right, n_grid=n_grid
                        ),
                    }
                )
    return pd.DataFrame(rows)


def interaction_strength(
    interaction_values,
    feature_names: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Off-diagonal interaction mass relative to main-effect mass.

    ``interaction_values`` has shape ``(n_samples, p, p)``. For each feature
    the strength is mean ``|phi_{i,j}|`` over ``j != i``, divided by
    mean ``|phi_{i,i}|``. Higher means interactions matter more than the
    diagonal main effect for that feature.
    """
    array = np.asarray(interaction_values, dtype=np.float64)
    if array.ndim != 3 or array.shape[1] != array.shape[2]:
        raise ValueError(
            "interaction_values must have shape (n_samples, n_features, n_features)"
        )
    n_features = array.shape[1]
    names = _feature_names(feature_names, n_features)
    abs_mean = np.mean(np.abs(array), axis=0)
    rows = []
    off_total = 0.0
    diag_total = 0.0
    for i, name in enumerate(names):
        main = float(abs_mean[i, i])
        off = float((np.sum(abs_mean[i, :]) - main) / max(n_features - 1, 1))
        off_total += float(np.sum(abs_mean[i, :]) - main)
        diag_total += main
        rows.append(
            {
                "feature": name,
                "main_effect": main,
                "mean_interaction": off,
                "interaction_ratio": float(off / main)
                if main > 0
                else float("nan"),
            }
        )
    overall = (
        float(off_total / (off_total + diag_total))
        if (off_total + diag_total) > 0
        else float("nan")
    )
    frame = pd.DataFrame(rows)
    frame["overall_interaction_fraction"] = overall
    return frame


def shap_sensitivity(
    x,
    shap_values,
    perturbed_shap_values,
    *,
    x_perturbed=None,
) -> dict[str, float]:
    """Max and average L2 change of SHAP under a set of perturbed explanations.

    ``perturbed_shap_values`` has shape ``(n_perturb, n_samples, n_features)``.
    Relative input stability (RIS) is returned when ``x_perturbed`` is given.
    Lower sensitivity / RIS is more robust.
    """
    base = _matrix(shap_values)
    perturbed = np.asarray(perturbed_shap_values, dtype=np.float64)
    if perturbed.ndim != 3 or perturbed.shape[1:] != base.shape:
        raise ValueError(
            "perturbed_shap_values must have shape (n_perturb, n_samples, n_features)"
        )
    deltas = np.linalg.norm(perturbed - base[None, :, :], axis=2)
    sample_max = deltas.max(axis=0)
    sample_mean = deltas.mean(axis=0)
    out = {
        "max_sensitivity": float(np.mean(sample_max)),
        "avg_sensitivity": float(np.mean(sample_mean)),
        "relative_input_stability": float("nan"),
    }
    if x_perturbed is None:
        return out
    features = _matrix(x)
    moved = np.asarray(x_perturbed, dtype=np.float64)
    if moved.ndim != 3 or moved.shape[1:] != features.shape:
        raise ValueError(
            "x_perturbed must have shape (n_perturb, n_samples, n_features)"
        )
    input_delta = np.linalg.norm(moved - features[None, :, :], axis=2)
    ratio = deltas / np.clip(input_delta, 1e-8, None)
    out["relative_input_stability"] = float(np.mean(ratio.max(axis=0)))
    return out


def shap_robustness(
    explain_fn: Callable[[np.ndarray], np.ndarray],
    x,
    shap_values,
    *,
    n_perturb: int = 10,
    noise: float = 0.05,
    seed: int = 0,
) -> dict[str, float]:
    """Recompute SHAP on noisy copies of ``x`` and report sensitivity / RIS.

    ``explain_fn`` maps a feature matrix to SHAP values of the same shape.
    """
    features = _matrix(x)
    base = _matrix(shap_values)
    if features.shape != base.shape:
        raise ValueError("x and SHAP must share the same shape")
    rng = np.random.default_rng(seed)
    scale = np.nanstd(features, axis=0)
    scale = np.where(scale > 0, scale, 1.0)
    moved = []
    explained = []
    for _ in range(int(n_perturb)):
        noise_matrix = rng.normal(0.0, 1.0, size=features.shape) * (
            noise * scale
        )
        probe = features + noise_matrix
        values = _matrix(explain_fn(probe))
        if values.shape != base.shape:
            raise ValueError(
                "explain_fn must return SHAP with the same shape as shap_values"
            )
        moved.append(probe)
        explained.append(values)
    return shap_sensitivity(
        features, base, np.stack(explained), x_perturbed=np.stack(moved)
    )


def global_curve_robustness(
    rebuild_fn: Callable[[np.ndarray], dict],
    x,
    reference: dict,
    *,
    n_bootstrap: int = 20,
    sample_fraction: float = 0.8,
    value_key: str = "ale",
    seed: int = 0,
) -> pd.DataFrame:
    """Bootstrap stability of ALE or PDP curves.

    ``rebuild_fn(x_subset)`` must return the same dict layout as ``reference``
    (ALE: DataFrame with columns ``x`` and ``ale``; PDP: dict with ``grid`` and
    ``pdp``). Distance is RMSE on the overlapping x range. Lower is stabler.
    """
    features = _matrix(x)
    n = features.shape[0]
    take = max(2, int(round(sample_fraction * n)))
    rng = np.random.default_rng(seed)
    rows = []
    for name, ref in reference.items():
        distances = []
        for _ in range(int(n_bootstrap)):
            index = rng.choice(n, size=take, replace=True)
            rebuilt = rebuild_fn(features[index])
            if name not in rebuilt:
                continue
            other = rebuilt[name]
            if value_key == "ale":
                distances.append(
                    curve_distance(
                        ref["x"], ref["ale"], other["x"], other["ale"]
                    )
                )
            else:
                distances.append(
                    curve_distance(
                        ref["grid"], ref["pdp"], other["grid"], other["pdp"]
                    )
                )
        rows.append(
            {
                "feature": str(name),
                "mean_distance": float(np.nanmean(distances))
                if distances
                else float("nan"),
                "max_distance": float(np.nanmax(distances))
                if distances
                else float("nan"),
                "n_bootstrap": int(len(distances)),
            }
        )
    return pd.DataFrame(rows)


def ale_robustness(
    model, x, ale_dict, names: Sequence[str], **kwargs
) -> pd.DataFrame:
    """Bootstrap ALE curve distances (resampling robustness)."""
    features = _matrix(x)
    feature_names = _feature_names(names, features.shape[1])

    def rebuild(subset):
        return ale_1d(model, subset, feature_names)

    return global_curve_robustness(
        rebuild, features, ale_dict, value_key="ale", **kwargs
    )


def pdp_robustness(
    model, x, ice_pdp_dict, names: Sequence[str], **kwargs
) -> pd.DataFrame:
    """Bootstrap PDP curve distances (resampling robustness)."""
    features = _matrix(x)
    feature_names = _feature_names(names, features.shape[1])

    def rebuild(subset):
        return ice_pdp(model, subset, feature_names)

    return global_curve_robustness(
        rebuild, features, ice_pdp_dict, value_key="pdp", **kwargs
    )


def ice_heterogeneity(ice_pdp_dict) -> pd.DataFrame:
    """Locality-related score for PDP: mean absolute ICE deviation from PDP.

    A large value means the global PDP hides strong individual heterogeneity.
    """
    rows = []
    for name, payload in ice_pdp_dict.items():
        ice = np.asarray(payload["ice"], dtype=np.float64)
        pdp = np.asarray(payload["pdp"], dtype=np.float64)
        if ice.ndim != 2 or pdp.ndim != 1 or ice.shape[1] != len(pdp):
            raise ValueError(f"ICE/PDP layout for {name!r} is inconsistent")
        deviation = np.mean(np.abs(ice - pdp[None, :]))
        rows.append(
            {
                "feature": str(name),
                "ice_mad": float(deviation),
                "ice_std": float(np.mean(np.std(ice, axis=0))),
            }
        )
    return pd.DataFrame(rows)


def ale_local_faithfulness(
    model, x, ale_dict, names: Sequence[str]
) -> pd.DataFrame:
    """Correlate ALE interval slopes with local finite differences of the model.

    For each ALE bin the mean prediction change when the feature is moved from
    the left to the right edge (holding other features) is compared with the
    ALE step. Higher Spearman correlation means a more faithful global curve.
    """
    features = _matrix(x)
    feature_names = _feature_names(names, features.shape[1])
    name_to_index = {name: index for index, name in enumerate(feature_names)}
    rows = []
    for name, curve in ale_dict.items():
        key = str(name)
        if key not in name_to_index:
            continue
        j = name_to_index[key]
        centers = np.asarray(curve["x"], dtype=np.float64)
        values = np.asarray(curve["ale"], dtype=np.float64)
        if len(centers) < 2:
            rows.append({"feature": key, "faithfulness": float("nan")})
            continue
        ale_steps = np.diff(values)
        model_steps = []
        for lo, hi in zip(centers[:-1], centers[1:]):
            mid = 0.5 * (lo + hi)
            span = 0.5 * abs(hi - lo) + 1e-12
            mask = np.abs(features[:, j] - mid) <= span
            if not np.any(mask):
                mask = np.ones(len(features), dtype=bool)
            left = features[mask].copy()
            right = features[mask].copy()
            left[:, j] = lo
            right[:, j] = hi
            model_steps.append(
                float(np.mean(_predict(model, right) - _predict(model, left)))
            )
        model_steps = np.asarray(model_steps, dtype=np.float64)
        ok = np.isfinite(ale_steps) & np.isfinite(model_steps)
        if (
            ok.sum() < 2
            or np.allclose(ale_steps[ok], ale_steps[ok][0])
            or np.allclose(model_steps[ok], model_steps[ok][0])
        ):
            score = float("nan")
        else:
            score = float(
                pd.Series(ale_steps[ok]).corr(
                    pd.Series(model_steps[ok]), method="spearman"
                )
            )
        rows.append({"feature": key, "faithfulness": score})
    return pd.DataFrame(rows)


def evaluate_shap_quality(
    model,
    x,
    shap_values,
    *,
    y=None,
    expected_value=None,
    explain_fn: Callable[[np.ndarray], np.ndarray] | None = None,
    k_neighbors: int = 5,
    top_k: int | None = None,
    n_faithfulness: int = 50,
    n_infidelity: int = 30,
    n_sensitivity: int = 40,
    sensitivity_n: int = 2,
    n_robust: int = 8,
    noise: float = 0.05,
    complexity_mass: float = 0.9,
    seed: int = 0,
) -> dict[str, float]:
    """Bundle faithfulness, locality, complexity and optional robustness for SHAP."""
    scores = {
        "additivity_rmse": shap_additivity(
            model, x, shap_values, expected_value=expected_value
        ),
        "faithfulness_correlation": shap_faithfulness_correlation(
            model, x, shap_values, n_perturb=n_faithfulness, seed=seed
        ),
        "infidelity": shap_infidelity(
            model,
            x,
            shap_values,
            n_perturb=n_infidelity,
            noise=noise,
            seed=seed + 1,
        ),
        "sensitivity_n": shap_sensitivity_n(
            model,
            x,
            shap_values,
            n=sensitivity_n,
            n_draws=n_sensitivity,
            seed=seed + 4,
        ),
        "sparseness": attribution_sparseness(shap_values),
    }
    scores.update(effective_complexity(shap_values, mass=complexity_mass))
    scores.update(
        shap_prediction_gaps(model, x, shap_values, top_k=top_k, seed=seed)
    )
    scores.update(
        explanation_locality(x, shap_values, k=k_neighbors, seed=seed + 2)
    )
    if y is not None:
        scores.update(explanation_contrastivity(shap_values, y))
    else:
        scores["contrastivity"] = float("nan")
        scores["contrast_cosine"] = float("nan")
        scores["n_low"] = float("nan")
        scores["n_high"] = float("nan")
    if explain_fn is not None:
        scores.update(
            shap_robustness(
                explain_fn,
                x,
                shap_values,
                n_perturb=n_robust,
                noise=noise,
                seed=seed + 3,
            )
        )
    else:
        scores["max_sensitivity"] = float("nan")
        scores["avg_sensitivity"] = float("nan")
        scores["relative_input_stability"] = float("nan")
    return scores


def evaluate_global_quality(
    model,
    x,
    *,
    ale_dict: dict | None = None,
    ice_pdp_dict: dict | None = None,
    shap_values=None,
    interaction_values=None,
    feature_names: Sequence[str] | None = None,
    n_bootstrap: int = 15,
    sample_fraction: float = 0.8,
    n_bins: int = 20,
    seed: int = 0,
) -> dict[str, pd.DataFrame]:
    """Faithfulness, robustness, agreement and interaction scores for global XAI."""
    out: dict[str, pd.DataFrame] = {}
    features = _matrix(x)
    names = _feature_names(feature_names, features.shape[1])
    if ale_dict is not None:
        out["ale_faithfulness"] = ale_local_faithfulness(
            model, features, ale_dict, names
        )
        out["ale_robustness"] = ale_robustness(
            model,
            features,
            ale_dict,
            names,
            n_bootstrap=n_bootstrap,
            sample_fraction=sample_fraction,
            seed=seed,
        )
    if ice_pdp_dict is not None:
        out["pdp_robustness"] = pdp_robustness(
            model,
            features,
            ice_pdp_dict,
            names,
            n_bootstrap=n_bootstrap,
            sample_fraction=sample_fraction,
            seed=seed + 1,
        )
        out["ice_heterogeneity"] = ice_heterogeneity(ice_pdp_dict)
    if (
        sum(item is not None for item in (ale_dict, ice_pdp_dict, shap_values))
        >= 2
    ):
        out["method_agreement"] = global_method_agreement(
            ale_dict=ale_dict,
            ice_pdp_dict=ice_pdp_dict,
            x=features,
            shap_values=shap_values,
            feature_names=names,
            n_bins=n_bins,
        )
    if interaction_values is not None:
        out["interaction_strength"] = interaction_strength(
            interaction_values, names
        )
    if not out:
        raise ValueError(
            "Provide ale_dict, ice_pdp_dict, shap_values and/or interaction_values"
        )
    return out
