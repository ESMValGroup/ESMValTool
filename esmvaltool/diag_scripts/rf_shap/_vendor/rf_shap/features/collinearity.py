from __future__ import annotations

from typing import Sequence

import numpy as np
from sklearn.linear_model import LinearRegression


def _corr_matrix(x: np.ndarray) -> np.ndarray:
    with np.errstate(invalid="ignore"):
        corr = np.corrcoef(x, rowvar=False)
    return np.nan_to_num(corr, nan=0.0)


def drop_by_correlation(x: np.ndarray, names: Sequence[str], threshold: float = 0.9) -> list[str]:
    corr = np.abs(_corr_matrix(x))
    np.fill_diagonal(corr, 0.0)
    keep = np.ones(x.shape[1], dtype=bool)
    order = np.argsort(-np.nanvar(x, axis=0))
    for i in order:
        if not keep[i]:
            continue
        conflicts = np.where(keep & (corr[i] > threshold))[0]
        keep[conflicts] = False
    return [names[i] for i, flag in enumerate(keep) if flag]


def vif_scores(x: np.ndarray) -> np.ndarray:
    n_features = x.shape[1]
    scores = np.full(n_features, np.inf)
    for i in range(n_features):
        others = np.delete(x, i, axis=1)
        if others.shape[1] == 0:
            scores[i] = 1.0
            continue
        model = LinearRegression()
        model.fit(others, x[:, i])
        pred = model.predict(others)
        ss_res = np.sum((x[:, i] - pred) ** 2)
        ss_tot = np.sum((x[:, i] - x[:, i].mean()) ** 2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
        scores[i] = np.inf if r2 >= 1 else 1.0 / (1.0 - r2)
    return scores


def drop_by_vif(x: np.ndarray, names: Sequence[str], threshold: float = 10.0) -> list[str]:
    keep = list(range(x.shape[1]))
    while len(keep) > 1:
        scores = vif_scores(x[:, keep])
        worst = int(np.argmax(scores))
        if scores[worst] <= threshold:
            break
        keep.pop(worst)
    return [names[i] for i in keep]


def drop_by_linreg(x: np.ndarray, names: Sequence[str], r2_threshold: float = 0.95) -> list[str]:
    """Drop a feature if it is almost perfectly predicted by the others."""
    keep = list(range(x.shape[1]))
    changed = True
    while changed and len(keep) > 1:
        changed = False
        block = x[:, keep]
        r2s = []
        for i in range(block.shape[1]):
            others = np.delete(block, i, axis=1)
            model = LinearRegression().fit(others, block[:, i])
            pred = model.predict(others)
            ss_res = np.sum((block[:, i] - pred) ** 2)
            ss_tot = np.sum((block[:, i] - block[:, i].mean()) ** 2)
            r2s.append(1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0)
        worst = int(np.argmax(r2s))
        if r2s[worst] >= r2_threshold:
            keep.pop(worst)
            changed = True
    return [names[i] for i in keep]


def drop_collinear(
    x: np.ndarray,
    names: Sequence[str],
    corr_threshold: float = 0.9,
    vif_threshold: float = 10.0,
    linreg_r2_threshold: float = 0.95,
) -> list[str]:
    kept = drop_by_correlation(x, names, corr_threshold)
    idx = [names.index(n) for n in kept]
    kept = drop_by_vif(x[:, idx], kept, vif_threshold)
    idx = [names.index(n) for n in kept]
    kept = drop_by_linreg(x[:, idx], kept, linreg_r2_threshold)
    return kept
