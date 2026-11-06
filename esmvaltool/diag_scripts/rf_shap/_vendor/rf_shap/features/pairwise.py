from __future__ import annotations

from typing import Sequence

import numpy as np
from sklearn.ensemble import ExtraTreesRegressor


def _importance(x: np.ndarray, y: np.ndarray, random_state: int = 42) -> np.ndarray:
    model = ExtraTreesRegressor(
        n_estimators=120,
        max_depth=10,
        n_jobs=-1,
        random_state=random_state,
    )
    model.fit(x, y)
    return np.asarray(model.feature_importances_)


def drop_weaker_of_pairs(
    x: np.ndarray,
    y: np.ndarray,
    names: Sequence[str],
    corr_threshold: float = 0.8,
    random_state: int = 42,
) -> list[str]:
    """If two features are highly correlated, drop the less important one."""
    names = list(names)
    if len(names) <= 1:
        return names
    corr = np.corrcoef(x, rowvar=False)
    corr = np.nan_to_num(np.abs(corr), nan=0.0)
    np.fill_diagonal(corr, 0.0)
    importance = _importance(x, y, random_state)
    keep = np.ones(len(names), dtype=bool)
    pairs = np.argwhere(np.triu(corr, k=1) > corr_threshold)
    pairs = sorted(pairs, key=lambda p: -corr[p[0], p[1]])
    for i, j in pairs:
        if not (keep[i] and keep[j]):
            continue
        drop = j if importance[i] >= importance[j] else i
        keep[drop] = False
    return [n for n, flag in zip(names, keep) if flag]
