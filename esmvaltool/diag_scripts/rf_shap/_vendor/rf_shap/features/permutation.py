from __future__ import annotations

from typing import Sequence

import numpy as np
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.inspection import permutation_importance


def permutation_elimination(
    x: np.ndarray,
    y: np.ndarray,
    names: Sequence[str],
    n_repeats: int = 5,
    drop_below: float = 0.0,
    max_rounds: int = 20,
    random_state: int = 42,
) -> list[str]:
    """Iteratively drop the feature with the lowest permutation importance."""
    keep = list(range(x.shape[1]))
    rng = random_state
    for _ in range(max_rounds):
        if len(keep) <= 1:
            break
        model = ExtraTreesRegressor(
            n_estimators=160,
            max_depth=12,
            n_jobs=-1,
            random_state=rng,
        )
        xx = x[:, keep]
        model.fit(xx, y)
        result = permutation_importance(
            model,
            xx,
            y,
            n_repeats=n_repeats,
            random_state=rng,
            n_jobs=-1,
        )
        scores = result.importances_mean
        worst = int(np.argmin(scores))
        if scores[worst] >= drop_below and np.all(scores > drop_below):
            break
        keep.pop(worst)
    return [names[i] for i in keep]
