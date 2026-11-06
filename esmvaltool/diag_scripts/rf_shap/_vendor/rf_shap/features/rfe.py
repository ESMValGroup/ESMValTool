from __future__ import annotations

from typing import Sequence

import numpy as np
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.feature_selection import RFE, RFECV
from sklearn.model_selection import KFold


def recursive_feature_elimination(
    x: np.ndarray,
    y: np.ndarray,
    names: Sequence[str],
    n_features_to_select: int | None = None,
    step: int = 1,
    cv: int | None = 3,
    random_state: int = 42,
) -> list[str]:
    estimator = ExtraTreesRegressor(
        n_estimators=120,
        max_depth=10,
        n_jobs=-1,
        random_state=random_state,
    )
    n_features = x.shape[1]
    if n_features <= 1:
        return list(names)
    if cv and n_features_to_select is None:
        selector = RFECV(
            estimator,
            step=step,
            cv=KFold(n_splits=min(cv, max(2, len(y) // 10)), shuffle=True, random_state=random_state),
            scoring="r2",
            n_jobs=-1,
        )
    else:
        k = n_features_to_select or max(1, n_features // 2)
        selector = RFE(estimator, n_features_to_select=k, step=step)
    selector.fit(x, y)
    mask = selector.support_
    return [n for n, flag in zip(names, mask) if flag]
