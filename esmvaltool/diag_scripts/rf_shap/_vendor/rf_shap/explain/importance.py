from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd
from sklearn.inspection import permutation_importance


def model_importance(model, names: Sequence[str]) -> pd.DataFrame:
    if not hasattr(model, "feature_importances_"):
        raise AttributeError("This model does not expose feature_importances_")
    values = np.asarray(model.feature_importances_)
    return pd.DataFrame({"feature": list(names), "importance": values}).sort_values("importance", ascending=False)


def permutation_scores(model, x, y, names: Sequence[str], n_repeats: int = 5, seed: int = 42) -> pd.DataFrame:
    result = permutation_importance(model, x, y, n_repeats=n_repeats, random_state=seed, n_jobs=-1)
    return pd.DataFrame(
        {
            "feature": list(names),
            "importance": result.importances_mean,
            "std": result.importances_std,
        }
    ).sort_values("importance", ascending=False)
