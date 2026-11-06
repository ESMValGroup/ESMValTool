"""Optional wrapper around ThrunGroup FastForest (MABSplit).

The repository is research code, not a NumPy-2 wheel. Clone
https://github.com/ThrunGroup/FastForest and put that directory on
``PYTHONPATH``. Do not install its ``requirements.txt`` into the
``esmvaltool`` environment: those pins replace NumPy, scikit-learn, and
pandas with versions that break ESMValTool and ``fast_tree_shap``.
"""

from __future__ import annotations

from typing import Any

import numpy as np

_IGNORED_SKLEARN_KEYS = {
    "n_jobs",
    "max_samples",
    "max_features",
    "min_samples_leaf",
    "bootstrap",
    "oob_score",
    "ccp_alpha",
    "warm_start",
    "criterion",
    "max_leaf_nodes",
}


def _load_fastforest():
    try:
        from data_structures.wrappers.histogram_random_forest_regressor import (
            HistogramRandomForestRegressor,
        )
        from utils.constants import MAB, MSE
    except ImportError as exc:
        raise ImportError(
            "fast_forest needs ThrunGroup/FastForest on PYTHONPATH "
            "(https://github.com/ThrunGroup/FastForest). "
            "Clone the repository and export PYTHONPATH to its root. "
            "Do not pip-install its requirements.txt into the esmvaltool environment."
        ) from exc
    return HistogramRandomForestRegressor, MAB, MSE


class FastForestRegressor:
    """sklearn-style fit/predict around histogram FastForest with MABSplit.

    ``supports_tree_shap`` is false: these trees are not scikit-learn,
    XGBoost, or LightGBM trees, so ``fast_tree_shap`` and ``shap.TreeExplainer``
    do not apply. ``compute_shap(..., method="auto")`` then uses a
    model-agnostic explainer. sklearn ``rf`` and ``fast_tree_shap`` are unchanged.
    """

    supports_tree_shap = False

    def __init__(
        self,
        n_estimators: int = 100,
        max_depth: int | None = 16,
        min_samples_split: int = 2,
        num_bins: int = 11,
        solver: str = "MAB",
        random_state: int = 0,
        budget: int | None = None,
        verbose: bool = False,
    ) -> None:
        self.n_estimators = int(n_estimators)
        self.max_depth = None if max_depth is None else int(max_depth)
        self.min_samples_split = int(min_samples_split)
        self.num_bins = int(num_bins)
        self.solver = solver
        self.random_state = int(random_state)
        self.budget = budget
        self.verbose = bool(verbose)
        self._model: Any = None
        self.n_features_in_: int | None = None

    def fit(self, x, y):
        forest_cls, mab, mse = _load_fastforest()
        features = np.asarray(x, dtype=np.float64)
        target = np.asarray(y, dtype=np.float64).reshape(-1)
        if features.ndim != 2:
            raise ValueError("X must be 2D (n_samples, n_features)")
        solver = mab if str(self.solver).upper() == "MAB" else str(self.solver)
        self.n_features_in_ = features.shape[1]
        self._model = forest_cls(
            n_estimators=self.n_estimators,
            max_depth=self.max_depth,
            min_samples_split=self.min_samples_split,
            num_bins=self.num_bins,
            solver=solver,
            criterion=mse,
            random_state=self.random_state,
            budget=self.budget,
            verbose=self.verbose,
        )
        self._model.fit(features, target)
        return self

    def predict(self, x) -> np.ndarray:
        if self._model is None:
            raise RuntimeError("FastForestRegressor is not fitted")
        features = np.asarray(x, dtype=np.float64)
        if features.ndim == 1:
            features = features.reshape(1, -1)
        return np.asarray(self._model.predict_batch(features), dtype=np.float64)

    @property
    def feature_importances_(self) -> np.ndarray:
        if self._model is None:
            raise RuntimeError("FastForestRegressor is not fitted")
        importance = np.asarray(self._model.calculate_mdi(), dtype=np.float64)
        total = float(importance.sum())
        if total > 0:
            importance = importance / total
        return importance


def fast_forest(
    n_estimators: int = 100,
    max_depth: int | None = 16,
    min_samples_split: int = 2,
    num_bins: int = 11,
    solver: str = "MAB",
    random_state: int = 0,
    budget: int | None = None,
    verbose: bool = False,
    **unused,
) -> FastForestRegressor:
    """Build an unfitted FastForest regressor.

    This is the fast training option from ThrunGroup/FastForest (Tiwari et al.,
    NeurIPS 2022, MABSplit). It bins each feature and allocates split
    evaluations with a multi-armed bandit, instead of scanning every sample
    at every candidate split. The published speedups are about the training
    step, up to roughly 100x in some of their settings.

    Extra scikit-learn keys such as ``n_jobs`` are ignored. Calling this
    raises ``ImportError`` until the FastForest repository is on ``PYTHONPATH``.
    """
    _load_fastforest()
    unused = {key: value for key, value in unused.items() if key not in _IGNORED_SKLEARN_KEYS}
    if unused:
        unknown = ", ".join(sorted(unused))
        raise TypeError(f"fast_forest() got unexpected arguments: {unknown}")
    return FastForestRegressor(
        n_estimators=n_estimators,
        max_depth=max_depth,
        min_samples_split=min_samples_split,
        num_bins=num_bins,
        solver=solver,
        random_state=random_state,
        budget=budget,
        verbose=verbose,
    )
