from __future__ import annotations

from typing import Sequence

import numpy as np

from rf_shap.cube.tabular import SampleTable
from rf_shap.features.collinearity import drop_collinear
from rf_shap.features.pairwise import drop_weaker_of_pairs
from rf_shap.features.permutation import permutation_elimination
from rf_shap.features.rfe import recursive_feature_elimination


def select_features(table: SampleTable, options: dict | None = None, seed: int = 42) -> SampleTable:
    options = options or {}
    methods: Sequence[str] = options.get("methods") or []
    names = list(table.feature_names)
    x, y = table.X, table.y
    history = []
    for method in methods:
        if len(names) <= 1:
            break
        idx = [table.feature_names.index(n) for n in names]
        xx = table.X[:, idx]
        if method == "collinearity":
            cfg = options.get("collinearity", {})
            names = drop_collinear(
                xx,
                names,
                corr_threshold=cfg.get("corr_threshold", 0.9),
                vif_threshold=cfg.get("vif_threshold", 10.0),
                linreg_r2_threshold=cfg.get("linreg_r2_threshold", 0.95),
            )
        elif method == "pairwise_importance":
            cfg = options.get("pairwise_importance", {})
            names = drop_weaker_of_pairs(
                xx,
                y,
                names,
                corr_threshold=cfg.get("corr_threshold", 0.8),
                random_state=seed,
            )
        elif method == "permutation":
            cfg = options.get("permutation", {})
            names = permutation_elimination(
                xx,
                y,
                names,
                n_repeats=int(cfg.get("n_repeats", 5)),
                drop_below=float(cfg.get("drop_below", 0.0)),
                max_rounds=int(cfg.get("max_rounds", 20)),
                random_state=seed,
            )
        elif method == "rfe":
            cfg = options.get("rfe", {})
            names = recursive_feature_elimination(
                xx,
                y,
                names,
                n_features_to_select=cfg.get("n_features_to_select"),
                step=int(cfg.get("step", 1)),
                cv=cfg.get("cv", 3),
                random_state=seed,
            )
        else:
            raise ValueError(f"Unknown feature selection method: {method}")
        history.append({"method": method, "kept": list(names)})
    selected = table.select_columns(names)
    selected.extras["selection_history"] = history
    return selected
