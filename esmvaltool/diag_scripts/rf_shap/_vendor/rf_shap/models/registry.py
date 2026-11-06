from __future__ import annotations

from typing import Any

from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor, RandomForestRegressor

from rf_shap.models.fast_forest import fast_forest


def _optional_import(name: str):
    try:
        return __import__(name)
    except ImportError as exc:
        raise ImportError(f"Model backend '{name}' is not installed") from exc


def available_models() -> list[str]:
    names = ["rf", "extra_trees", "hist_gbm"]
    for pkg, alias in [("xgboost", "xgboost"), ("lightgbm", "lightgbm")]:
        try:
            __import__(pkg)
            names.append(alias)
            if pkg == "lightgbm":
                names.append("lightgbm_rf")
        except ImportError:
            pass
    try:
        __import__("data_structures.wrappers.histogram_random_forest_regressor")
        names.append("fastforest")
    except ImportError:
        pass
    return names


def build_model(name: str = "extra_trees", params: dict[str, Any] | None = None, seed: int = 42):
    """Tree-model factory, including faster RF-style backends.

    Fast RF options (inspired by Kursa & Piwoński 2026 SoftwareX `fru`,
    plus common Python substitutes):
      extra_trees   : ExtraTrees, usually much faster than sklearn RF
      hist_gbm      : histogram boosting, very fast on large cubes
      lightgbm_rf   : LightGBM random-forest boosting mode
      rf            : sklearn RandomForest with multi-core defaults
      fastforest    : ThrunGroup FastForest (MABSplit), only if that repo is on PYTHONPATH
    """
    params = dict(params or {})
    name = name.lower()
    params.setdefault("random_state", seed)

    if name == "rf":
        params.setdefault("n_estimators", 300)
        params.setdefault("n_jobs", -1)
        params.setdefault("max_depth", 16)
        params.setdefault("min_samples_leaf", 2)
        params.setdefault("max_samples", 0.7)
        return RandomForestRegressor(**params)

    if name in {"extra_trees", "fast_rf", "et"}:
        params.setdefault("n_estimators", 300)
        params.setdefault("n_jobs", -1)
        params.setdefault("max_depth", 16)
        params.setdefault("min_samples_leaf", 2)
        return ExtraTreesRegressor(**params)

    if name in {"hist_gbm", "hgb"}:
        params.pop("n_jobs", None)
        params.setdefault("max_depth", 8)
        params.setdefault("learning_rate", 0.08)
        params.setdefault("max_iter", 300)
        return HistGradientBoostingRegressor(**params)

    if name == "xgboost":
        xgb = _optional_import("xgboost")
        params.setdefault("n_estimators", 400)
        params.setdefault("max_depth", 8)
        params.setdefault("n_jobs", -1)
        params.setdefault("tree_method", "hist")
        params.setdefault("learning_rate", 0.08)
        return xgb.XGBRegressor(**params)

    if name == "lightgbm":
        lgb = _optional_import("lightgbm")
        params.setdefault("n_estimators", 400)
        params.setdefault("n_jobs", -1)
        params.setdefault("learning_rate", 0.08)
        return lgb.LGBMRegressor(**params)

    if name == "lightgbm_rf":
        lgb = _optional_import("lightgbm")
        params.setdefault("boosting_type", "rf")
        params.setdefault("bagging_freq", 1)
        params.setdefault("bagging_fraction", 0.8)
        params.setdefault("feature_fraction", 0.8)
        params.setdefault("n_estimators", 400)
        params.setdefault("n_jobs", -1)
        params.pop("learning_rate", None)
        return lgb.LGBMRegressor(**params)

    if name in {"fastforest", "fast_forest", "mabforest"}:
        return fast_forest(**params)

    raise ValueError(f"Unknown model '{name}'. Available: {available_models()}")
