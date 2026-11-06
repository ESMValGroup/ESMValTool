from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import xarray as xr

from rf_shap.config import PipelineConfig, load_config
from rf_shap.cube.align import align_xy
from rf_shap.cube.features import build_x_features, pca_transform_table
from rf_shap.cube.io import load_cube
from rf_shap.cube.preprocess import preprocess_dataset
from rf_shap.cube.tabular import SampleTable, cube_to_table
from rf_shap.explain.ale import ale_1d
from rf_shap.explain.ice_pdp import ice_pdp
from rf_shap.explain.importance import model_importance, permutation_scores
from rf_shap.explain.shap_explain import compute_shap, shap_to_cube
from rf_shap.features.select import select_features
from rf_shap.models.registry import build_model
from rf_shap.models.sampling import sample_train
from rf_shap.models.split import SplitResult, split_table
from rf_shap.predict.metrics import regression_metrics
from rf_shap.predict.reconstruct import reconstruct_prediction


@dataclass
class PipelineResult:
    metrics: dict
    model: Any
    table: SampleTable
    split: SplitResult
    prediction: np.ndarray
    prediction_cube: xr.DataArray | None = None
    shap_cube: xr.Dataset | None = None
    artifacts: dict = field(default_factory=dict)

    def save(self, out_dir: str | Path) -> Path:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        pd.Series(self.metrics).to_csv(out / "metrics.csv")
        if self.prediction_cube is not None:
            self.prediction_cube.to_netcdf(out / "prediction.nc")
        if self.shap_cube is not None:
            self.shap_cube.to_netcdf(out / "shap.nc")
        if "rf_importance" in self.artifacts:
            self.artifacts["rf_importance"].to_csv(out / "rf_importance.csv", index=False)
        if "permutation_importance" in self.artifacts:
            self.artifacts["permutation_importance"].to_csv(out / "permutation_importance.csv", index=False)
        if "shap_table" in self.artifacts:
            self.artifacts["shap_table"].to_csv(out / "shap_samples.csv", index=False)
        return out


class DiagnosisPipeline:
    def __init__(self, config: PipelineConfig):
        self.config = config

    def load_xy(self) -> tuple[xr.Dataset, xr.Dataset]:
        data_cfg = self.config.raw.get("data", {})
        x_cfg = data_cfg.get("x", {})
        y_cfg = data_cfg.get("y", {})
        common = dict(
            time_name=data_cfg.get("time_name", "time"),
            lat_name=data_cfg.get("lat_name", "lat"),
            lon_name=data_cfg.get("lon_name", "lon"),
        )
        x = load_cube(x_cfg["path"], x_cfg.get("variables"), **common)
        y = load_cube(y_cfg["path"], y_cfg.get("variables"), **common)
        return align_xy(
            x,
            y,
            regrid=data_cfg.get("regrid", "nearest"),
            join_time=data_cfg.get("join_time", "inner"),
            target=data_cfg.get("target", "y"),
        )

    def prepare_table(self, x: xr.Dataset, y: xr.Dataset) -> SampleTable:
        prep_cfg = self.config.raw.get("preprocess", {})
        x = preprocess_dataset(x, prep_cfg.get("x"))
        y = preprocess_dataset(y, prep_cfg.get("y"))
        xt_cfg = self.config.raw.get("x_transform", {})
        x, meta = build_x_features(
            x,
            mode=xt_cfg.get("mode", "raw"),
            lags=xt_cfg.get("lags", [1, 2]),
            pca_n_components=xt_cfg.get("pca_n_components", 0.95),
            pca_whiten=bool(xt_cfg.get("pca_whiten", False)),
        )
        y, x = xr.align(y, x, join="inner")
        table = cube_to_table(x, y)
        table.extras["x_transform"] = meta
        if xt_cfg.get("mode") in {"pca", "lag_pca"} and xt_cfg.get("feature_pca", True):
            # additional table-level PCA is optional; cube EOF already applied
            pass
        if xt_cfg.get("table_pca"):
            z, names, model = pca_transform_table(
                table.X,
                n_components=xt_cfg.get("pca_n_components", 0.95),
                whiten=bool(xt_cfg.get("pca_whiten", False)),
                random_state=self.config.seed,
            )
            table = table.with_features(z, names)
            table.extras["pca_model"] = model
        return table

    def _fit_one(self, split: SplitResult) -> PipelineResult:
        seed = self.config.seed
        train = select_features(split.train, self.config.raw.get("feature_selection"), seed=seed)
        test = split.test.select_columns(train.feature_names)
        samp = self.config.raw.get("sampling", {})
        names = list(train.index.names)
        idx = sample_train(
            train.X,
            train.y,
            method=samp.get("method", "random"),
            n_samples=samp.get("n_samples"),
            n_bins=int(samp.get("n_bins", 10)),
            seed=seed,
            lat=np.asarray(train.index.get_level_values("lat")) if "lat" in names else None,
            lon=np.asarray(train.index.get_level_values("lon")) if "lon" in names else None,
            time=np.asarray(train.index.get_level_values("time")) if "time" in names else None,
            area_weights=samp.get("area_weights"),
            time_weights=samp.get("time_weights"),
            strata=samp.get("strata"),
            stratum_columns=samp.get("stratum_columns"),
            categorical_columns=samp.get("categorical_columns"),
            categorical=samp.get("categorical"),
            balance=samp.get("balance", "equal"),
            n_clusters=samp.get("n_clusters"),
            extreme_fraction=float(samp.get("extreme_fraction", 0.5)),
            extreme_quantile=float(samp.get("extreme_quantile", 0.9)),
            extreme_side=samp.get("extreme_side", "both"),
            threshold=samp.get("threshold"),
            threshold_band=samp.get("threshold_band"),
            threshold_column=samp.get("threshold_column"),
            n_iter=int(samp.get("n_iter", 2000)),
            active_rounds=int(samp.get("active_rounds", 3)),
            active_batch=samp.get("active_batch"),
            dependence=samp.get("dependence", "auto"),
            objective_weights=samp.get("objective_weights"),
        )
        train_s = train.take(idx)
        model_cfg = self.config.raw.get("model", {})
        model = build_model(model_cfg.get("name", "extra_trees"), model_cfg.get("params"), seed=seed)
        model.fit(train_s.X, train_s.y)
        pred = model.predict(test.X)
        metrics = regression_metrics(test.y, pred)
        metrics["n_train"] = int(len(train_s.y))
        metrics["n_test"] = int(len(test.y))
        metrics["n_features"] = int(len(train.feature_names))
        pred_cube = None
        template = train.extras.get("template", split.train.extras.get("template"))
        if template is not None:
            train.extras["template"] = template
            test.extras["template"] = template
        if self.config.get("predict", "restore_cube", default=True) and template is not None:
            pred_cube = reconstruct_prediction(pred, test)
        artifacts: dict[str, Any] = {
            "features": list(train.feature_names),
            "selection": train.extras.get("selection_history"),
            "sampled_index": train_s.index,
        }
        explain_cfg = self.config.raw.get("explain", {})
        if explain_cfg.get("rf_importance", True):
            try:
                artifacts["rf_importance"] = model_importance(model, train.feature_names)
            except AttributeError:
                pass
        if explain_cfg.get("permutation_importance", True):
            artifacts["permutation_importance"] = permutation_scores(
                model,
                test.X,
                test.y,
                train.feature_names,
                n_repeats=int(explain_cfg.get("permutation_n_repeats", 5)),
                seed=seed,
            )
        shap_cube = None
        shap_cfg = explain_cfg.get("shap", {})
        if shap_cfg.get("enabled", True):
            shap_table = test
            if template is not None:
                shap_table.extras["template"] = template
            shap_res = compute_shap(
                model,
                shap_table,
                method=shap_cfg.get("method", "auto"),
                algorithm=shap_cfg.get("algorithm", "v2"),
                max_samples=shap_cfg.get("max_samples", 2000),
                background_size=int(shap_cfg.get("background_size", 100)),
                n_jobs=int(shap_cfg.get("n_jobs", -1)),
                seed=seed,
                background=train_s.X,
                feature_perturbation=shap_cfg.get("feature_perturbation", "tree_path_dependent"),
            )
            artifacts["shap_result"] = shap_res
            artifacts["shap_table"] = shap_res.to_frame()
            shap_cube = shap_to_cube(shap_res, shap_table)
        if explain_cfg.get("ale", True):
            artifacts["ale"] = ale_1d(model, test.X, train.feature_names)
        if explain_cfg.get("ice_pdp", True):
            artifacts["ice_pdp"] = ice_pdp(
                model,
                test.X,
                train.feature_names,
                n_grid=int(explain_cfg.get("pdp_grid", 20)),
                ice_n_samples=int(explain_cfg.get("ice_n_samples", 50)),
                seed=seed,
            )
        used_split = SplitResult(train_s, test, split.method, split.extras)
        return PipelineResult(metrics, model, train, used_split, pred, pred_cube, shap_cube, artifacts)

    def run(self, x: xr.Dataset | None = None, y: xr.Dataset | None = None) -> PipelineResult | list[PipelineResult]:
        if x is None or y is None:
            x, y = self.load_xy()
        table = self.prepare_table(x, y)
        split = split_table(table, self.config.raw.get("split"), seed=self.config.seed)
        if isinstance(split, list):
            results = [self._fit_one(s) for s in split]
            return results
        return self._fit_one(split)


def run_pipeline(config_path=None, override=None, x=None, y=None):
    cfg = load_config(config_path, override)
    result = DiagnosisPipeline(cfg).run(x=x, y=y)
    out_dir = cfg.output_dir
    if isinstance(result, list):
        for i, item in enumerate(result):
            item.save(out_dir / f"fold_{i+1}")
    else:
        result.save(out_dir)
    return result
