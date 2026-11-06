"""Helpers shared by the ESMValTool diagnostic and the local runner."""

from __future__ import annotations

from pathlib import Path

import xarray as xr

from rf_shap.config import PipelineConfig
from rf_shap.cube.align import align_xy
from rf_shap.pipeline import DiagnosisPipeline, PipelineResult
from rf_shap.plot import (
    plot_ale,
    plot_error_map,
    plot_ice_pdp,
    plot_importance_bars,
    plot_pred_scatter,
    plot_shap_maps,
    plot_shap_summary,
)
from rf_shap.predict.reconstruct import reconstruct_prediction

_LAT = ("lat", "latitude", "nav_lat")
_LON = ("lon", "longitude", "nav_lon")
_TIME = ("time", "valid_time")


def rename_axes(ds: xr.Dataset) -> xr.Dataset:
    mapping = {}
    for names, dest in ((_TIME, "time"), (_LAT, "lat"), (_LON, "lon")):
        for name in names:
            if name in ds.dims or name in ds.coords:
                if name != dest:
                    mapping[name] = dest
                break
    if mapping:
        ds = ds.rename(mapping)
    return ds


def open_field(path: str | Path, name: str | None = None) -> xr.DataArray:
    ds = xr.open_dataset(path)
    ds = rename_axes(ds)
    if name and name in ds.data_vars:
        da = ds[name]
    else:
        data_vars = [v for v in ds.data_vars if {"lat", "lon"} <= set(ds[v].dims)]
        if not data_vars:
            raise ValueError(f"No lat/lon variable in {path}")
        da = ds[data_vars[0]]
    extra = [d for d in da.dims if d not in ("time", "lat", "lon")]
    if extra:
        da = da.isel({d: 0 for d in extra})
    if "time" in da.dims:
        da = da.transpose("time", "lat", "lon", missing_dims="ignore")
    else:
        da = da.transpose("lat", "lon", missing_dims="ignore")
    if name:
        da = da.rename(name)
    return da.load()


def dataset_from_fields(fields: dict[str, xr.DataArray]) -> xr.Dataset:
    ds = xr.Dataset(fields)
    ds, *_ = xr.align(ds, join="inner", copy=False)
    return ds


def load_xy_from_files(x_files: dict, y_file: str | Path, y_name: str = "y") -> tuple[xr.Dataset, xr.Dataset]:
    x_fields = {name: open_field(path, name) for name, path in dict(x_files).items()}
    y = xr.Dataset({y_name: open_field(y_file, y_name)})
    return dataset_from_fields(x_fields), y


def default_pipeline_cfg() -> dict:
    return {
        "seed": 42,
        "x_transform": {"mode": "raw"},
        "preprocess": {
            "x": {"deseasonalize": False, "detrend": False},
            "y": {"deseasonalize": False, "detrend": False},
        },
        "feature_selection": {"methods": ["collinearity", "pairwise_importance"]},
        "split": {"method": "time", "test_size": 0.2},
        "sampling": {"method": "random", "n_samples": 8000},
        "model": {
            "name": "extra_trees",
            "params": {"n_estimators": 120, "max_depth": 12, "n_jobs": -1},
        },
        "predict": {"restore_cube": True},
        "explain": {
            "rf_importance": True,
            "permutation_importance": True,
            "permutation_n_repeats": 3,
            "shap": {
                "enabled": True,
                "method": "auto",
                "max_samples": 1500,
                "background_size": 80,
            },
            "ale": True,
            "ice_pdp": True,
            "ice_n_samples": 40,
            "pdp_grid": 15,
        },
    }


def run_pipeline_on_cubes(
    x: xr.Dataset,
    y: xr.Dataset,
    override: dict | None = None,
    regrid: str = "nearest",
    join_time: str = "inner",
    align_target: str = "y",
) -> PipelineResult:
    pipe_cfg = PipelineConfig.from_dict(default_pipeline_cfg())
    if override:
        pipe_cfg = pipe_cfg.updated(override)
    x, y = align_xy(x, y, regrid=regrid, join_time=join_time, target=align_target)
    result = DiagnosisPipeline(pipe_cfg).run(x=x, y=y)
    if isinstance(result, list):
        return result[0]
    return result


def save_result_files(result: PipelineResult, work_dir: str | Path, plot_dir: str | Path | None = None) -> None:
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    plots = Path(plot_dir or work / "plots")
    plots.mkdir(parents=True, exist_ok=True)

    import pandas as pd

    pd.Series(result.metrics).to_frame("value").to_csv(work / "metrics.csv")
    if result.prediction_cube is not None:
        result.prediction_cube.to_dataset(name="prediction").to_netcdf(work / "prediction.nc")
        try:
            obs_cube = reconstruct_prediction(result.split.test.y, result.split.test, name="obs")
            plot_error_map(obs_cube, result.prediction_cube).savefig(plots / "prediction_maps.png")
        except Exception:
            pass
    if result.shap_cube is not None:
        result.shap_cube.to_netcdf(work / "shap.nc")
        try:
            plot_shap_maps(result.shap_cube).savefig(plots / "shap_maps.png")
        except Exception:
            pass
        shap_res = result.artifacts.get("shap_result")
        if shap_res is not None:
            xs = result.split.test.X[shap_res.sample_positions]
            plot_shap_summary(shap_res.values, xs, shap_res.feature_names).savefig(plots / "shap_summary.png")
    if "rf_importance" in result.artifacts:
        result.artifacts["rf_importance"].to_csv(work / "rf_importance.csv", index=False)
        plot_importance_bars(result.artifacts["rf_importance"], title="tree importance").savefig(
            plots / "rf_importance.png"
        )
    if "permutation_importance" in result.artifacts:
        result.artifacts["permutation_importance"].to_csv(work / "permutation_importance.csv", index=False)
        plot_importance_bars(
            result.artifacts["permutation_importance"],
            title="permutation importance",
        ).savefig(plots / "permutation_importance.png")
    if "ale" in result.artifacts:
        try:
            plot_ale(result.artifacts["ale"]).savefig(plots / "ale.png")
        except Exception:
            pass
    if "ice_pdp" in result.artifacts:
        try:
            plot_ice_pdp(result.artifacts["ice_pdp"]).savefig(plots / "ice_pdp.png")
        except Exception:
            pass
    plot_pred_scatter(result.split.test.y, result.prediction, result.metrics).savefig(plots / "pred_vs_obs.png")
    import matplotlib.pyplot as plt

    plt.close("all")
