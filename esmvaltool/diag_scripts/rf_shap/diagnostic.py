"""ESMValTool diagnostic wrapper around the rf_shap pipeline.

The preprocessor writes one NetCDF per variable. This script groups those
files into predictor cubes X and a target cube Y, runs DiagnosisPipeline,
and writes prediction/SHAP cubes plus figures into the ESMValTool work and
plot directories.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from esmvaltool.diag_scripts.rf_shap._bootstrap import ensure_rf_shap_on_path

ensure_rf_shap_on_path()

# isort: off
from esmvaltool.diag_scripts.rf_shap.compare import (  # noqa: E402
    compare_pipeline_results,
)
from esmvaltool.diag_scripts.shared import (  # noqa: E402
    ProvenanceLogger,
    get_diagnostic_filename,
    run_diagnostic,
    save_figure,
    select_metadata,
)
from rf_shap.esmval_bridge import (  # noqa: E402
    dataset_from_fields,
    load_xy_from_files,
    open_field,
    run_pipeline_on_cubes,
)
from rf_shap.plot import (  # noqa: E402
    plot_ale,
    plot_ale_overlay,
    plot_comparison_bars,
    plot_cube_difference,
    plot_error_map,
    plot_ice_pdp,
    plot_importance_bars,
    plot_importance_overlay,
    plot_metrics_compare,
    plot_pdp_overlay,
    plot_pred_scatter,
    plot_shap_dependence,
    plot_shap_dependence_overlay,
    plot_shap_maps,
    plot_shap_summary,
    write_three_column_html,
)
from rf_shap.predict.reconstruct import reconstruct_prediction  # noqa: E402
from rf_shap.synthetic import make_synthetic_cubes  # noqa: E402
# isort: on

logger = logging.getLogger(Path(__file__).stem)


def _provenance(_cfg, ancestors, caption: str) -> dict:
    """Build an ESMValTool provenance record for one diagnostic product."""
    return {
        "caption": caption,
        "authors": ["huang_feini"],
        "references": ["acknow_project"],
        "ancestors": list(ancestors),
        "domains": ["global"],
        "plot_types": ["other"],
        "statistics": ["other"],
    }


def _save_nc(basename: str, dataset, cfg, provenance: dict) -> Path:
    path = Path(get_diagnostic_filename(basename, cfg, extension="nc"))
    path.parent.mkdir(parents=True, exist_ok=True)
    encoding = {
        name: {"dtype": "float32", "zlib": True, "complevel": 1, "_FillValue": None}
        for name in dataset.data_vars
    }
    tmp = path.with_suffix(path.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()
    dataset.to_netcdf(tmp, encoding=encoding)
    tmp.replace(path)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(str(path), provenance)
    logger.info("Wrote %s", path)
    return path


def _save_csv(basename: str, frame, cfg, provenance: dict) -> Path:
    path = Path(get_diagnostic_filename(basename, cfg, extension="csv"))
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(str(path), provenance)
    logger.info("Wrote %s", path)
    return path


def load_xy_from_esmvaltool(cfg):
    """Build X and Y cubes from preprocessor output listed in cfg['input_data']."""
    records = list(cfg["input_data"].values())
    predictors = list(cfg.get("predictors") or [])
    target = cfg.get("target")
    if not target:
        msg = "Recipe script must set 'target' (Y short_name or variable_group)"
        raise ValueError(msg)
    if not predictors:
        names = sorted({r.get("short_name") for r in records})
        predictors = [n for n in names if n != target]
    if not predictors:
        msg = "No predictor variables left after removing the target"
        raise ValueError(msg)

    pred_dataset = cfg.get("predictor_dataset")
    target_dataset = cfg.get("target_dataset")

    x_fields = {}
    ancestors = []
    for name in predictors:
        hits = select_metadata(records, short_name=name) or select_metadata(
            records,
            variable_group=name,
        )
        if pred_dataset:
            hits = [h for h in hits if h.get("dataset") == pred_dataset]
        if not hits:
            msg = f"No preprocessed file for predictor {name!r}"
            raise KeyError(msg)
        item = hits[0]
        x_fields[name] = open_field(item["filename"], name)
        ancestors.append(item["filename"])

    y_hits = select_metadata(records, short_name=target) or select_metadata(
        records,
        variable_group=target,
    )
    if target_dataset:
        y_hits = [h for h in y_hits if h.get("dataset") == target_dataset]
    if not y_hits:
        msg = f"No preprocessed file for target {target!r}"
        raise KeyError(msg)
    y_item = y_hits[0]
    y_name = str(cfg.get("target_name") or target)
    y_da = open_field(y_item["filename"], y_name)
    ancestors.append(y_item["filename"])

    x = dataset_from_fields(x_fields)
    y = xr_dataset({y_name: y_da})
    return x, y, ancestors


def xr_dataset(fields):
    """Wrap a mapping of DataArrays as an xarray Dataset."""
    return xr.Dataset(fields)


def _plot_path(cfg, basename: str) -> Path:
    ext = cfg.get("output_file_type") or "png"
    if isinstance(ext, (list, tuple)):
        ext = ext[0]
    return Path(cfg["plot_dir"]) / str(ext) / f"{basename}.{ext}"


def _save_fig(basename: str, caption: str, cfg, ancestors, fig) -> Path | None:
    if fig is None:
        return None
    save_figure(basename, _provenance(cfg, ancestors, caption), cfg, figure=fig)
    return _plot_path(cfg, basename)


def _shap_xy(result):
    shap_res = result.artifacts.get("shap_result")
    if shap_res is None:
        return None
    x = result.split.test.X[shap_res.sample_positions]
    # ShapResult.values is a numpy array attribute, not pandas.Series.values.
    shap_values = np.asarray(getattr(shap_res, "values"), dtype=np.float64)
    return x, shap_values, list(shap_res.feature_names)


def _source_caption(base: str, source_label: str) -> str:
    return f"{base} ({source_label})" if source_label else base


def _try_fig(plots, key, basename, caption, cfg, ancestors, factory) -> None:
    try:
        path = _save_fig(basename, caption, cfg, ancestors, factory())
        if path:
            plots[key] = path
    except Exception:
        logger.exception("Could not plot %s", key)


def _write_prediction_outputs(result, cfg, ancestors, tag, source_label, plots, prov) -> None:
    if result.prediction_cube is None:
        return
    _save_nc(
        f"{tag}prediction",
        result.prediction_cube.to_dataset(name="prediction"),
        cfg,
        prov,
    )

    def _factory():
        y_cube = reconstruct_prediction(
            result.split.test.y,
            result.split.test,
            name="obs",
        )
        return plot_error_map(y_cube, result.prediction_cube)

    _try_fig(
        plots,
        "prediction_maps",
        f"{tag}prediction_maps",
        _source_caption("RF predicted versus true target (time mean)", source_label),
        cfg,
        ancestors,
        _factory,
    )


def _write_shap_outputs(result, cfg, ancestors, tag, source_label, plots) -> None:
    if result.shap_cube is not None:
        shap_cap = _source_caption(
            "SHAP restored to a (time, lat, lon) cube; time-mean of test-period cells",
            source_label,
        )
        shap_prov = _provenance(cfg, ancestors, shap_cap)
        _save_nc(f"{tag}shap", result.shap_cube, cfg, shap_prov)
        _try_fig(
            plots,
            "shap_maps",
            f"{tag}shap_maps",
            shap_cap,
            cfg,
            ancestors,
            lambda: plot_shap_maps(result.shap_cube),
        )

    shap_xy = _shap_xy(result)
    if shap_xy is None:
        return
    xs, shap_vals, names = shap_xy
    _try_fig(
        plots,
        "shap_summary",
        f"{tag}shap_summary",
        _source_caption("SHAP summary (colour is the feature value)", source_label),
        cfg,
        ancestors,
        lambda: plot_shap_summary(shap_vals, xs, names),
    )
    _try_fig(
        plots,
        "shap_dependence",
        f"{tag}shap_dependence",
        _source_caption(
            "SHAP dependence: x is the driver, y is SHAP of that driver",
            source_label,
        ),
        cfg,
        ancestors,
        lambda: plot_shap_dependence(xs, shap_vals, names),
    )


def _write_importance_outputs(result, cfg, ancestors, tag, source_label, plots, prov) -> None:
    if "rf_importance" in result.artifacts:
        _save_csv(f"{tag}rf_importance", result.artifacts["rf_importance"], cfg, prov)
        _try_fig(
            plots,
            "rf_importance",
            f"{tag}rf_importance",
            _source_caption("Tree feature importance", source_label),
            cfg,
            ancestors,
            lambda: plot_importance_bars(
                result.artifacts["rf_importance"],
                title="tree importance",
            ),
        )
    if "permutation_importance" in result.artifacts:
        frame = result.artifacts["permutation_importance"]
        _save_csv(f"{tag}permutation_importance", frame, cfg, prov)
        _try_fig(
            plots,
            "permutation_importance",
            f"{tag}permutation_importance",
            _source_caption("Permutation importance", source_label),
            cfg,
            ancestors,
            lambda: plot_importance_bars(frame, title="permutation importance"),
        )


def _write_effect_outputs(result, cfg, ancestors, tag, source_label, plots) -> None:
    if "ale" in result.artifacts:
        _try_fig(
            plots,
            "ale",
            f"{tag}ale",
            _source_caption("Accumulated Local Effects", source_label),
            cfg,
            ancestors,
            lambda: plot_ale(result.artifacts["ale"]),
        )
    if "ice_pdp" in result.artifacts:
        _try_fig(
            plots,
            "ice_pdp",
            f"{tag}ice_pdp",
            _source_caption("ICE (thin) and PDP (thick)", source_label),
            cfg,
            ancestors,
            lambda: plot_ice_pdp(result.artifacts["ice_pdp"]),
        )


def write_outputs(
    result,
    cfg,
    ancestors: list[str],
    prefix: str = "",
    source_label: str = "",
) -> dict:
    """Write cubes, tables and figures for one fitted mapping. Returns plot paths."""
    tag = f"{prefix}_" if prefix else ""
    plots: dict[str, Path] = {}
    caption = _source_caption(
        cfg.get("caption", "rf_shap diagnostic of collocated X versus Y"),
        source_label,
    )
    prov = _provenance(cfg, ancestors, caption)
    metrics = pd.Series(result.metrics).to_frame("value")
    _save_csv(f"{tag}metrics", metrics, cfg, prov)
    if prefix == "cmip":
        _save_csv("metrics", metrics, cfg, prov)

    _write_prediction_outputs(result, cfg, ancestors, tag, source_label, plots, prov)
    _write_shap_outputs(result, cfg, ancestors, tag, source_label, plots)
    _write_importance_outputs(result, cfg, ancestors, tag, source_label, plots, prov)
    _write_effect_outputs(result, cfg, ancestors, tag, source_label, plots)

    _try_fig(
        plots,
        "pred_vs_obs",
        f"{tag}pred_vs_obs",
        _source_caption("Predicted versus true target on the test sample", source_label),
        cfg,
        ancestors,
        lambda: plot_pred_scatter(
            result.split.test.y,
            result.prediction,
            result.metrics,
        ),
    )
    label = f" [{source_label}]" if source_label else ""
    logger.info("Metrics%s: %s", label, result.metrics)
    return plots


def fit_from_cfg(cfg):
    """Load X/Y the way this diagnostic does and fit one pipeline."""
    if cfg.get("use_synthetic"):
        x, y = make_synthetic_cubes(
            n_time=int(cfg.get("synthetic_n_time", 24)),
            n_lat=int(cfg.get("synthetic_n_lat", 8)),
            n_lon=int(cfg.get("synthetic_n_lon", 10)),
            seed=int(cfg.get("seed", 7)),
        )
        ancestors = []
    elif cfg.get("x_files") and cfg.get("y_file"):
        y_name = str(cfg.get("target_name") or cfg.get("target") or "y")
        x, y = load_xy_from_files(cfg["x_files"], cfg["y_file"], y_name=y_name)
        ancestors = [*list(cfg["x_files"].values()), cfg["y_file"]]
    else:
        x, y, ancestors = load_xy_from_esmvaltool(cfg)
    logger.info("X variables: %s", list(x.data_vars))
    logger.info("Y variables: %s", list(y.data_vars))
    result = run_pipeline_on_cubes(
        x,
        y,
        override=cfg.get("rf_shap") or {},
        regrid=cfg.get("regrid", "nearest"),
        join_time=cfg.get("join_time", "inner"),
        align_target=cfg.get("align_target", "y"),
    )
    return result, ancestors


_COMPARISON_PLOTS = {
    "ale_distance": {
        "value_col": "distance",
        "caption": "ALE curve distance of each driver: process model versus observations",
        "xlabel": "ALE RMSE (model vs observations)",
    },
    "dependence_distance": {
        "value_col": "distance",
        "caption": "SHAP dependence-curve distance of each driver: process model versus observations",
        "xlabel": "dependence-curve RMSE (model vs observations)",
    },
    "slope_difference": {
        "value_col": "difference",
        "caption": "Dependence-SHAP slope difference (process model minus observations)",
        "xlabel": "slope_model - slope_obs",
    },
    "interaction_distance": {
        "value_col": "distance",
        "caption": "SHAP interaction-curve distance: process model versus observations",
        "xlabel": "interaction-curve RMSE (model vs observations)",
    },
}


def _side_labels(cfg) -> tuple[str, str]:
    model = str(cfg.get("predictor_dataset") or cfg.get("target_dataset") or "CMIP")
    ref = cfg.get("reference") or {}
    obs_pred = ref.get("predictor_dataset") or "observations"
    obs_tgt = ref.get("target_dataset") or obs_pred
    obs = str(obs_pred) if obs_pred == obs_tgt else f"{obs_tgt} / {obs_pred}"
    return obs, model


def write_explanation_comparison(reference_result, other_result, cfg, ancestors: list[str]) -> dict:
    """Write curve distances and slope differences. Reference is the baseline."""
    reference_cfg = cfg.get("reference") or {}
    tables = compare_pipeline_results(
        reference_result,
        other_result,
        interaction=bool(reference_cfg.get("interaction", False)),
        n_bins=int(reference_cfg.get("n_bins", 20)),
    )
    plots: dict[str, Path] = {}
    for name, frame in tables.items():
        spec = _COMPARISON_PLOTS.get(
            name,
            {
                "value_col": "distance" if "distance" in frame.columns else frame.columns[-1],
                "caption": name.replace("_", " "),
                "xlabel": None,
            },
        )
        _save_csv(name, frame, cfg, _provenance(cfg, ancestors, spec["caption"]))
        _try_fig(
            plots,
            name,
            name,
            spec["caption"],
            cfg,
            ancestors,
            lambda frame=frame, spec=spec: plot_comparison_bars(
                frame,
                spec["value_col"],
                spec["caption"],
                xlabel=spec.get("xlabel"),
            ),
        )
    return plots


def write_paired_overlays(
    obs_result,
    cmip_result,
    cfg,
    ancestors: list[str],
    obs_label: str,
    cmip_label: str,
) -> dict:
    """Overlay and difference figures for the comparison column."""
    plots: dict[str, Path] = {}
    n_bins = int((cfg.get("reference") or {}).get("n_bins", 20))
    obs_xy, cmip_xy = _shap_xy(obs_result), _shap_xy(cmip_result)
    if obs_xy is not None and cmip_xy is not None:
        _try_fig(
            plots,
            "compare_shap_dependence",
            "compare_shap_dependence",
            "SHAP dependence curves: x is the driver, y is SHAP of that driver (observations vs CMIP)",
            cfg,
            ancestors,
            lambda: plot_shap_dependence_overlay(
                obs_xy[0],
                obs_xy[1],
                obs_xy[2],
                cmip_xy[0],
                cmip_xy[1],
                cmip_xy[2],
                obs_label=obs_label,
                cmip_label=cmip_label,
                n_bins=n_bins,
            ),
        )
    if obs_result.shap_cube is not None and cmip_result.shap_cube is not None:
        _try_fig(
            plots,
            "compare_shap_maps",
            "compare_shap_maps",
            "Time-mean SHAP (CMIP minus observations) on shared drivers",
            cfg,
            ancestors,
            lambda: plot_cube_difference(
                obs_result.shap_cube,
                cmip_result.shap_cube,
                title="SHAP time mean (CMIP - observations)",
            ),
        )
    if "ale" in obs_result.artifacts and "ale" in cmip_result.artifacts:
        _try_fig(
            plots,
            "compare_ale",
            "compare_ale",
            "ALE curves for observations and CMIP on shared drivers",
            cfg,
            ancestors,
            lambda: plot_ale_overlay(
                obs_result.artifacts["ale"],
                cmip_result.artifacts["ale"],
                obs_label=obs_label,
                cmip_label=cmip_label,
            ),
        )
    if "ice_pdp" in obs_result.artifacts and "ice_pdp" in cmip_result.artifacts:
        _try_fig(
            plots,
            "compare_pdp",
            "compare_pdp",
            "Partial dependence of observations and CMIP on shared drivers",
            cfg,
            ancestors,
            lambda: plot_pdp_overlay(
                obs_result.artifacts["ice_pdp"],
                cmip_result.artifacts["ice_pdp"],
                obs_label=obs_label,
                cmip_label=cmip_label,
            ),
        )
    if "rf_importance" in obs_result.artifacts and "rf_importance" in cmip_result.artifacts:
        _try_fig(
            plots,
            "compare_rf_importance",
            "compare_rf_importance",
            "Tree importance of observations versus CMIP",
            cfg,
            ancestors,
            lambda: plot_importance_overlay(
                obs_result.artifacts["rf_importance"],
                cmip_result.artifacts["rf_importance"],
                obs_label=obs_label,
                cmip_label=cmip_label,
            ),
        )
    if (
        "permutation_importance" in obs_result.artifacts
        and "permutation_importance" in cmip_result.artifacts
    ):
        _try_fig(
            plots,
            "compare_permutation_importance",
            "compare_permutation_importance",
            "Permutation importance of observations versus CMIP",
            cfg,
            ancestors,
            lambda: plot_importance_overlay(
                obs_result.artifacts["permutation_importance"],
                cmip_result.artifacts["permutation_importance"],
                obs_label=obs_label,
                cmip_label=cmip_label,
            ),
        )
    _try_fig(
        plots,
        "compare_skill",
        "compare_skill",
        "RF mapping skill on the test sample: observations versus CMIP",
        cfg,
        ancestors,
        lambda: plot_metrics_compare(
            obs_result.metrics,
            cmip_result.metrics,
            obs_label=obs_label,
            cmip_label=cmip_label,
        ),
    )
    return plots


def _write_gallery(cfg, obs_label, cmip_label, obs_plots, cmip_plots, compare_plots) -> None:
    rows = [
        {
            "title": "Predicted vs true target",
            "caption": "Scatter of ExtraTrees predictions on the held-out times.",
            "obs": obs_plots.get("pred_vs_obs"),
            "cmip": cmip_plots.get("pred_vs_obs"),
            "compare": compare_plots.get("compare_skill"),
        },
        {
            "title": "Prediction maps",
            "caption": "Time-mean true field, RF prediction, and prediction error.",
            "obs": obs_plots.get("prediction_maps"),
            "cmip": cmip_plots.get("prediction_maps"),
            "compare": compare_plots.get("compare_skill"),
        },
        {
            "title": "SHAP restored to cube",
            "caption": "Time-mean SHAP on test-period cells.",
            "obs": obs_plots.get("shap_maps"),
            "cmip": cmip_plots.get("shap_maps"),
            "compare": compare_plots.get("compare_shap_maps"),
        },
        {
            "title": "SHAP summary",
            "caption": "Each point is a sample; colour is the driver value.",
            "obs": obs_plots.get("shap_summary"),
            "cmip": cmip_plots.get("shap_summary"),
            "compare": compare_plots.get("dependence_distance"),
        },
        {
            "title": "SHAP dependence",
            "caption": (
                "X-axis is the driver; y-axis is SHAP of that driver. "
                "The comparison column overlays binned curves."
            ),
            "obs": obs_plots.get("shap_dependence"),
            "cmip": cmip_plots.get("shap_dependence"),
            "compare": compare_plots.get("compare_shap_dependence"),
        },
        {
            "title": "ALE",
            "caption": "Accumulated local effects of each driver.",
            "obs": obs_plots.get("ale"),
            "cmip": cmip_plots.get("ale"),
            "compare": compare_plots.get("compare_ale") or compare_plots.get("ale_distance"),
        },
        {
            "title": "ICE / PDP",
            "caption": "Individual conditional expectation (thin) and partial dependence (thick).",
            "obs": obs_plots.get("ice_pdp"),
            "cmip": cmip_plots.get("ice_pdp"),
            "compare": compare_plots.get("compare_pdp"),
        },
        {
            "title": "Tree importance",
            "obs": obs_plots.get("rf_importance"),
            "cmip": cmip_plots.get("rf_importance"),
            "compare": compare_plots.get("compare_rf_importance"),
        },
        {
            "title": "Permutation importance",
            "obs": obs_plots.get("permutation_importance"),
            "cmip": cmip_plots.get("permutation_importance"),
            "compare": compare_plots.get("compare_permutation_importance"),
        },
        {
            "title": "ALE curve distance",
            "caption": "RMSE between observational and CMIP ALE curves.",
            "compare": compare_plots.get("ale_distance"),
        },
        {
            "title": "SHAP dependence-curve distance",
            "caption": "RMSE between observational and CMIP SHAP dependence curves.",
            "compare": compare_plots.get("dependence_distance"),
        },
        {
            "title": "Dependence slope difference",
            "caption": "Slope of SHAP versus the driver (CMIP minus observations).",
            "compare": compare_plots.get("slope_difference"),
        },
    ]
    plot_dir = Path(cfg["plot_dir"])
    intro = (
        "Left: observational mapping (reference). Middle: CMIP mapping. "
        "Right: overlay or distance between the two."
    )
    gallery = write_three_column_html(
        plot_dir / "gallery.html",
        title=cfg.get("caption") or "rf_shap observations versus CMIP",
        obs_label=obs_label,
        cmip_label=cmip_label,
        rows=rows,
        intro=intro,
    )
    root = plot_dir.resolve().parents[2]
    extra = write_three_column_html(
        root / "rf_shap_gallery.html",
        title=cfg.get("caption") or "rf_shap observations versus CMIP",
        obs_label=obs_label,
        cmip_label=cmip_label,
        rows=rows,
        intro=(
            f"{intro} ESMValTool's index.html still lists every figure; "
            "this page is the three-column layout."
        ),
    )
    logger.info("Wrote three-column gallery %s and %s", gallery, extra)


def main(cfg) -> None:
    """Run the rf_shap diagnostic for one recipe task."""
    logger.info("Starting rf_shap ESMValTool diagnostic")
    result, ancestors = fit_from_cfg(cfg)
    reference = cfg.get("reference")
    if not reference:
        write_outputs(result, cfg, ancestors)
        return
    ref_cfg = dict(cfg)
    ref_cfg.pop("reference", None)
    ref_cfg.update(reference)
    logger.info("Fitting the reference (observational) explanation")
    reference_result, reference_ancestors = fit_from_cfg(ref_cfg)
    obs_label, cmip_label = _side_labels(cfg)
    all_anc = [*list(ancestors), *list(reference_ancestors)]
    obs_plots = write_outputs(
        reference_result,
        cfg,
        all_anc,
        prefix="obs",
        source_label=obs_label,
    )
    cmip_plots = write_outputs(
        result,
        cfg,
        all_anc,
        prefix="cmip",
        source_label=cmip_label,
    )
    compare_plots = write_explanation_comparison(reference_result, result, cfg, all_anc)
    compare_plots.update(
        write_paired_overlays(
            reference_result,
            result,
            cfg,
            all_anc,
            obs_label,
            cmip_label,
        ),
    )
    _write_gallery(cfg, obs_label, cmip_label, obs_plots, cmip_plots, compare_plots)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
