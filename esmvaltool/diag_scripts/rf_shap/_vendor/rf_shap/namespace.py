"""Expose rf_shap under another package name (e.g. esmvaltool.rf_shap)."""

from __future__ import annotations

import pkgutil
import sys
from types import ModuleType


def expose_package(dest_name: str) -> None:
    """Register ``dest_name.*`` as aliases of ``rf_shap.*`` (same module objects).

    Existing modules under ``dest_name`` (for example ``diagnostic``) are kept.
    """
    import rf_shap

    for info in pkgutil.walk_packages(rf_shap.__path__, prefix="rf_shap."):
        alias = dest_name + info.name[len("rf_shap") :]
        if alias in sys.modules:
            continue
        try:
            __import__(info.name)
        except Exception:
            continue
        sys.modules[alias] = sys.modules[info.name]


def install_namespace(dest: ModuleType) -> None:
    """Fill an ESMValTool-facing package with the rf_shap public API."""
    from rf_shap import (
        DiagnosisPipeline,
        PipelineConfig,
        fast_forest,
        fast_tree_shap,
        linear_tree_shap,
        load_config,
        regrid_curvilinear_nearest_to_1deg,
        regrid_file_to_1deg,
        regrid_to_1deg,
        run_pipeline,
    )
    from rf_shap.explain import (
        ale_curve_distance,
        compare_explanations,
        curve_fidelity,
        dependence_curve_distance,
        dependence_fidelity,
        dependence_slope_difference,
        dependence_slopes,
        evaluate_global_quality,
        evaluate_shap_quality,
        explanation_locality,
        global_method_agreement,
        interaction_curve_distance,
        interaction_shap,
        interaction_strength,
        seed_stability,
        run_attribution_verification,
    )
    from rf_shap import cube, explain, features, models, pipeline, predict
    from rf_shap import config as config_mod
    from rf_shap import plot, synthetic

    expose_package(dest.__name__)
    dest.DiagnosisPipeline = DiagnosisPipeline
    dest.PipelineConfig = PipelineConfig
    dest.load_config = load_config
    dest.run_pipeline = run_pipeline
    dest.fast_tree_shap = fast_tree_shap
    dest.linear_tree_shap = linear_tree_shap
    dest.fast_forest = fast_forest
    dest.regrid_to_1deg = regrid_to_1deg
    dest.regrid_curvilinear_nearest_to_1deg = (
        regrid_curvilinear_nearest_to_1deg
    )
    dest.regrid_file_to_1deg = regrid_file_to_1deg
    dest.compare_explanations = compare_explanations
    dest.ale_curve_distance = ale_curve_distance
    dest.curve_fidelity = curve_fidelity
    dest.dependence_curve_distance = dependence_curve_distance
    dest.dependence_fidelity = dependence_fidelity
    dest.dependence_slopes = dependence_slopes
    dest.dependence_slope_difference = dependence_slope_difference
    dest.interaction_curve_distance = interaction_curve_distance
    dest.interaction_shap = interaction_shap
    dest.evaluate_shap_quality = evaluate_shap_quality
    dest.evaluate_global_quality = evaluate_global_quality
    dest.explanation_locality = explanation_locality
    dest.global_method_agreement = global_method_agreement
    dest.interaction_strength = interaction_strength
    dest.seed_stability = seed_stability
    dest.run_attribution_verification = run_attribution_verification
    dest.cube = cube
    dest.models = models
    dest.explain = explain
    dest.features = features
    dest.predict = predict
    dest.pipeline = pipeline
    dest.config = config_mod
    dest.plot = plot
    dest.synthetic = synthetic
    dest.__all__ = [
        "DiagnosisPipeline",
        "PipelineConfig",
        "load_config",
        "run_pipeline",
        "fast_tree_shap",
        "linear_tree_shap",
        "fast_forest",
        "regrid_to_1deg",
        "regrid_curvilinear_nearest_to_1deg",
        "regrid_file_to_1deg",
        "compare_explanations",
        "ale_curve_distance",
        "curve_fidelity",
        "dependence_curve_distance",
        "dependence_fidelity",
        "dependence_slopes",
        "dependence_slope_difference",
        "interaction_curve_distance",
        "interaction_shap",
        "evaluate_shap_quality",
        "evaluate_global_quality",
        "explanation_locality",
        "global_method_agreement",
        "interaction_strength",
        "seed_stability",
        "run_attribution_verification",
        "cube",
        "models",
        "explain",
        "features",
        "predict",
        "pipeline",
        "config",
        "plot",
        "synthetic",
    ]
