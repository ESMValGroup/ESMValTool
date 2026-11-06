"""CMIP6-like cube diagnosis with tree models and SHAP."""

from rf_shap.config import PipelineConfig, load_config
from rf_shap.cube.regrid import (
    regrid_curvilinear_nearest_to_1deg,
    regrid_file_to_1deg,
    regrid_to_1deg,
)
from rf_shap.explain.shap_explain import fast_tree_shap, linear_tree_shap
from rf_shap.models.fast_forest import fast_forest
from rf_shap.pipeline import DiagnosisPipeline, run_pipeline

__all__ = [
    "PipelineConfig",
    "load_config",
    "DiagnosisPipeline",
    "run_pipeline",
    "fast_tree_shap",
    "linear_tree_shap",
    "fast_forest",
    "regrid_to_1deg",
    "regrid_curvilinear_nearest_to_1deg",
    "regrid_file_to_1deg",
]
__version__ = "0.1.0"
