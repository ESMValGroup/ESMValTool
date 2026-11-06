"""Compare two fitted explanations inside the ESMValTool diagnostic.

The reference result is the baseline (typically observations). The other
result is the explanation being judged (typically a process model).
"""

from __future__ import annotations

import numpy as np

from rf_shap.explain.evaluate import compare_explanations, interaction_shap


def _xy_shap(result):
    shap_res = result.artifacts.get("shap_result")
    if shap_res is None:
        return None
    positions = shap_res.sample_positions
    x = np.asarray(result.split.test.X[positions], dtype=np.float64)
    values = np.asarray(shap_res.values, dtype=np.float64)
    names = [str(name) for name in shap_res.feature_names]
    return x, values, names


def _shared_features(names_ref, names_other):
    common = [name for name in names_ref if name in names_other]
    if not common:
        raise ValueError(
            "The two explanations share no features: "
            f"{list(names_ref)} versus {list(names_other)}"
        )
    return common


def _columns(matrix, names, keep):
    index = [names.index(name) for name in keep]
    return np.asarray(matrix, dtype=np.float64)[:, index]


def compare_pipeline_results(reference, other, *, interaction: bool = False, n_bins: int = 20):
    """Curve distances and dependence-slope differences between two pipeline fits.

    ``reference`` is the baseline explanation. ``other`` is compared against it.
    Interactive SHAP is computed only when ``interaction`` is true.
    Features that appear in only one fit (for example after collinearity
    screening) are dropped; comparison uses the shared names.
    """
    kwargs = {}
    ale_ref = reference.artifacts.get("ale")
    ale_other = other.artifacts.get("ale")

    ref_view = _xy_shap(reference)
    other_view = _xy_shap(other)
    keep = None
    if ref_view is not None and other_view is not None:
        x_ref, shap_ref, names_ref = ref_view
        x_other, shap_other, names_other = other_view
        keep = _shared_features(names_ref, names_other)
        kwargs["x_ref"] = _columns(x_ref, names_ref, keep)
        kwargs["shap_ref"] = _columns(shap_ref, names_ref, keep)
        kwargs["x_other"] = _columns(x_other, names_other, keep)
        kwargs["shap_other"] = _columns(shap_other, names_other, keep)
        kwargs["feature_names"] = keep
        if interaction:
            kwargs["interaction_ref"] = interaction_shap(reference.model, kwargs["x_ref"])
            kwargs["interaction_other"] = interaction_shap(other.model, kwargs["x_other"])
    if ale_ref is not None and ale_other is not None:
        if keep is None:
            keep = _shared_features(list(ale_ref), list(ale_other))
        kwargs["ale_ref"] = {name: ale_ref[name] for name in keep if name in ale_ref}
        kwargs["ale_other"] = {name: ale_other[name] for name in keep if name in ale_other}
    return compare_explanations(n_bins=n_bins, **kwargs)
