from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import pandas as pd
import xarray as xr

from rf_shap.cube.tabular import SampleTable, values_to_cubes


@dataclass
class ShapResult:
    values: np.ndarray
    expected_value: float | np.ndarray
    method: str
    index: pd.MultiIndex
    feature_names: list[str]
    sample_positions: np.ndarray
    algorithm: str = ""

    def to_frame(self) -> pd.DataFrame:
        frame = pd.DataFrame(self.values, columns=self.feature_names, index=self.index[self.sample_positions])
        return frame.reset_index()


def _try_import_shap():
    try:
        import shap
        return shap
    except ImportError as exc:
        raise ImportError("Install shap to compute SHAP values: pip install shap") from exc


def _prepare_numpy_for_fasttreeshap():
    """FastTreeSHAP 0.1.6 still calls ``np.obj2sctype``, removed in NumPy 2."""
    import numpy as np

    if hasattr(np, "obj2sctype"):
        return

    def obj2sctype(rep, default=None):
        try:
            if isinstance(rep, type) and issubclass(rep, np.generic):
                return rep
        except TypeError:
            pass
        try:
            return np.dtype(rep).type
        except (TypeError, ValueError):
            return default

    np.obj2sctype = obj2sctype


def _try_import_fasttreeshap():
    _prepare_numpy_for_fasttreeshap()
    try:
        import fasttreeshap
        return fasttreeshap
    except ImportError:
        return None


def _is_tree_model(model) -> bool:
    if getattr(model, "supports_tree_shap", None) is False:
        return False
    name = type(model).__name__.lower()
    tree_bits = (
        "forest",
        "tree",
        "xgb",
        "lgbm",
        "lightgbm",
        "histgradient",
        "extratrees",
        "gbm",
        "catboost",
    )
    return any(bit in name for bit in tree_bits)


def _background(x: np.ndarray, size: int, seed: int) -> np.ndarray:
    shap = _try_import_shap()
    n = min(size, len(x))
    if n <= 0:
        return x
    try:
        return shap.kmeans(x, n).data
    except Exception:
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(x), size=n, replace=False)
        return x[idx]


def _subset(x: np.ndarray, max_samples: int | None, seed: int) -> np.ndarray:
    if max_samples is None or max_samples >= len(x):
        return np.arange(len(x))
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(len(x), size=int(max_samples), replace=False))


def fast_tree_shap(
    model,
    x: np.ndarray,
    *,
    algorithm: str = "v2",
    n_jobs: int = -1,
    background: np.ndarray | None = None,
    feature_perturbation: str | None = None,
    memory_tolerance: float = -1,
    feature_names: Sequence[str] | None = None,
    check_additivity: bool = False,
) -> tuple[np.ndarray, float | np.ndarray]:
    """Fast SHAP for a random forest or other tree model.

    Calls LinkedIn FastTreeSHAP (Yang, 2021, https://github.com/linkedin/fasttreeshap)
    instead of copying that implementation. On scikit-learn random forests this is
    the fast option: stock ``shap.TreeExplainer`` explains sklearn trees on one
    core, while FastTreeSHAP v1/v2 runs in parallel. Published single-core gains
    are about 1.5x (v1, same memory) and 2.5x (v2, more memory); with all cores
    the sklearn-RF speedups in that paper are roughly 7–15x.

    ``algorithm`` is ``v0`` (original TreeSHAP), ``v1``, ``v2``, or ``auto``.
    The default perturbation is path-dependent, which needs no background sample
    and does not slow down as the training set grows. Pass a small ``background``
    and ``feature_perturbation="interventional"`` when features are strongly
    dependent and that causal variant is required.

    Install the optional backend with ``pip install fasttreeshap``. On macOS,
    ``brew install libomp`` is required if the build cannot find OpenMP.
    """
    fast = _try_import_fasttreeshap()
    if fast is None:
        raise ImportError(
            "fast_tree_shap needs the LinkedIn FastTreeSHAP package: pip install fasttreeshap"
        )
    algo = algorithm if algorithm in {"v0", "v1", "v2", "auto"} else "v2"
    if feature_perturbation is None:
        feature_perturbation = "interventional" if background is not None else "tree_path_dependent"
    data = None
    if feature_perturbation == "interventional":
        if background is None:
            raise ValueError("interventional FastTreeSHAP requires a background sample")
        data = np.asarray(background, dtype=np.float64)
    explainer = fast.TreeExplainer(
        model,
        data=data,
        feature_perturbation=feature_perturbation,
        algorithm=algo,
        n_jobs=-1 if n_jobs is None else n_jobs,
        memory_tolerance=memory_tolerance,
        feature_names=None if feature_names is None else list(feature_names),
    )
    values = explainer.shap_values(np.asarray(x, dtype=np.float64), check_additivity=check_additivity)
    if isinstance(values, list):
        values = values[-1]
    return _normalize_values(values), explainer.expected_value


def _sklearn_tree_estimators(model) -> list:
    if hasattr(model, "estimators_"):
        return [estimator for estimator in np.ravel(model.estimators_) if hasattr(estimator, "tree_")]
    if hasattr(model, "tree_"):
        return [model]
    return []


def _root_prediction(estimator) -> float:
    return float(np.asarray(estimator.tree_.value[0]).reshape(-1)[0])


def linear_tree_shap(
    model,
    x: np.ndarray,
    *,
    version: str = "v2",
) -> tuple[np.ndarray, float]:
    """Exact linear-time SHAP for one scikit-learn tree or a forest of them.

    Calls yupbank/linear_tree_shap (https://github.com/yupbank/linear_tree_shap)
    instead of copying it. For one tree of depth D and L leaves the cost is
    O(LD) and does not store an extra polynomial per sample. ``version="v2"``
    uses Chebyshev interpolation so trees deeper than about 12 do not overflow
    the way the original polynomial form can. A random forest is the mean of
    its trees, so the forest SHAP values are the mean of the per-tree values.

    This does not replace ``fast_tree_shap``. XGBoost and LightGBM are not
    scikit-learn trees and are left to those other explainers.
    """
    try:
        import linear_tree_shap
    except ImportError as exc:
        raise ImportError(
            "linear_tree_shap needs https://github.com/yupbank/linear_tree_shap "
            "installed with: pip install --no-deps --no-build-isolation "
            "git+https://github.com/yupbank/linear_tree_shap.git"
        ) from exc
    estimators = _sklearn_tree_estimators(model)
    if not estimators:
        raise TypeError("linear_tree_shap only supports scikit-learn trees and forests")
    features = np.asarray(x, dtype=np.float64)
    if features.ndim == 1:
        features = features.reshape(1, -1)
    use_v2 = str(version).lower() != "v1"
    parts = []
    for estimator in estimators:
        tree = estimator.tree_
        weighted = getattr(tree, "weighted_n_node_samples", None)
        saved_counts = None
        if weighted is not None and not np.array_equal(weighted, tree.n_node_samples):
            # Bootstrap forests keep duplicate-sample counts in weighted_n_node_samples.
            # linear_tree_shap.copy_tree reads n_node_samples, which drops that multiplicity.
            saved_counts = np.array(tree.n_node_samples, copy=True)
            tree.n_node_samples[:] = np.rint(weighted).astype(saved_counts.dtype, copy=False)
        try:
            explainer = linear_tree_shap.TreeExplainer(estimator)
            if use_v2:
                values = explainer.shap_values_v2(features)
            else:
                values = explainer.shap_values(features)
        except ImportError:
            rows = [
                explainer.py_shap_values_v2(row) if use_v2 else explainer.py_shap_values(row)
                for row in features
            ]
            values = np.vstack(rows)
        finally:
            if saved_counts is not None:
                tree.n_node_samples[:] = saved_counts
        parts.append(np.asarray(values, dtype=np.float64))
    values = np.mean(parts, axis=0) if len(parts) > 1 else parts[0]
    expected = float(np.mean([_root_prediction(estimator) for estimator in estimators]))
    return _normalize_values(values), expected


def _tree_shap(model, x_explain, background, n_jobs: int):
    shap = _try_import_shap()
    kwargs = {}
    if background is None:
        explainer = shap.TreeExplainer(model)
    else:
        explainer = shap.TreeExplainer(model, data=background, feature_perturbation="interventional")
    if hasattr(explainer, "shap_values"):
        try:
            values = explainer.shap_values(x_explain, check_additivity=False)
        except TypeError:
            values = explainer.shap_values(x_explain)
    else:
        values = explainer(x_explain).values
    return np.asarray(values), explainer.expected_value, "treeshap"


def _kernel_shap(model, x_explain, background):
    shap = _try_import_shap()
    explainer = shap.KernelExplainer(model.predict, background)
    values = explainer.shap_values(x_explain, nsamples="auto")
    return np.asarray(values), explainer.expected_value, "kernel"


def _permutation_shap(model, x_explain, background):
    shap = _try_import_shap()
    explainer = shap.PermutationExplainer(model.predict, background)
    explanation = explainer(x_explain)
    return np.asarray(explanation.values), explanation.base_values, "permutation"


def _sampling_shap(model, x_explain, background):
    shap = _try_import_shap()
    explainer = shap.SamplingExplainer(model.predict, background)
    values = explainer.shap_values(x_explain)
    return np.asarray(values), explainer.expected_value, "sampling"


def _normalize_values(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values)
    if arr.ndim == 3:
        # classification-style (n, features, outputs) or (outputs, n, features)
        if arr.shape[-1] == 1:
            arr = arr[..., 0]
        elif arr.shape[0] == 1:
            arr = arr[0]
        else:
            arr = arr[..., -1]
    return arr


def _row_sample(x: np.ndarray, size: int, seed: int) -> np.ndarray:
    x = np.asarray(x)
    n = min(int(size), len(x))
    if n >= len(x):
        return x
    rng = np.random.default_rng(seed)
    return x[rng.choice(len(x), size=n, replace=False)]


def compute_shap(
    model,
    table: SampleTable,
    method: str = "auto",
    algorithm: str = "v2",
    max_samples: int | None = 2000,
    background_size: int = 100,
    n_jobs: int = -1,
    seed: int = 42,
    background: np.ndarray | None = None,
    feature_perturbation: str = "tree_path_dependent",
) -> ShapResult:
    """Compute SHAP with several speed/quality trade-offs.

    Efficiency options:
      auto/fasttreeshap : ``fast_tree_shap``, the fast RF + SHAP option
      treeshap          : path-dependent TreeSHAP (no background)
      linear_treeshap   : exact O(LD) SHAP for sklearn trees, v2 by default
      kernel/permutation/sampling : model-agnostic, use a small background
      max_samples       : only explain a subset, then restore those cells to the cube

    ``auto`` uses FastTreeSHAP v2 on tree models when the package is installed,
    and otherwise path-dependent TreeSHAP. Explicit ``fasttreeshap`` does not
    fall back: install ``fasttreeshap`` first.
    """
    x = np.asarray(table.X, dtype=np.float64)
    positions = _subset(x, max_samples, seed)
    x_explain = x[positions]
    method = method.lower()
    perturbation = feature_perturbation.lower()
    values = expected = used = None
    used_algorithm = ""

    if method in {"auto", "fasttreeshap"} and _is_tree_model(model):
        fast_background = None
        if perturbation == "interventional":
            raw_bg = x if background is None else background
            fast_background = _row_sample(raw_bg, background_size, seed)
        try:
            values, expected = fast_tree_shap(
                model,
                x_explain,
                algorithm=algorithm,
                n_jobs=n_jobs,
                background=fast_background,
                feature_perturbation=perturbation,
                feature_names=list(table.feature_names),
            )
            used = "fasttreeshap"
            used_algorithm = algorithm if algorithm in {"v0", "v1", "v2", "auto"} else "v2"
        except ImportError:
            if method == "fasttreeshap":
                raise

    if values is None and method in {"auto", "treeshap"} and _is_tree_model(model):
        # path-dependent TreeSHAP does not grow with a background set
        values, expected, used = _tree_shap(model, x_explain, None, n_jobs)

    if values is None and method in {"linear_treeshap", "lineartreeshap"}:
        linear_version = algorithm if algorithm in {"v1", "v2"} else "v2"
        values, expected = linear_tree_shap(model, x_explain, version=linear_version)
        used = "linear_treeshap"
        used_algorithm = linear_version

    if values is None:
        bg = background if background is not None else _background(x, background_size, seed)
        if len(bg) > background_size:
            bg = _row_sample(bg, background_size, seed)
        if method in {"kernel"}:
            values, expected, used = _kernel_shap(model, x_explain, bg)
        elif method in {"sampling"}:
            values, expected, used = _sampling_shap(model, x_explain, bg)
        else:
            values, expected, used = _permutation_shap(model, x_explain, bg)

    values = _normalize_values(values)
    if values.shape[0] != len(positions):
        raise ValueError(f"SHAP values have unexpected shape {values.shape}")
    return ShapResult(
        values=values,
        expected_value=expected,
        method=used or method,
        index=table.index,
        feature_names=list(table.feature_names),
        sample_positions=positions,
        algorithm=used_algorithm,
    )


def shap_to_cube(result: ShapResult, table: SampleTable) -> xr.Dataset:
    template = table.extras.get("template")
    if template is None:
        raise ValueError("SampleTable is missing a cube template; cannot restore SHAP")
    index = table.index[result.sample_positions]
    ds = values_to_cubes(result.values, index, template, result.feature_names)
    ds.attrs["shap_method"] = result.method
    ds.attrs["shap_algorithm"] = result.algorithm
    ds.attrs["expected_value"] = str(result.expected_value)
    return ds
