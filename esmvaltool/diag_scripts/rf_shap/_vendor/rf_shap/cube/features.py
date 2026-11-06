from __future__ import annotations

from typing import Sequence

import numpy as np
import xarray as xr
from sklearn.decomposition import PCA


def add_lags(ds: xr.Dataset, lags: Sequence[int]) -> xr.Dataset:
    out = ds.copy()
    for lag in lags:
        if int(lag) <= 0:
            raise ValueError("lags must be positive integers")
        shifted = ds.shift(time=int(lag))
        rename = {name: f"{name}_lag{lag}" for name in ds.data_vars}
        out = out.merge(shifted.rename(rename))
    return out.isel(time=slice(max(int(l) for l in lags), None))


def pca_transform_table(x: np.ndarray, n_components=0.95, whiten: bool = False, random_state: int = 42):
    model = PCA(n_components=n_components, whiten=whiten, random_state=random_state)
    transformed = model.fit_transform(x)
    names = [f"pc{i+1}" for i in range(transformed.shape[1])]
    return transformed, names, model


def eof_compress(ds: xr.Dataset, n_components: int = 5) -> xr.Dataset:
    """Optional spatial EOF compression: each variable becomes leading PC fields reconstructed.

    This keeps the cube shape so X/Y still share time/lat/lon.
    """
    out = {}
    for name, da in ds.data_vars.items():
        values = da.transpose("time", "lat", "lon").values
        n_t, n_y, n_x = values.shape
        flat = values.reshape(n_t, -1)
        col_mean = np.nanmean(flat, axis=0)
        filled = np.where(np.isfinite(flat), flat, col_mean)
        k = min(n_components, n_t, filled.shape[1])
        model = PCA(n_components=k)
        scores = model.fit_transform(filled)
        recon = model.inverse_transform(scores).reshape(n_t, n_y, n_x)
        recon[~np.isfinite(values)] = np.nan
        out[name] = xr.DataArray(recon, coords=da.coords, dims=da.dims, name=name)
    return xr.Dataset(out)


def build_x_features(
    x: xr.Dataset,
    mode: str = "raw",
    lags: Sequence[int] | None = None,
    pca_n_components=0.95,
    pca_whiten: bool = False,
) -> tuple[xr.Dataset, dict]:
    """Build X options while preserving a (time, lat, lon) cube.

    Modes:
      raw: original predictors
      lag: original + time lags
      pca: spatial EOF reconstruction of each variable (cube-preserving)
      lag_pca: lag then EOF reconstruction
    Feature-space PCA on the flattened table is applied later in the pipeline.
    """
    meta = {"mode": mode, "lags": list(lags or [])}
    cube = x
    if mode in {"lag", "lag_pca"}:
        cube = add_lags(cube, lags or [1])
    if mode in {"pca", "lag_pca"}:
        n = pca_n_components if isinstance(pca_n_components, int) else 8
        cube = eof_compress(cube, n_components=int(n))
        meta["eof_n_components"] = int(n)
    meta["feature_pca"] = pca_n_components if mode in {"pca", "lag_pca"} else None
    meta["pca_whiten"] = pca_whiten
    return cube, meta
