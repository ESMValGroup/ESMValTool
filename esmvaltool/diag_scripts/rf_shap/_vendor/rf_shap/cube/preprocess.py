from __future__ import annotations

from typing import Literal

import numpy as np
import xarray as xr
from sklearn.linear_model import LinearRegression

Climatology = Literal["month", "dayofyear"]


def _group_index(da: xr.DataArray, how: Climatology) -> xr.DataArray:
    if how == "month":
        return da.time.dt.month
    return da.time.dt.dayofyear


def detrend_along_time(da: xr.DataArray) -> xr.DataArray:
    time = np.arange(da.sizes["time"], dtype=np.float64)
    values = da.transpose("time", "lat", "lon").values
    n_t, n_y, n_x = values.shape
    flat = values.reshape(n_t, -1)
    out = np.full_like(flat, np.nan, dtype=np.float64)
    t = time.reshape(-1, 1)
    model = LinearRegression()
    for i in range(flat.shape[1]):
        col = flat[:, i]
        mask = np.isfinite(col)
        if mask.sum() < 3:
            continue
        model.fit(t[mask], col[mask])
        out[mask, i] = col[mask] - model.predict(t[mask])
    return xr.DataArray(out.reshape(n_t, n_y, n_x), coords=da.coords, dims=("time", "lat", "lon"), name=da.name)


def deseasonalize(da: xr.DataArray, how: Climatology = "month") -> xr.DataArray:
    try:
        idx = _group_index(da, how)
        clim = da.groupby(idx).mean("time")
        return da.groupby(idx) - clim
    except (TypeError, AttributeError, ValueError):
        n = 12 if how == "month" else min(365, da.sizes["time"])
        idx = np.arange(da.sizes["time"]) % n
        grouped = da.assign_coords(_season=("time", idx)).groupby("_season")
        return grouped - grouped.mean("time")


def anomaly(da: xr.DataArray) -> xr.DataArray:
    return da - da.mean("time")


def standardize(da: xr.DataArray) -> xr.DataArray:
    std = da.std("time")
    std = xr.where(std == 0, np.nan, std)
    return (da - da.mean("time")) / std


def preprocess_dataarray(
    da: xr.DataArray,
    *,
    detrend: bool = False,
    deseasonalize_flag: bool = False,
    anomaly_flag: bool = False,
    standardize_flag: bool = False,
    climatology: Climatology = "month",
) -> xr.DataArray:
    out = da
    if detrend:
        out = detrend_along_time(out)
    if deseasonalize_flag:
        out = deseasonalize(out, climatology)
    if anomaly_flag:
        out = anomaly(out)
    if standardize_flag:
        out = standardize(out)
    out.name = da.name
    return out


def preprocess_dataset(ds: xr.Dataset, options: dict | None = None) -> xr.Dataset:
    options = options or {}
    processed = {}
    for name, da in ds.data_vars.items():
        processed[name] = preprocess_dataarray(
            da,
            detrend=bool(options.get("detrend", False)),
            deseasonalize_flag=bool(options.get("deseasonalize", False)),
            anomaly_flag=bool(options.get("anomaly", False)),
            standardize_flag=bool(options.get("standardize", False)),
            climatology=options.get("climatology", "month"),
        )
    return xr.Dataset(processed)
