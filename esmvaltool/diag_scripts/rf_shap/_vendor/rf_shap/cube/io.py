from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
import xarray as xr

from rf_shap.utils import check_names

STANDARD_DIMS = ("time", "lat", "lon")


def _rename_dims(ds: xr.Dataset, time_name: str, lat_name: str, lon_name: str) -> xr.Dataset:
    mapping = {}
    for src, dst in [(time_name, "time"), (lat_name, "lat"), (lon_name, "lon")]:
        if src in ds.dims and src != dst:
            mapping[src] = dst
        elif src in ds.coords and src != dst and dst not in ds.dims:
            mapping[src] = dst
    if mapping:
        ds = ds.rename(mapping)
    missing = [name for name in STANDARD_DIMS if name not in ds.dims and name not in ds.coords]
    if missing:
        raise ValueError(f"Cube must contain time/lat/lon. Missing: {missing}. Dims={list(ds.dims)}")
    return ds


def _ensure_cube(ds: xr.Dataset, variables: Sequence[str] | None = None) -> xr.Dataset:
    keep = list(variables) if variables else [v for v in ds.data_vars]
    check_names(keep, ds.data_vars, "cube")
    ds = ds[keep]
    for name in keep:
        da = ds[name]
        if set(STANDARD_DIMS) - set(da.dims):
            raise ValueError(f"Variable {name} must have dims time, lat, lon; got {da.dims}")
        ds[name] = da.transpose("time", "lat", "lon")
    return ds


def dataset_from_arrays(
    arrays: Mapping[str, np.ndarray],
    time,
    lat,
    lon,
) -> xr.Dataset:
    data_vars = {}
    for name, values in arrays.items():
        arr = np.asarray(values)
        if arr.ndim != 3:
            raise ValueError(f"{name} must be (time, lat, lon), got {arr.shape}")
        data_vars[name] = (("time", "lat", "lon"), arr)
    return xr.Dataset(
        data_vars,
        coords={"time": np.asarray(time), "lat": np.asarray(lat), "lon": np.asarray(lon)},
    )


def table_to_dataset(
    frame: pd.DataFrame,
    variables: Sequence[str],
    time_name: str = "time",
    lat_name: str = "lat",
    lon_name: str = "lon",
) -> xr.Dataset:
    needed = [time_name, lat_name, lon_name, *variables]
    check_names(needed, frame.columns, "table")
    work = frame[needed].copy()
    work = work.rename(columns={time_name: "time", lat_name: "lat", lon_name: "lon"})
    work = work.set_index(["time", "lat", "lon"]).sort_index()
    if work.index.duplicated().any():
        work = work.groupby(level=[0, 1, 2]).mean()
    ds = work.to_xarray()
    return _ensure_cube(ds, variables)


def load_netcdf(
    path: str | Path,
    variables: Sequence[str] | None = None,
    time_name: str = "time",
    lat_name: str = "lat",
    lon_name: str = "lon",
) -> xr.Dataset:
    ds = xr.open_dataset(path)
    ds = _rename_dims(ds, time_name, lat_name, lon_name)
    return _ensure_cube(ds, variables)


def load_npy(
    path: str | Path,
    variables: Sequence[str] | None = None,
    time_name: str = "time",
    lat_name: str = "lat",
    lon_name: str = "lon",
) -> xr.Dataset:
    path = Path(path)
    if path.suffix == ".npz":
        payload = np.load(path, allow_pickle=True)
        keys = list(payload.files)
        time = payload[time_name]
        lat = payload[lat_name]
        lon = payload[lon_name]
        var_names = variables or [k for k in keys if k not in {time_name, lat_name, lon_name}]
        arrays = {name: payload[name] for name in var_names}
        return _ensure_cube(dataset_from_arrays(arrays, time, lat, lon), var_names)

    if path.is_dir():
        time = np.load(path / f"{time_name}.npy")
        lat = np.load(path / f"{lat_name}.npy")
        lon = np.load(path / f"{lon_name}.npy")
        files = sorted(p for p in path.glob("*.npy") if p.stem not in {time_name, lat_name, lon_name})
        var_names = list(variables) if variables else [p.stem for p in files]
        arrays = {name: np.load(path / f"{name}.npy") for name in var_names}
        return _ensure_cube(dataset_from_arrays(arrays, time, lat, lon), var_names)

    arr = np.load(path)
    if arr.ndim != 3:
        raise ValueError("Single npy file must be shaped (time, lat, lon)")
    name = variables[0] if variables else path.stem
    n_t, n_y, n_x = arr.shape
    ds = dataset_from_arrays(
        {name: arr},
        np.arange(n_t),
        np.linspace(-90, 90, n_y),
        np.linspace(0, 360, n_x, endpoint=False),
    )
    return _ensure_cube(ds, [name])


def load_table(
    path: str | Path,
    variables: Sequence[str] | None = None,
    time_name: str = "time",
    lat_name: str = "lat",
    lon_name: str = "lon",
) -> xr.Dataset:
    path = Path(path)
    if path.suffix.lower() in {".parquet", ".pq"}:
        frame = pd.read_parquet(path)
    else:
        frame = pd.read_csv(path)
    if time_name in frame.columns:
        frame[time_name] = pd.to_datetime(frame[time_name], errors="ignore")
    vars_ = list(variables) if variables else [c for c in frame.columns if c not in {time_name, lat_name, lon_name}]
    return table_to_dataset(frame, vars_, time_name, lat_name, lon_name)


def load_cube(
    path: str | Path,
    variables: Sequence[str] | None = None,
    time_name: str = "time",
    lat_name: str = "lat",
    lon_name: str = "lon",
) -> xr.Dataset:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {".nc", ".nc4", ".cdf"} or suffix == "":
        if path.is_dir():
            return load_npy(path, variables, time_name, lat_name, lon_name)
        return load_netcdf(path, variables, time_name, lat_name, lon_name)
    if suffix in {".npy", ".npz"}:
        return load_npy(path, variables, time_name, lat_name, lon_name)
    if suffix in {".csv", ".parquet", ".pq", ".txt"}:
        return load_table(path, variables, time_name, lat_name, lon_name)
    raise ValueError(f"Unsupported data path: {path}")
