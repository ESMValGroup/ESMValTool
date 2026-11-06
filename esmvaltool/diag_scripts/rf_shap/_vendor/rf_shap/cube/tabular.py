from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
import pandas as pd
import xarray as xr


@dataclass
class SampleTable:
    X: np.ndarray
    y: np.ndarray
    feature_names: list[str]
    index: pd.MultiIndex
    y_name: str = "y"
    extras: dict = field(default_factory=dict)

    def __post_init__(self):
        self.X = np.asarray(self.X, dtype=np.float64)
        self.y = np.asarray(self.y, dtype=np.float64).reshape(-1)
        if self.X.ndim != 2:
            raise ValueError("X must be 2D (n_samples, n_features)")
        if len(self.X) != len(self.y) or len(self.X) != len(self.index):
            raise ValueError("X, y and index must share the same sample count")

    @property
    def n_samples(self) -> int:
        return len(self.y)

    def subset(self, mask) -> "SampleTable":
        mask = np.asarray(mask, dtype=bool)
        extras = dict(self.extras)
        extras["parent_mask"] = mask
        return SampleTable(self.X[mask], self.y[mask], list(self.feature_names), self.index[mask], self.y_name, extras)

    def take(self, idx) -> "SampleTable":
        idx = np.asarray(idx)
        return SampleTable(self.X[idx], self.y[idx], list(self.feature_names), self.index[idx], self.y_name, dict(self.extras))

    def with_features(self, X: np.ndarray, names: Sequence[str]) -> "SampleTable":
        return SampleTable(X, self.y, list(names), self.index, self.y_name, dict(self.extras))

    def select_columns(self, names: Sequence[str]) -> "SampleTable":
        pos = [self.feature_names.index(n) for n in names]
        return self.with_features(self.X[:, pos], names)

    def to_frame(self) -> pd.DataFrame:
        frame = pd.DataFrame(self.X, columns=self.feature_names, index=self.index)
        frame[self.y_name] = self.y
        return frame.reset_index()


def cube_to_table(
    x: xr.Dataset,
    y: xr.DataArray | xr.Dataset,
    y_name: str | None = None,
) -> SampleTable:
    if isinstance(y, xr.Dataset):
        if y_name is None:
            if len(y.data_vars) != 1:
                raise ValueError("Y dataset must contain a single variable or y_name")
            y_name = list(y.data_vars)[0]
        y_da = y[y_name]
    else:
        y_da = y
        y_name = y_da.name or "y"

    x, y_da = xr.align(x, y_da, join="inner")
    x_da = x.to_array(dim="feature").transpose("time", "lat", "lon", "feature")
    y_da = y_da.transpose("time", "lat", "lon")
    X = np.asarray(x_da.values, dtype=np.float64).reshape(-1, x_da.sizes["feature"])
    y_vals = np.asarray(y_da.values, dtype=np.float64).reshape(-1)
    times, lats, lons = np.meshgrid(
        np.asarray(x_da.time.values),
        np.asarray(x_da.lat.values),
        np.asarray(x_da.lon.values),
        indexing="ij",
    )
    index = pd.MultiIndex.from_arrays(
        [times.reshape(-1), lats.reshape(-1), lons.reshape(-1)],
        names=["time", "lat", "lon"],
    )
    names = [str(v) for v in x_da["feature"].values]
    table = SampleTable(X, y_vals, names, index, y_name=str(y_name))
    finite = np.isfinite(table.X).all(axis=1) & np.isfinite(table.y)
    extras = {
        "template": y_da,
        "n_raw": len(table.y),
        "n_valid": int(finite.sum()),
    }
    return SampleTable(table.X[finite], table.y[finite], names, table.index[finite], table.y_name, extras)


def table_to_cube(
    values: np.ndarray,
    index: pd.MultiIndex,
    template: xr.DataArray,
    name: str = "prediction",
) -> xr.DataArray:
    series = pd.Series(np.asarray(values, dtype=np.float64).reshape(-1), index=index)
    da = series[~series.index.duplicated(keep="last")].to_xarray()
    da = da.reindex(time=template.time, lat=template.lat, lon=template.lon)
    da = da.transpose("time", "lat", "lon")
    da.name = name
    return da


def values_to_cubes(
    matrix: np.ndarray,
    index: pd.MultiIndex,
    template: xr.DataArray,
    names: Sequence[str],
) -> xr.Dataset:
    data_vars = {}
    matrix = np.asarray(matrix)
    for i, name in enumerate(names):
        data_vars[str(name)] = table_to_cube(matrix[:, i], index, template, name=str(name))
    return xr.Dataset(data_vars)
