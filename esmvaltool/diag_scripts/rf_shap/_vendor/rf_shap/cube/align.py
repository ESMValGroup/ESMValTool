from __future__ import annotations

from typing import Literal

import numpy as np
import pandas as pd
import xarray as xr


RegridMode = Literal["none", "nearest", "linear"]
JoinTime = Literal["inner", "nearest"]


def stamp_year_month(data: xr.Dataset | xr.DataArray) -> xr.Dataset | xr.DataArray:
    """Replace the time axis with the first day of each calendar month.

    Yearly observational files often use mid-month dates or timestamps with a
    few stray nanoseconds. Mapping them onto ``YYYY-MM-01`` lets cubes join
    without interpolation. Duplicate months are averaged.
    """
    if "time" not in data.coords and "time" not in data.dims:
        raise ValueError("Need a time coordinate to stamp year-month")
    stamped = pd.to_datetime(data.time.values).to_period("M").to_timestamp()
    out = data.assign_coords(time=("time", np.asarray(stamped, dtype="datetime64[ns]")))
    if pd.Index(out.time.values).duplicated().any():
        out = out.groupby("time").mean(keep_attrs=True)
    return out.sortby("time")


def align_month_cubes(*cubes: xr.Dataset, join: str = "inner") -> tuple[xr.Dataset, ...]:
    """Inner (or outer) join of already-gridded cubes on year-month time."""
    stamped = [stamp_year_month(cube) for cube in cubes]
    return tuple(xr.align(*stamped, join=join, copy=False))


def _interp_to(source: xr.Dataset, target: xr.Dataset, method: RegridMode) -> xr.Dataset:
    if method == "none":
        return source
    kwargs = dict(lat=target.lat, lon=target.lon)
    if method == "nearest":
        return source.interp(lat=target.lat, lon=target.lon, method="nearest")
    return source.interp(**kwargs, method="linear")


def _align_time(x: xr.Dataset, y: xr.Dataset, how: JoinTime) -> tuple[xr.Dataset, xr.Dataset]:
    if how == "inner":
        common = np.intersect1d(x.time.values, y.time.values)
        if common.size == 0:
            raise ValueError("X and Y have no overlapping time stamps")
        return x.sel(time=common), y.sel(time=common)
    y_aligned = y.reindex(time=x.time, method="nearest")
    return x, y_aligned


def align_xy(
    x: xr.Dataset,
    y: xr.Dataset,
    regrid: RegridMode = "nearest",
    join_time: JoinTime = "inner",
    target: str = "y",
) -> tuple[xr.Dataset, xr.Dataset]:
    """Align X and Y onto the same (time, lat, lon) grid.

    By default Y (observations) is the target grid, matching typical
    CMIP6-vs-obs diagnosis. Set target='x' to regrid observations to the model grid.
    """
    if target not in {"x", "y"}:
        raise ValueError("target must be 'x' or 'y'")
    ref, other = (x, y) if target == "x" else (y, x)
    other_r = _interp_to(other, ref, regrid)
    if target == "x":
        x_a, y_a = ref, other_r
    else:
        y_a, x_a = ref, other_r
    x_a, y_a = _align_time(x_a, y_a, join_time)
    x_a, y_a = xr.align(x_a, y_a, join="inner", copy=False)
    if x_a.sizes != y_a.sizes:
        raise ValueError(f"Aligned cubes still differ: X={dict(x_a.sizes)} Y={dict(y_a.sizes)}")
    return x_a, y_a
