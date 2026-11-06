from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr


def make_synthetic_cubes(
    n_time: int = 24,
    n_lat: int = 8,
    n_lon: int = 10,
    seed: int = 0,
) -> tuple[xr.Dataset, xr.Dataset]:
    rng = np.random.default_rng(seed)
    time = pd.date_range("2000-01-01", periods=n_time, freq="MS")
    lat = np.linspace(10, 50, n_lat)
    lon = np.linspace(70, 130, n_lon)
    t = np.arange(n_time)[:, None, None]
    la = lat[None, :, None]
    lo = lon[None, None, :]
    tas = 10 + 8 * np.sin(2 * np.pi * t / 12) + 0.02 * la + 0.1 * rng.normal(size=(n_time, n_lat, n_lon))
    pr = 2 + 1.2 * np.cos(2 * np.pi * t / 12) + 0.01 * (lo - 100) + 0.1 * rng.normal(size=(n_time, n_lat, n_lon))
    huss = 5 + 0.3 * tas / 10 + 0.05 * rng.normal(size=(n_time, n_lat, n_lon))
    # observation as a nonlinear mix of model-like inputs
    obs = 0.6 * tas + 0.8 * pr + 0.2 * tas * pr / 20 - 0.15 * huss + 0.2 * rng.normal(size=tas.shape)
    x = xr.Dataset(
        {
            "tas": (("time", "lat", "lon"), tas),
            "pr": (("time", "lat", "lon"), pr),
            "huss": (("time", "lat", "lon"), huss),
        },
        coords={"time": time, "lat": lat, "lon": lon},
    )
    y = xr.Dataset({"obs": (("time", "lat", "lon"), obs)}, coords=x.coords)
    return x, y
