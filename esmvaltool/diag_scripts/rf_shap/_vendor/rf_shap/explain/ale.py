from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd


def _grid(x_col: np.ndarray, n_grid: int) -> np.ndarray:
    qs = np.linspace(0.02, 0.98, n_grid)
    return np.unique(np.nanquantile(x_col, qs))


def ale_1d(model, x: np.ndarray, names: Sequence[str], n_grid: int = 20) -> dict[str, pd.DataFrame]:
    """1D Accumulated Local Effects for each feature."""
    x = np.asarray(x, dtype=np.float64)
    out = {}
    n = len(x)
    for j, name in enumerate(names):
        grid = _grid(x[:, j], n_grid)
        if len(grid) < 2:
            continue
        effects = [0.0]
        centers = []
        for lo, hi in zip(grid[:-1], grid[1:]):
            mask = (x[:, j] >= lo) & (x[:, j] <= hi)
            if mask.sum() == 0:
                effects.append(effects[-1])
                centers.append(0.5 * (lo + hi))
                continue
            left = x[mask].copy()
            right = x[mask].copy()
            left[:, j] = lo
            right[:, j] = hi
            delta = model.predict(right) - model.predict(left)
            effects.append(effects[-1] + float(np.mean(delta)))
            centers.append(0.5 * (lo + hi))
        values = np.asarray(effects[1:])
        values = values - np.mean(values)
        out[str(name)] = pd.DataFrame({"x": centers, "ale": values})
    return out
