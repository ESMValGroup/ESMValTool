from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd


def ice_pdp(
    model,
    x: np.ndarray,
    names: Sequence[str],
    n_grid: int = 20,
    ice_n_samples: int = 50,
    seed: int = 42,
) -> dict[str, dict]:
    x = np.asarray(x, dtype=np.float64)
    rng = np.random.default_rng(seed)
    n = len(x)
    ice_idx = rng.choice(n, size=min(ice_n_samples, n), replace=False)
    out = {}
    for j, name in enumerate(names):
        grid = np.unique(np.nanquantile(x[:, j], np.linspace(0.02, 0.98, n_grid)))
        ice = np.zeros((len(ice_idx), len(grid)))
        for k, value in enumerate(grid):
            xx = x[ice_idx].copy()
            xx[:, j] = value
            ice[:, k] = model.predict(xx)
        pdp = ice.mean(axis=0)
        out[str(name)] = {
            "grid": grid,
            "pdp": pdp,
            "ice": ice,
            "ice_index": ice_idx,
        }
    return out
