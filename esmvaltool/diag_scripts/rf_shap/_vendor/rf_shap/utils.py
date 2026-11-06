from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
import pandas as pd


def ensure_1d(arr) -> np.ndarray:
    return np.asarray(arr).reshape(-1)


def as_float(arr) -> np.ndarray:
    return np.asarray(arr, dtype=np.float64)


def valid_mask(*arrays: np.ndarray) -> np.ndarray:
    mask = np.ones(len(arrays[0]), dtype=bool)
    for arr in arrays:
        values = np.asarray(arr)
        if values.ndim == 1:
            mask &= np.isfinite(values)
        else:
            mask &= np.isfinite(values).all(axis=1)
    return mask


def subset_index(index: pd.MultiIndex, mask: np.ndarray) -> pd.MultiIndex:
    return index[np.asarray(mask, dtype=bool)]


def check_names(names: Sequence[str], available: Iterable[str], kind: str) -> list[str]:
    available_set = set(available)
    missing = [name for name in names if name not in available_set]
    if missing:
        raise KeyError(f"{kind} missing variables: {missing}")
    return list(names)
