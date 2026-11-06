from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from rf_shap.cube.tabular import SampleTable


@dataclass
class SplitResult:
    train: SampleTable
    test: SampleTable
    method: str
    extras: dict


def _mask_split(table: SampleTable, train_mask: np.ndarray, method: str, extras: dict | None = None) -> SplitResult:
    train_mask = np.asarray(train_mask, dtype=bool)
    if train_mask.all() or (~train_mask).all():
        raise ValueError("Split produced empty train or test set")
    return SplitResult(table.subset(train_mask), table.subset(~train_mask), method, extras or {})


def _sorted_times(index: pd.MultiIndex) -> np.ndarray:
    times = np.asarray(index.get_level_values("time"))
    uniq = np.unique(times)
    try:
        return np.sort(uniq)
    except TypeError:
        return uniq


def time_split(table: SampleTable, test_size: float = 0.2, train_end=None) -> SplitResult:
    times = pd.Series(table.index.get_level_values("time"))
    uniq = pd.Index(_sorted_times(table.index))
    if train_end is not None:
        cutoff = pd.Timestamp(train_end) if np.issubdtype(uniq.to_numpy().dtype, np.datetime64) else train_end
        train_mask = times <= cutoff
    else:
        cut_i = max(1, int(round(len(uniq) * (1 - test_size)))) - 1
        cutoff = uniq[cut_i]
        train_mask = times <= cutoff
    return _mask_split(table, train_mask.to_numpy(), "time", {"cutoff": cutoff})


def random_split(table: SampleTable, test_size: float = 0.2, seed: int = 42) -> SplitResult:
    rng = np.random.default_rng(seed)
    n = table.n_samples
    test_n = max(1, int(round(n * test_size)))
    test_idx = rng.choice(n, size=test_n, replace=False)
    train_mask = np.ones(n, dtype=bool)
    train_mask[test_idx] = False
    return _mask_split(table, train_mask, "random")


def space_split(table: SampleTable, mode: str = "lon_band", n_blocks: int = 4, test_size: float = 0.2, seed: int = 42) -> SplitResult:
    lat = np.asarray(table.index.get_level_values("lat"), dtype=float)
    lon = np.asarray(table.index.get_level_values("lon"), dtype=float)
    mode = mode.lower()
    if mode == "lon_band":
        q = np.quantile(lon, 1 - test_size)
        train_mask = lon < q
    elif mode == "lat_band":
        q = np.quantile(lat, 1 - test_size)
        train_mask = lat < q
    elif mode == "checkerboard":
        lat_bin = np.digitize(lat, np.linspace(lat.min(), lat.max(), n_blocks + 1)[1:-1])
        lon_bin = np.digitize(lon, np.linspace(lon.min(), lon.max(), n_blocks + 1)[1:-1])
        train_mask = ((lat_bin + lon_bin) % 2) == 0
    elif mode == "block":
        lat_bin = np.digitize(lat, np.linspace(lat.min(), lat.max(), n_blocks + 1)[1:-1])
        lon_bin = np.digitize(lon, np.linspace(lon.min(), lon.max(), n_blocks + 1)[1:-1])
        block_id = lat_bin * n_blocks + lon_bin
        uniq = np.unique(block_id)
        rng = np.random.default_rng(seed)
        n_test = max(1, int(round(len(uniq) * test_size)))
        test_blocks = set(rng.choice(uniq, size=n_test, replace=False))
        train_mask = np.array([b not in test_blocks for b in block_id])
    else:
        raise ValueError(f"Unknown spatial split mode: {mode}")
    return _mask_split(table, train_mask, f"space:{mode}")


def extreme_split(table: SampleTable, quantile: float = 0.9, hold_extremes_as: str = "test") -> SplitResult:
    y = table.y
    lo, hi = np.nanquantile(y, [1 - quantile, quantile])
    extreme = (y <= lo) | (y >= hi)
    if hold_extremes_as == "test":
        train_mask = ~extreme
    elif hold_extremes_as == "train":
        train_mask = extreme
    else:
        raise ValueError("hold_extremes_as must be 'test' or 'train'")
    return _mask_split(table, train_mask, "extreme", {"lo": float(lo), "hi": float(hi)})


def time_window_splits(
    table: SampleTable,
    n_splits: int = 3,
    test_span: int = 12,
    expanding: bool = True,
) -> list[SplitResult]:
    uniq = _sorted_times(table.index)
    times = np.asarray(table.index.get_level_values("time"))
    if len(uniq) <= test_span + 1:
        raise ValueError("Not enough time steps for a sliding-window split")
    results = []
    max_start = len(uniq) - test_span
    starts = np.linspace(test_span, max_start, n_splits).astype(int)
    for start in starts:
        test_times = uniq[start : start + test_span]
        if expanding:
            train_times = uniq[:start]
        else:
            train_times = uniq[max(0, start - test_span) : start]
        train_mask = np.isin(times, train_times)
        test_mask = np.isin(times, test_times)
        if train_mask.sum() == 0 or test_mask.sum() == 0:
            continue
        results.append(
            SplitResult(
                table.subset(train_mask),
                table.subset(test_mask),
                "time_window",
                {"test_start": test_times[0], "test_end": test_times[-1], "expanding": expanding},
            )
        )
    if not results:
        raise ValueError("Sliding window produced no valid folds")
    return results


def split_table(table: SampleTable, options: dict | None = None, seed: int = 42):
    options = options or {}
    method = str(options.get("method", "time")).lower()
    test_size = float(options.get("test_size", 0.2))
    if method == "time":
        cfg = options.get("time", {})
        return time_split(table, test_size=test_size, train_end=cfg.get("train_end"))
    if method == "random":
        return random_split(table, test_size=test_size, seed=seed)
    if method == "space":
        cfg = options.get("space", {})
        return space_split(
            table,
            mode=cfg.get("mode", "lon_band"),
            n_blocks=int(cfg.get("n_blocks", 4)),
            test_size=test_size,
            seed=seed,
        )
    if method == "extreme":
        cfg = options.get("extreme", {})
        return extreme_split(
            table,
            quantile=float(cfg.get("quantile", 0.9)),
            hold_extremes_as=cfg.get("hold_extremes_as", "test"),
        )
    if method in {"time_window", "sliding"}:
        cfg = options.get("time_window", {})
        return time_window_splits(
            table,
            n_splits=int(cfg.get("n_splits", 3)),
            test_span=int(cfg.get("test_span", 12)),
            expanding=bool(cfg.get("expanding", True)),
        )
    raise ValueError(f"Unknown split method: {method}")
