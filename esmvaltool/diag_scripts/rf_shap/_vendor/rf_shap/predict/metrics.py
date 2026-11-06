from __future__ import annotations

import numpy as np
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def regression_metrics(y_true, y_pred) -> dict[str, float]:
    yt = np.asarray(y_true, dtype=np.float64).reshape(-1)
    yp = np.asarray(y_pred, dtype=np.float64).reshape(-1)
    mask = np.isfinite(yt) & np.isfinite(yp)
    yt, yp = yt[mask], yp[mask]
    if len(yt) == 0:
        raise ValueError("No finite prediction pairs")
    rmse = float(np.sqrt(mean_squared_error(yt, yp)))
    bias = float(np.mean(yp - yt))
    ubrmse = float(np.sqrt(max(rmse**2 - bias**2, 0.0)))
    nse_den = np.sum((yt - yt.mean()) ** 2)
    nse = float(1 - np.sum((yp - yt) ** 2) / nse_den) if nse_den > 0 else np.nan
    r = _corr(yt, yp)
    alpha = float(yp.std() / yt.std()) if yt.std() > 0 else np.nan
    beta = float(yp.mean() / yt.mean()) if yt.mean() != 0 else np.nan
    kge = float(1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)) if np.isfinite(r) else np.nan
    mape_den = np.where(yt == 0, np.nan, yt)
    mape = float(np.nanmean(np.abs((yp - yt) / mape_den)))
    return {
        "n": int(len(yt)),
        "r2": float(r2_score(yt, yp)),
        "rmse": rmse,
        "mae": float(mean_absolute_error(yt, yp)),
        "bias": bias,
        "ubrmse": ubrmse,
        "nse": nse,
        "kge": kge,
        "pearson_r": r,
        "mape": mape,
    }
