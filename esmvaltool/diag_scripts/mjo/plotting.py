"""Plotting utilities for the MJO Hovmöller diagnostic."""

import logging

import matplotlib.pyplot as plt
import numpy as np

logger = logging.getLogger(__name__)


def plot_hovmoeller(
    regression,
    longitude,
    lags,
    cfg,
    title=None,
):
    """Create and return a longitude-lag Hovmöller figure."""
    regression = np.ma.masked_invalid(regression)
    longitude = np.asarray(longitude)
    lags = np.asarray(lags)

    if regression.shape != (lags.size, longitude.size):
        msg = (
            "Regression shape must be (lag, longitude), but got "
            f"{regression.shape} for {lags.size} lags and "
            f"{longitude.size} longitudes."
        )
        raise ValueError(msg)

    # Ensure longitudes and corresponding data are ordered.
    order = np.argsort(longitude)
    longitude = longitude[order]
    regression = regression[:, order]

    vmax = cfg.get("vmax")
    if vmax is None:
        vmax = np.ma.max(np.ma.abs(regression))
        if np.ma.is_masked(vmax) or not np.isfinite(vmax) or vmax == 0:
            vmax = 1.0

    nlevels = int(cfg.get("contour_levels", 21))
    if nlevels < 3:
        msg = "contour_levels must be at least 3"
        raise ValueError(msg)

    levels = np.linspace(-float(vmax), float(vmax), nlevels)

    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)

    contour = ax.contourf(
        longitude,
        lags,
        regression,
        levels=levels,
        cmap=cfg.get("colormap", "RdYlBu"),
        extend="both",
    )

    ax.axhline(0.0, color="black", linewidth=0.8)

    lon0, lon1 = cfg["reference_longitudes"]
    for ref_lon in (lon0, lon1):
        ax.axvline(
            ref_lon,
            color="black",
            linestyle="--",
            linewidth=0.8,
        )

    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Lag (days)")
    ax.set_xlim(cfg.get("longitude_limits", [0.0, 360.0]))

    ax.set_title(title or cfg.get("plot_title", "MJO Hovmöller diagram"))

    colorbar = fig.colorbar(contour, ax=ax, pad=0.02)
    colorbar.set_label(
        cfg.get(
            "colorbar_label",
            "Precipitation regression coefficient",
        ),
    )

    return fig
