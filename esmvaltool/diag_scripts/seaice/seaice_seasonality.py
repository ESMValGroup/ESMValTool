# Copyright (C) 2026 ESMValTool contributors.
"""Map Antarctic sea ice advance, retreat, and season duration."""

import logging
from itertools import pairwise

import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    run_diagnostic,
)

logger = logging.getLogger(__name__)

FIELD_NAMES = ("advance", "retreat", "duration")
TITLES = ("Advance", "Retreat", "Season duration")
N_DIMENSIONS = 3
CONCENTRATION_TOLERANCE = 1e-5


def validate_daily_ice_year(cube):
    """Require one complete 15 February to 14 February daily time axis."""
    if cube.ndim != N_DIMENSIONS or cube.coord_dims("time") != (0,):
        message = "siconc must have dimensions (time, latitude, longitude)"
        raise ValueError(message)
    dates = cube.coord("time").units.num2date(cube.coord("time").points)
    first, last = dates[0], dates[-1]
    if (first.month, first.day) != (2, 15) or (
        last.year,
        last.month,
        last.day,
    ) != (first.year + 1, 2, 14):
        message = "siconc must span 15 February to 14 February"
        raise ValueError(message)
    if any(
        not np.isclose((right - left).total_seconds(), 86400)
        for left, right in pairwise(dates)
    ):
        message = "siconc must have one sample on every day"
        raise ValueError(message)
    return len(dates)


def concentration_fraction(cube):
    """Return masked concentration fractions from CMOR percentages."""
    data = np.ma.masked_invalid(np.ma.asarray(cube.data, dtype=float))
    units = str(cube.units).strip().lower()
    if units in ("%", "percent"):
        data = data / 100.0
    elif units != "1":
        message = f"siconc must have units '%' or '1', got {cube.units}"
        raise ValueError(message)
    values = data.compressed()
    if np.any(
        (values < -CONCENTRATION_TOLERANCE)
        | (values > 1 + CONCENTRATION_TOLERANCE),
    ):
        message = "siconc has values outside the physical range 0-1"
        raise ValueError(message)
    return data


def seasonality_fields(concentration, threshold=0.15, consecutive_days=5):
    """Compute day of advance, retreat, and duration for one ice year.

    Advance is the first day in the first sustained run above threshold.
    Retreat is the first day below threshold after the last ice-covered
    day, capped at the end of the ice year. Days are one-based relative
    to 15 February. Cells without a sustained advance or with missing
    daily samples are masked in all three outputs.
    """
    if concentration.ndim != N_DIMENSIONS:
        message = "concentration must be (time, latitude, longitude)"
        raise ValueError(message)
    if not 0 < threshold < 1 or consecutive_days < 1:
        message = "threshold and consecutive_days must be positive"
        raise ValueError(message)
    n_days = concentration.shape[0]
    if consecutive_days > n_days:
        message = "consecutive_days exceeds the length of the ice year"
        raise ValueError(message)

    data = np.ma.masked_invalid(np.ma.asarray(concentration, dtype=float))
    valid = ~np.ma.getmaskarray(data).any(axis=0)
    above = np.ma.filled(data >= threshold, fill_value=False)
    run = np.zeros(above.shape[1:], dtype=np.int16)
    advance = np.full(above.shape[1:], np.nan)
    last_ice = np.zeros(above.shape[1:], dtype=np.int16)
    for day, covered in enumerate(above, start=1):
        run = np.where(covered, run + 1, 0)
        new_advance = np.isnan(advance) & (run == consecutive_days)
        advance[new_advance] = day - consecutive_days + 1
        last_ice[covered] = day

    retreat = np.minimum(last_ice + 1, n_days).astype(float)
    perennial = above.all(axis=0)
    advance[perennial] = 1
    retreat[perennial] = n_days
    mask = ~valid | np.isnan(advance)
    advance = np.ma.masked_array(advance, mask=mask)
    retreat = np.ma.masked_array(retreat, mask=mask)
    duration = retreat - advance
    return advance, retreat, duration


def output_cube(source, field, name):
    """Preserve the regridded horizontal coordinates on a 2D output."""
    cube = source[0].copy(data=field)
    cube.remove_coord("time")
    cube.standard_name = None
    cube.var_name = f"siseason_{name}"
    cube.long_name = f"sea ice {name} day relative to 15 February"
    cube.units = "days"
    return cube


def plot_comparison(results, cfg, ancestors):
    """Plot every model on matching polar maps and colour scales."""
    n_models = len(results)
    fig, axes = plt.subplots(
        n_models,
        3,
        figsize=(14, 4.3 * n_models),
        squeeze=False,
        subplot_kw={"projection": ccrs.SouthPolarStereo()},
    )
    for row, (label, cube, fields) in enumerate(results):
        lon, lat = np.meshgrid(
            cube.coord("longitude").points,
            cube.coord("latitude").points,
        )
        for col, (name, title, field) in enumerate(
            zip(FIELD_NAMES, TITLES, fields, strict=True),
        ):
            ax = axes[row, col]
            ax.set_extent([-180, 180, -90, -50], ccrs.PlateCarree())
            ax.gridlines(linewidth=0.4, alpha=0.5)
            mesh = ax.pcolormesh(
                lon,
                lat,
                field,
                transform=ccrs.PlateCarree(),
                shading="auto",
                cmap="Blues" if name == "duration" else "viridis",
                vmin=0,
                vmax=366,
            )
            ax.set_title(f"{label}: {title}")
            fig.colorbar(
                mesh,
                ax=ax,
                shrink=0.7,
                label="Days since 15 February",
            )
    fig.suptitle("Antarctic sea ice seasonality", fontsize=14)
    path = get_plot_filename("seaice_seasonality_comparison", cfg)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    with ProvenanceLogger(cfg) as provenance:
        provenance.log(
            path,
            {
                "caption": "Model comparison of Antarctic sea ice advance, "
                "retreat, and season duration.",
                "statistics": ["other"],
                "domains": ["shpolar"],
                "plot_types": ["polar"],
                "authors": ["beucher_romain"],
                "references": ["massom13plosone"],
                "ancestors": ancestors,
            },
        )


def main(cfg):
    """Run the seasonality analysis for every dataset in the recipe."""
    threshold = float(cfg.get("concentration_threshold", 0.15))
    consecutive_days = int(cfg.get("consecutive_days", 5))
    results = []
    ancestors = []
    metadata = sorted(
        cfg["input_data"].values(),
        key=lambda item: (item["dataset"], item.get("ensemble", "")),
    )
    for item in metadata:
        cube = iris.load_cube(item["filename"])
        validate_daily_ice_year(cube)
        concentration = concentration_fraction(cube)
        fields = seasonality_fields(
            concentration,
            threshold=threshold,
            consecutive_days=consecutive_days,
        )
        label = f"{item['dataset']} {item.get('ensemble', '')}".strip()
        logger.info("Computed sea ice seasonality for %s", label)
        output = iris.cube.CubeList(
            output_cube(cube, field, name)
            for name, field in zip(FIELD_NAMES, fields, strict=True)
        )
        path = get_diagnostic_filename(
            f"seaice_seasonality_{item['dataset']}_{item.get('ensemble', '')}",
            cfg,
        )
        iris.save(output, path)
        with ProvenanceLogger(cfg) as provenance:
            provenance.log(
                path,
                {
                    "caption": f"Antarctic sea ice seasonality for {label}.",
                    "statistics": ["other"],
                    "domains": ["shpolar"],
                    "plot_types": ["geo"],
                    "authors": ["beucher_romain"],
                    "references": ["massom13plosone"],
                    "ancestors": [item["filename"]],
                },
            )
        results.append((label, cube, fields))
        ancestors.append(item["filename"])
    plot_comparison(results, cfg, ancestors)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
