"""Antarctic sea-ice seasonality (COSIMA Sea_Ice_Seasonality_Statistics).

Translation of ``03-Advanced-Recipes/Sea_Ice_Seasonality_Statistics.ipynb``
onto CMORised ``SIday`` sea-ice concentration, following Massom et al.
(2013):

- a cell is *advancing* on the first day of the first run of five
  consecutive days at or above 15% concentration;
- it is *retreating* on the last day at or above 15%;
- the *season duration* is the difference.

A sea-ice year runs 15 February to 14 February, so the recipe selects
that window and the day numbers below are relative to 15 February.

Cells that never reach the threshold are masked; cells that never drop
below it get the full year.  The five-day rule is what distinguishes a
genuine advance from a transient excursion, and the diagnostic keeps it
rather than the simpler "first day above threshold" it is sometimes
confused with.
"""

import logging
import os

import iris
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from esmvaltool.diag_scripts.shared import (  # noqa: E402
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
)

from esmvaltool.diag_scripts.cosima import cosima_common as cc

logger = logging.getLogger(os.path.basename(__file__))

PANELS = (
    ("advance", "Sea ice advance day", "twilight"),
    ("retreat", "Sea ice retreat day", "twilight"),
    ("duration", "Sea ice season duration", "Blues"),
)


def seasonality_fields(cube, threshold=0.15, consecutive_days=5):
    """Advance, retreat and duration from a daily concentration cube."""
    concentration = cc.fraction_from_percent(cube)
    if concentration.ndim != 3:
        raise ValueError(
            f"expected (time, lat, lon) daily concentration, got shape "
            f"{concentration.shape}")
    return cc.seaice_seasonality(
        concentration, threshold=threshold,
        consecutive_days=consecutive_days)


def _plot_panel(figure, position, field, cube, title, cmap, max_lat):
    """One circumpolar panel."""
    import cartopy.crs as ccrs

    axes = cc.south_polar_axes(figure, position, max_lat=max_lat)
    lon = cc.lon_2d(cube)
    lat = cc.lat_2d(cube)
    mesh = axes.pcolormesh(lon, lat, field, cmap=cmap,
                           transform=ccrs.PlateCarree(), shading="auto")
    axes.coastlines(resolution="110m", linewidth=0.5)
    axes.set_title(title, fontsize=11)
    figure.colorbar(mesh, ax=axes, orientation="vertical", shrink=0.65,
                    pad=0.06, label="Days since 15 February")
    return axes


def main(cfg):
    """Run the diagnostic."""
    threshold = float(cfg.get("concentration_threshold", 0.15))
    consecutive = int(cfg.get("consecutive_days", 5))
    max_lat = float(cfg.get("max_latitude", -50.0))
    grouped = group_metadata(cfg["input_data"].values(), "dataset")

    for dataset, metadata in sorted(grouped.items()):
        meta = metadata[0]
        cube = cc.load_cube(meta["filename"], meta["short_name"])
        try:
            advance, retreat, duration = seasonality_fields(
                cube, threshold=threshold, consecutive_days=consecutive)
        except ValueError as error:
            logger.warning("%s: %s", dataset, error)
            continue

        fields = {"advance": advance, "retreat": retreat,
                  "duration": duration}

        figure = plt.figure(figsize=(16, 5.5))
        for index, (name, title, cmap) in enumerate(PANELS):
            _plot_panel(figure, (1, 3, index + 1), fields[name], cube,
                        f"{title}\n{dataset}", cmap, max_lat)
        figure.suptitle(
            f"Antarctic sea-ice seasonality (Massom et al. 2013): {dataset}"
            f"\nthreshold {threshold:.0%}, {consecutive} consecutive days",
            fontsize=12)
        plot_path = get_plot_filename(
            f"cosima_seaice_seasonality_{dataset}", cfg)
        figure.savefig(plot_path, bbox_inches="tight", dpi=150)
        plt.close(figure)
        logger.info("Wrote %s", plot_path)

        lat = cc.lat_1d(cube)
        lon = cc.lon_1d(cube)
        cubes = iris.cube.CubeList([
            cc.make_cube(
                np.ma.filled(fields[name], np.nan),
                [(lat, "latitude", "degrees_north"),
                 (lon, "longitude", "degrees_east")],
                f"siseason_{name}", "days",
                f"sea ice {name} day relative to 15 February")
            for name, _, _ in PANELS])
        path = get_diagnostic_filename(
            f"cosima_seaice_seasonality_{dataset}", cfg)
        iris.save(cubes, path)
        record = {
            "caption": f"Sea-ice advance, retreat and season duration for "
                       f"{dataset}, after Massom et al. (2013).",
            "statistics": ["other"],
            "domains": ["shpolar"],
            "plot_types": ["polar"],
            "authors": ["beucher_romain"],
            "references": ["cosima_recipes", "massom13plosone"],
            "ancestors": [meta["filename"]],
        }
        with ProvenanceLogger(cfg) as logger_:
            logger_.log(path, record)
            logger_.log(plot_path, record)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
