"""MJO lag-regression Hovmöller diagnostic.

Description
-----------
This diagnostic computes lag-regression Hovmöller diagrams of the
Madden-Julian Oscillation (MJO). Its input is a daily time-longitude
field of tropical precipitation anomalies: already averaged over a
latitude band and with the annual cycle removed by the recipe's
preprocessor (`extract_region` + `meridional_statistics` +
`anomalies: {period: day}`). These anomalies are Lanczos band-pass
filtered to the MJO period range. A reference index is built by
averaging the filtered field over a longitude sector, and the full
filtered field is then regressed against this index at a range of lags.
The resulting lag-longitude regression field is saved as NetCDF and
plotted as a Hovmöller diagram, in which an eastward-propagating
diagonal band is the signature of MJO convection.

Caveats
-------
The diagnostic assumes that every input dataset has already been
reduced to a single daily time step, averaged over the tropical
latitude band, and had its day-of-year climatology removed by the
recipe's preprocessor (e.g., via `daily_statistics: mean`,
`meridional_statistics: mean`, and `anomalies: {period: day}`); none of
this is checked on the cube itself. If a dataset is not truly daily, the
band-pass filter periods and the lag axis will be silently
misinterpreted. The latitude band used in captions and provenance
attributes is read directly from the cube's scalar `latitude` coordinate
(as left behind by `meridional_statistics`), so it always matches what
the preprocessor actually did. Exactly one preprocessed file is expected
per dataset (grouped by the `alias` facet).

Author
------
Jurij Schönfeld (DLR, Germany)

Configuration options in recipe
--------------------------------
colorbar_label: str, optional (default: 'Precipitation regression
    coefficient')
    Label for the figure's colorbar.
colormap: str, optional (default: 'RdYlBu')
    Matplotlib colormap used for the Hovmöller contour plot.
contour_levels: int, optional (default: 21)
    Number of contour levels in the Hovmöller plot. Must be at least 3.
high_period: int
    Upper period cutoff (in days) of the Lanczos band-pass filter.
lanczos_weights: int
    Number of weights of the Lanczos band-pass filter. Must be an odd
    integer greater than 1.
longitude_limits: list of float, optional (default: [0.0, 360.0])
    Longitude axis limits of the Hovmöller plot.
low_period: int
    Lower period cutoff (in days) of the Lanczos band-pass filter.
max_lag: int
    Maximum lag (in days, in both directions) computed by the lag
    regression.
plot_title: str, optional (default: 'MJO Hovmöller diagram')
    Title of the Hovmöller plot.
reference_longitudes: list of float
    Longitude sector (`[lon0, lon1]`) used to build the MJO reference
    index that the filtered field is regressed against.
"""

import logging
from pathlib import Path

import iris
import iris.analysis

from esmvaltool.diag_scripts.mjo.filtering import lanczos_bandpass
from esmvaltool.diag_scripts.mjo.plotting import plot_hovmoeller
from esmvaltool.diag_scripts.mjo.regression import (
    lag_regression,
    regression_cube,
)
from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
)

logger = logging.getLogger(Path(__file__).stem)


def compute_mjo_index(cube, longitude_bounds):
    """Average the tropical field over the reference longitude sector."""
    region = cube.intersection(longitude=tuple(longitude_bounds))
    return region.collapsed("longitude", iris.analysis.MEAN).data


def process_tropical_cube(cube, cfg):
    """Process a daily time-longitude tropical anomaly field."""
    filtered = lanczos_bandpass(
        cube,
        low_period=cfg["low_period"],
        high_period=cfg["high_period"],
        weights=cfg["lanczos_weights"],
    )

    index = compute_mjo_index(filtered, cfg["reference_longitudes"])
    regression, lags = lag_regression(
        filtered,
        index,
        max_lag=cfg["max_lag"],
    )

    result = regression_cube(
        regression,
        lags,
        filtered.coord("longitude"),
        filtered.units,
        latitude_coord=filtered.coord("latitude"),
    )
    result.attributes.update(
        {
            "reference_longitudes": str(cfg["reference_longitudes"]),
            "filter_period": (
                f"{cfg['low_period']}-{cfg['high_period']} days"
            ),
            "lanczos_weights": cfg["lanczos_weights"],
        },
    )
    return result


def get_provenance_record(result, label, ancestors, cfg):
    """Create provenance shared by the NetCDF and figure."""
    lon0, lon1 = cfg["reference_longitudes"]
    lat0, lat1 = result.coord("latitude").bounds[0]

    return {
        "caption": (
            f"{label}: lagged regression of {cfg['low_period']}- to "
            f"{cfg['high_period']}-day filtered precipitation against the "
            f"precipitation index averaged over {lon0:g}-{lon1:g}E and "
            f"{abs(lat0):g}S-{abs(lat1):g}N."
        ),
        "statistics": ["mean", "other"],
        "domains": ["trop"],
        "plot_types": ["zonal"],
        "authors": ["schoenfeld_jurij"],
        "references": ["hannah20james"],
        "ancestors": ancestors,
    }


def save_result(result, label, ancestors, cfg):
    """Save one regression field and its Hovmöller figure."""
    safe_label = label.replace(" ", "_").replace("/", "_")
    basename = f"{safe_label}_mjo_hovmoeller"

    result.attributes["dataset"] = label
    provenance = get_provenance_record(result, label, ancestors, cfg)

    save_data(basename, provenance, cfg, result)

    figure = plot_hovmoeller(
        result.data,
        result.coord("longitude").points,
        result.coord("lag").points,
        cfg,
        title=f"{label}: MJO precipitation lag regression",
    )

    save_figure(
        basename,
        provenance,
        cfg,
        figure=figure,
        dpi=300,
        bbox_inches="tight",
    )


def main(cfg):
    """Run the diagnostic."""
    groups = group_metadata(
        cfg["input_data"].values(),
        "alias",
        sort="dataset",
    )

    for label, entries in groups.items():
        if len(entries) != 1:
            msg = f"Expected one input file for {label!r}, got {len(entries)}"
            raise ValueError(msg)

        metadata = entries[0]
        logger.info("Processing %s", label)

        cube = iris.load_cube(metadata["filename"])
        result = process_tropical_cube(cube, cfg)

        save_result(
            result,
            label,
            [metadata["filename"]],
            cfg,
        )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
