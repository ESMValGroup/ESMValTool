# (C) Crown Copyright 2026, the Met Office.
"""Cloud-type feedbacks from anomaly regression.

This diagnostic computes cloud-fraction feedbacks by cloud type as the linear
regression slope of deseasonalized `clisccp` anomalies (at each grid point)
against deseasonalized global-mean `tas` anomalies.
"""

import logging
import os

import cartopy.crs as ccrs
import iris
import iris.analysis
import iris.analysis.cartography
import iris.plot as iplt
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    run_diagnostic,
)

logger = logging.getLogger(os.path.basename(__file__))


def _get_input_metadata(cfg, short_name):
    """Return the single metadata entry for a variable short name."""
    matches = [
        d for d in cfg["input_data"].values() if d["short_name"] == short_name
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one '{short_name}' dataset, got {len(matches)}"
        )
    return matches[0]


def _extract_region_mean(cube, lat_min, lat_max, lon_min=None, lon_max=None):
    """Return area-weighted mean over a lat/lon box.

    Parameters
    ----------
    cube : iris.cube.Cube
        Input cube with latitude and longitude dimension coordinates.
    lat_min, lat_max : float
        Latitude bounds (degrees, south positive up to north).
    lon_min, lon_max : float or None
        Longitude bounds.  When both are ``None`` no longitude constraint
        is applied (i.e. a global zonal band is returned).

    Returns
    -------
    iris.cube.Cube
        Cube with latitude and longitude collapsed, leaving only any
        remaining (cloud-type) dimensions.
    """
    lat_con = iris.Constraint(
        latitude=lambda cell: lat_min <= cell.point <= lat_max
    )
    if lon_min is not None and lon_max is not None:
        lon_con = iris.Constraint(
            longitude=lambda cell: lon_min <= cell.point <= lon_max
        )
        region = cube.extract(lat_con & lon_con)
    else:
        region = cube.extract(lat_con)
    if region is None:
        raise ValueError(
            f"No data found in region lat=[{lat_min}, {lat_max}]"
            + (f" lon=[{lon_min}, {lon_max}]" if lon_min is not None else "")
        )
    lat = region.coord("latitude")
    lon = region.coord("longitude")
    if not lat.has_bounds():
        lat.guess_bounds()
    if not lon.has_bounds():
        lon.guess_bounds()
    return region.collapsed(
        ["latitude", "longitude"],
        iris.analysis.MEAN,
        weights=iris.analysis.cartography.area_weights(region),
    )


def _plot_regional_cloud_type_averages(feedback_cube, regions, cfg):
    """Plot area-weighted regional averages for all cloud types.

    Produces a single figure with one subplot per region.  Each subplot
    shows all cloud-type bins together:

    * Two cloud-type dims (e.g. tau × pressure) → ``pcolormesh`` heatmap.
    * One cloud-type dim → bar chart.

    Parameters
    ----------
    feedback_cube : iris.cube.Cube
        Feedback cube with shape ``(cloud_type_dims..., lat, lon)``.
    regions : list of dict
        Each dict must have ``name``, ``lat_min``, ``lat_max`` and
        optionally ``lon_min``, ``lon_max``.
    cfg : dict
        ESMValTool configuration dictionary.

    Returns
    -------
    str
        Path of the saved figure.
    """
    colormap_range = cfg.get("colormap_range_histogram")

    cloud_type_coords = [
        c
        for c in feedback_cube.coords(dim_coords=True)
        if c.name() not in {"latitude", "longitude"}
    ]

    n_regions = len(regions)
    fig, axes = plt.subplots(
        1,
        n_regions,
        figsize=(4 * n_regions, 4),
        squeeze=False,
    )

    for col, region in enumerate(regions):
        ax = axes[0, col]
        name = region.get("name", f"region_{col}")
        lat_min = region["lat_min"]
        lat_max = region["lat_max"]
        lon_min = region.get("lon_min")
        lon_max = region.get("lon_max")

        try:
            regional = _extract_region_mean(
                feedback_cube, lat_min, lat_max, lon_min, lon_max
            )
        except ValueError as exc:
            logger.warning("Skipping region '%s': %s", name, exc)
            ax.set_visible(True)
            ax.set_title(f"{name}\n(no data)")
            continue

        data = np.ma.filled(regional.data, np.nan)
        if colormap_range:
            vmin, vmax = colormap_range[0], colormap_range[1]
        else:
            max_abs = float(np.nanmax(np.abs(data))) or 1.0
            vmin, vmax = -max_abs, max_abs

        if regional.ndim == 2 and len(cloud_type_coords) >= 2:
            # data is (tau, plev); transpose so tau is on x and plev on y.
            im = ax.pcolormesh(data.T, cmap="RdBu_r", vmin=vmin, vmax=vmax)
            plt.colorbar(
                im,
                ax=ax,
                label=str(regional.units),
                orientation="vertical",
                extend="both",
            )
            x_coord = cloud_type_coords[0]
            y_coord = cloud_type_coords[1]
            ax.set_xlabel(x_coord.var_name or x_coord.name())
            ax.set_ylabel(y_coord.var_name or y_coord.name())
            ax.set_xticks(np.arange(len(x_coord.points)) + 0.5)
            ax.set_yticks(np.arange(len(y_coord.points)) + 0.5)
            ax.set_xticklabels(
                [f"{v:.0f}" for v in x_coord.points], rotation=45, ha="right"
            )
            ax.set_yticklabels([f"{v:.0f}" for v in y_coord.points])
        elif regional.ndim == 1 and len(cloud_type_coords) >= 1:
            coord = cloud_type_coords[0]
            ax.bar(range(len(coord.points)), data, color="steelblue")
            ax.axhline(0, color="k", linewidth=0.8)
            ax.set_ylim(vmin, vmax)
            ax.set_ylabel(str(regional.units))
            ax.set_xticks(range(len(coord.points)))
            ax.set_xticklabels(
                [f"{v:.0f}" for v in coord.points], rotation=45, ha="right"
            )
            ax.set_xlabel(coord.var_name or coord.name())
        else:
            ax.text(
                0.5,
                0.5,
                f"{float(np.nanmean(data)):.4f} {regional.units}",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )

        ax.set_title(name)

    fig.suptitle("Regional cloud-type feedbacks (clisccp)")
    plt.tight_layout()
    filename = get_plot_filename("clisccp_regional_cloud_type_averages", cfg)
    plt.savefig(filename)
    plt.close(fig)
    return filename


def _regression_slope_against_index(field_cube, index_cube):
    """Regress field(time, ...) against index(time) and return slope cube."""
    if field_cube.shape[0] != index_cube.shape[0]:
        raise ValueError(
            "Mismatched time dimension lengths for field and index: "
            f"{field_cube.shape[0]} vs {index_cube.shape[0]}"
        )

    y_data = np.ma.asarray(field_cube.data)
    x_data = np.ma.asarray(index_cube.data)
    time_len = y_data.shape[0]
    x_broadcast = np.ma.array(
        np.broadcast_to(
            x_data.reshape((time_len,) + (1,) * (y_data.ndim - 1)),
            y_data.shape,
        )
    )

    combined_mask = np.ma.getmaskarray(y_data) | np.ma.getmaskarray(
        x_broadcast
    )
    y = np.ma.array(y_data, mask=combined_mask)
    x = np.ma.array(x_broadcast, mask=combined_mask)

    x_mean = x.mean(axis=0)
    y_mean = y.mean(axis=0)
    cov_xy = ((x - x_mean) * (y - y_mean)).mean(axis=0)
    var_x = ((x - x_mean) ** 2).mean(axis=0)
    slope = np.ma.masked_where(np.abs(var_x) < 1.0e-20, cov_xy / var_x)

    out_cube = field_cube[0].copy(data=slope)
    out_cube.var_name = "clisccp_feedback"
    out_cube.standard_name = None
    out_cube.long_name = (
        "Cloud fraction feedback by cloud type (anomaly regression)"
    )
    field_units = str(field_cube.units)
    if field_units in {"1", "no_unit", "unknown"}:
        out_cube.units = "K-1"
    else:
        out_cube.units = f"{field_units} K-1"
    return out_cube


def _plot_robinson_map(cube_2d, title, filename, vmin=None, vmax=None):
    """Plot a 2D lat-lon cube on Robinson projection."""
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    if vmin is not None and vmax is not None:
        levels = np.linspace(vmin, vmax, 21)
        contour = iplt.contourf(
            cube_2d, axes=ax, cmap="RdBu_r", levels=levels, extend="both"
        )
    else:
        contour = iplt.contourf(cube_2d, axes=ax, cmap="RdBu_r")
    ax.coastlines()
    ax.set_global()
    ax.set_title(title)
    cbar = plt.colorbar(contour, ax=ax, orientation="horizontal", pad=0.08)
    cbar.set_label(cube_2d.units)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close(fig)


def _plot_all_cloud_types(feedback_cube, cfg):
    """Create Robinson maps for each cloud-type slice if present."""
    colormap_range = cfg.get("colormap_range")
    vmin = colormap_range[0] if colormap_range else None
    vmax = colormap_range[1] if colormap_range else None

    dim_coords = feedback_cube.coords(dim_coords=True)
    extra_coords = [
        c for c in dim_coords if c.name() not in {"latitude", "longitude"}
    ]

    plot_files = []
    records = []
    if not extra_coords:
        filename = get_plot_filename("clisccp_feedback_map", cfg)
        _plot_robinson_map(
            feedback_cube, "clisccp feedback", filename, vmin=vmin, vmax=vmax
        )
        plot_files.append(filename)
        records.append("clisccp feedback")
        return plot_files, records

    dim_sizes = [len(c.points) for c in extra_coords]
    for idx_tuple in np.ndindex(*dim_sizes):
        indexers = [slice(None)] * feedback_cube.ndim
        parts = []
        for coord, idx in zip(extra_coords, idx_tuple, strict=True):
            dim_idx = feedback_cube.coord_dims(coord)[0]
            indexers[dim_idx] = idx
            parts.append(f"{coord.var_name or coord.name()}_{idx}")

        sub_cube = feedback_cube[tuple(indexers)]
        if sub_cube.ndim != 2:
            continue

        suffix = "_".join(parts)
        title = f"clisccp feedback {suffix}"
        filename = get_plot_filename(f"clisccp_feedback_map_{suffix}", cfg)
        _plot_robinson_map(sub_cube, title, filename, vmin=vmin, vmax=vmax)
        plot_files.append(filename)
        records.append(title)

    return plot_files, records


def main(config):
    """Run the diagnostic."""
    clisccp_meta = _get_input_metadata(config, "clisccp")
    tas_meta = _get_input_metadata(config, "tas")

    clisccp_cube = iris.load_cube(clisccp_meta["filename"])
    tas_cube = iris.load_cube(tas_meta["filename"])

    feedback_cube = _regression_slope_against_index(clisccp_cube, tas_cube)
    nc_file = get_diagnostic_filename("clisccp_anomaly_feedback", config)
    iris.save(feedback_cube, nc_file)

    ancestors = [clisccp_meta["filename"], tas_meta["filename"]]
    record = {
        "caption": (
            "Cloud-fraction feedback by cloud type computed from linear "
            "regression of deseasonalized clisccp anomalies against "
            "deseasonalized global-mean tas anomalies."
        ),
        "statistics": ["anomaly", "mean"],
        "domains": ["global"],
        "plot_types": ["map"],
        "authors": ["bodas-salcedo_alejandro"],
        "references": [],
        "ancestors": ancestors,
    }

    with ProvenanceLogger(config) as provenance_logger:
        provenance_logger.log(nc_file, record)

        create_plots = config.get("create_plots", True)
        if create_plots:
            plot_files, titles = _plot_all_cloud_types(feedback_cube, config)
            for plot_file, title in zip(plot_files, titles, strict=True):
                plot_record = dict(record)
                plot_record["caption"] = title
                provenance_logger.log(plot_file, plot_record)

        regions = config.get("regions", [])
        if regions:
            regional_file = _plot_regional_cloud_type_averages(
                feedback_cube, regions, config
            )
            regional_record = dict(record)
            regional_record["caption"] = (
                "Regional cloud-type average feedbacks"
            )
            provenance_logger.log(regional_file, regional_record)


if __name__ == "__main__":
    with run_diagnostic() as configuration:
        main(configuration)
