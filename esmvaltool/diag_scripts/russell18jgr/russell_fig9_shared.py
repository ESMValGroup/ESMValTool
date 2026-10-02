"""Shared machinery for russell18jgr figures 9a, 9b and 9c.

Scatter plots of westerly wind band width, integrated heat uptake and
integrated carbon uptake south of 30S, with a line of best fit.
Port of the common parts of russell18jgr-fig9[abc].ncl.
"""

import logging
import os
import sys

import iris
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import russell_common as rc

logger = logging.getLogger(__name__)

CARBON_FACTOR = 0.000031536  # kg s-1 -> Pg yr-1
HEAT_FACTOR = 1.0e-15  # W -> PW


def band_width(tauu_meta):
    """Westerly band width from the zonal mean wind stress."""
    cube = rc.time_mean(
        rc.load_cube(tauu_meta["filename"], tauu_meta["short_name"])
    )
    zonal = rc.zonal_mean(cube.data)
    lat = rc.coord_1d_or_2d(cube, "latitude")
    if lat.ndim == 2:
        lat = np.mean(lat, axis=1)
    return rc.westerly_band_width(zonal, lat)


def area_for(meta, area_data, data_shape, cube):
    """Cell areas: areacello fx file, or the manual regular-grid calc."""
    area_meta = rc.match_by_dataset(area_data, meta["dataset"])
    if rc.has_regular_grid(cube):
        # some models provide the flux on a different grid than
        # areacello, so the area is computed manually (as in NCL)
        return rc.ncl_cell_area(rc.lat_1d(cube), rc.lon_1d(cube), data_shape)
    if area_meta is None:
        raise ValueError(
            f"areacello file of {meta['dataset']} not found in the "
            "recipe. If available, please copy the dataset name to the "
            "additional dataset section of areacello."
        )
    return np.ma.masked_invalid(
        rc.load_cube(area_meta["filename"], "areacello").data
    )


def integrated_flux(meta, area_data, factor):
    """Integrated flux (sum flux*area*factor) south of 30S."""
    cube = rc.time_mean(rc.load_cube(meta["filename"], meta["short_name"]))
    data = np.ma.masked_invalid(cube.data)
    area = area_for(meta, area_data, data.shape, cube)
    lat = rc.lat_1d(cube)
    return rc.southern_ocean_flux_sum(data, area, lat, factor)


def heat_flux(meta, area_data):
    """Integrated heat flux south of 30S in PW (always uses areacello)."""
    cube = rc.time_mean(rc.load_cube(meta["filename"], meta["short_name"]))
    data = np.ma.masked_invalid(cube.data)
    area_meta = rc.match_by_dataset(area_data, meta["dataset"])
    if area_meta is None:
        raise ValueError(
            f"areacello file of {meta['dataset']} not found in the "
            "recipe. If available, please copy the dataset name to the "
            "additional dataset section of areacello."
        )
    area = np.ma.masked_invalid(
        rc.load_cube(area_meta["filename"], "areacello").data
    )
    lat = rc.lat_1d(cube)
    return rc.southern_ocean_flux_sum(data, area, lat, HEAT_FACTOR)


def scatter_with_regression(
    cfg,
    xvals,
    yvals,
    datasets,
    xlabel,
    ylabel,
    title,
    xlim,
    xtick_spacing,
    plot_name,
):
    """Scatter plot with line of best fit and per-model markers."""
    from esmvaltool.diag_scripts.shared import get_plot_filename

    xvals = np.asarray(xvals, dtype=float)
    yvals = np.asarray(yvals, dtype=float)

    fig, axes = plt.subplots(figsize=(9, 6))
    slope, intercept = np.polyfit(xvals, yvals, 1)
    xline = np.linspace(xlim[0], xlim[1], len(xvals))
    yline = intercept + slope * xline
    axes.plot(xline, yline, color="black", linewidth=1)

    for x, y, dataset in zip(xvals, yvals, datasets, strict=True):
        style = rc.style_for(dataset, cfg.get("styleset", "CMIP5"))
        axes.plot(
            x,
            y,
            linestyle="none",
            marker=style["mark"],
            color=style["color"],
            markersize=8,
            label=dataset,
        )

    axes.set_xlim(*xlim)
    if xtick_spacing is not None:
        axes.set_xticks(
            np.arange(xlim[0], xlim[1] + xtick_spacing / 2, xtick_spacing)
        )
    ymin = np.floor((yvals.min() - 0.1) * 10) / 10
    ymax = np.ceil((yvals.max() + 0.1) * 10) / 10
    axes.set_ylim(ymin, ymax)
    axes.grid(color="grey", linewidth=0.5)
    axes.set_xlabel(xlabel, fontsize=9)
    axes.set_ylabel(ylabel, fontsize=9)
    axes.set_title(title, fontsize=11)
    axes.tick_params(labelsize=8)
    axes.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=7,
        frameon=False,
    )
    fig.tight_layout()

    plot_file = get_plot_filename(plot_name, cfg)
    fig.savefig(plot_file, bbox_inches="tight", dpi=200)
    plt.close(fig)
    logger.info("Wrote %s", plot_file)
    return plot_file, (xline, yline)


def save_pairs(
    cfg,
    basename,
    var_name,
    description,
    datasets,
    metas,
    pair_values,
    regline,
    plot_file,
    caption,
    ancestor_lists,
):
    """Per-dataset netCDF output + provenance, as the NCL scripts did."""
    from esmvaltool.diag_scripts.shared import (
        ProvenanceLogger,
        get_diagnostic_filename,
    )

    xline, yline = regline
    with ProvenanceLogger(cfg) as prov:
        plot_logged = False
        for i, dataset in enumerate(datasets):
            meta = metas[i]
            out = iris.cube.Cube(
                np.array(pair_values[i]),
                var_name=var_name,
                long_name=f"{description} for dataset {dataset}",
                attributes={
                    "regline_y_coord": yline.astype(float),
                    "regline_x_coord": xline.astype(float),
                },
            )
            nc_name = get_diagnostic_filename(
                f"{basename}_{dataset}_{meta['start_year']}-"
                f"{meta['end_year']}",
                cfg,
            )
            iris.save(out, nc_name)
            record = rc.provenance_record(
                caption, ancestor_lists[i], plot_types=("scatter",)
            )
            prov.log(nc_name, record)
            if not plot_logged:
                prov.log(plot_file, record)
                plot_logged = True
