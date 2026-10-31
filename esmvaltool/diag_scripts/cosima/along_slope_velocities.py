"""Antarctic slope current (COSIMA Along-slope-velocities).

Translation of ``03-Advanced-Recipes/Along-slope-velocities.ipynb`` onto
CMORised data.  The horizontal velocity is projected onto the local
bathymetric contour to separate the along-slope flow, which is the
Antarctic Slope Current, from the cross-slope flow that exchanges water
with the shelf:

    u_along = u * h_y/|grad h| - v * h_x/|grad h|
    v_cross = u * h_x/|grad h| + v * h_y/|grad h|

with ``h`` the ocean depth ``deptho``.  The notebook takes the
bathymetric gradient with ``xgcm`` on MOM5's B-grid; on the regular
grid that ``regrid`` produces, the gradient is a central difference
with spherical metrics, and ``uo``/``vo`` are already colocated with
``deptho``.

The along-slope component is also averaged into depth bins to give the
classic slope-current profile against bathymetry.
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

DEFAULT_ISOBATH_BINS = np.array(
    [0.0, 500.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 6000.0])


def slope_components(uo, vo, deptho):
    """Along- and cross-slope surface velocity components."""
    u_cube = cc.surface_field(cc.time_mean(uo))
    v_cube = cc.surface_field(cc.time_mean(vo))
    u_field = cc.masked_data(u_cube)
    v_field = cc.masked_data(v_cube)
    depth = cc.masked_data(deptho)
    along, cross = cc.project_along_slope(u_field, v_field, depth, u_cube)
    return along, cross, depth, u_cube


def isobath_profile(field, depth, bins, area):
    """Area-weighted mean of ``field`` in each bathymetry bin."""
    weighted = cc.bin_sum(field * area, depth, bins)
    weights = cc.bin_sum(
        np.ma.masked_where(np.ma.getmaskarray(field), area), depth, bins)
    with np.errstate(divide="ignore", invalid="ignore"):
        profile = np.where(weights > 0, weighted / weights, np.nan)
    return np.ma.masked_invalid(profile)


def main(cfg):
    """Run the diagnostic."""
    bins = np.asarray(cfg.get("isobath_bins", DEFAULT_ISOBATH_BINS),
                      dtype=float)
    max_lat = float(cfg.get("max_latitude", -59.0))
    grouped = group_metadata(cfg["input_data"].values(), "dataset")

    figure, axes_list = plt.subplots(1, 2, figsize=(15, 5.5))
    plotted = 0

    for index, (dataset, metadata) in enumerate(sorted(grouped.items())):
        by_variable = group_metadata(metadata, "short_name")
        missing = [n for n in ("uo", "vo", "deptho") if n not in by_variable]
        if missing:
            logger.warning("%s: missing %s, skipping",
                           dataset, ", ".join(missing))
            continue

        uo = cc.load_cube(by_variable["uo"][0]["filename"], "uo")
        vo = cc.load_cube(by_variable["vo"][0]["filename"], "vo")
        deptho = cc.load_cube(by_variable["deptho"][0]["filename"], "deptho")
        along, cross, depth, grid = slope_components(uo, vo, deptho)

        area = cc.cell_area(grid)
        lat = cc.lat_2d(grid)
        southern = lat < max_lat
        along_south = np.ma.masked_where(~southern, along)
        depth_south = np.ma.masked_where(~southern, depth)

        profile = isobath_profile(along_south, depth_south, bins, area)
        centres = cc.bin_centres(bins)
        style = cc.style_for(dataset, index)
        axes_list[1].plot(profile, centres, color=style["color"],
                          linestyle=style["dash"], linewidth=2.0,
                          marker="o", label=dataset)
        plotted += 1

        if plotted == 1:
            scale = np.ma.max(np.ma.abs(along_south)) or 0.2
            mesh = axes_list[0].pcolormesh(
                cc.lon_1d(grid), cc.lat_1d(grid), along_south,
                cmap="RdBu_r", shading="auto", vmin=-scale, vmax=scale)
            axes_list[0].set_title(
                f"{dataset}: along-slope surface velocity", fontsize=11)
            axes_list[0].set_xlabel("Longitude")
            axes_list[0].set_ylabel("Latitude")
            figure.colorbar(mesh, ax=axes_list[0], shrink=0.8, pad=0.03,
                            label="m s$^{-1}$")

        coords = [(cc.lat_1d(grid), "latitude", "degrees_north"),
                  (cc.lon_1d(grid), "longitude", "degrees_east")]
        cubes = iris.cube.CubeList([
            cc.make_cube(np.ma.filled(along, np.nan), coords, "u_along",
                         "m s-1", "along-slope velocity component"),
            cc.make_cube(np.ma.filled(cross, np.nan), coords, "v_cross",
                         "m s-1", "cross-slope velocity component"),
        ])
        path = get_diagnostic_filename(
            f"cosima_along_slope_{dataset}", cfg)
        iris.save(cubes, path)
        record = {
            "caption": f"Along- and cross-slope surface velocity components "
                       f"of {dataset}, projected onto the bathymetry.",
            "statistics": ["mean"],
            "domains": ["shpolar"],
            "plot_types": ["map"],
            "authors": ["beucher_romain"],
            "references": ["cosima_recipes"],
            "ancestors": [m["filename"] for m in metadata],
        }
        with ProvenanceLogger(cfg) as logger_:
            logger_.log(path, record)

    if not plotted:
        logger.warning("No dataset had uo, vo and deptho")
        plt.close(figure)
        return

    axes_list[1].axvline(0.0, color="black", linewidth=0.6)
    axes_list[1].invert_yaxis()
    axes_list[1].set_xlabel("Along-slope velocity (m s$^{-1}$)")
    axes_list[1].set_ylabel("Bathymetry (m)")
    axes_list[1].set_title(
        f"Along-slope flow by isobath, south of {max_lat:.0f}$^\\circ$",
        fontsize=11)
    axes_list[1].legend(fontsize=9, frameon=False)
    axes_list[1].grid(alpha=0.3)

    plot_path = get_plot_filename("cosima_along_slope_velocities", cfg)
    figure.savefig(plot_path, bbox_inches="tight", dpi=150)
    plt.close(figure)
    with ProvenanceLogger(cfg) as logger_:
        logger_.log(plot_path, {
            **record,
            "caption": "Along- and cross-slope velocities by model.",
            "ancestors": [m["filename"] for m in cfg["input_data"].values()],
        })
    logger.info("Wrote %s", plot_path)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
