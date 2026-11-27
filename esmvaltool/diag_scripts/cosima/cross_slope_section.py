"""Cross-slope hydrographic section (COSIMA Cross-slope_section).

Translation of ``03-Advanced-Recipes/Cross-slope_section.ipynb`` onto
CMORised data.  The notebook picks a point on a stored 1000 m isobath
file, takes the direction normal to the bathymetry there, and uses
``metpy.interpolate.cross_section`` to cut a temperature section from
the shelf, across the slope, to the open ocean.

Two of those pieces are replaced:

- the stored ``Antarctic_slope_contour_1000m*.nc`` files are specific to
  the ACCESS-OM2/OM3 grids, so the contour is extracted from the
  model's own ``deptho`` with
  :func:`cosima_common.extract_contour`, which is the
  ``Create_contour.ipynb`` algorithm;
- ``metpy`` is not an ESMValTool dependency, so the section is sampled
  with the ESMValCore ``extract_transect`` preprocessor if the recipe
  asks for it, or here with a straight bilinear interpolation along the
  great-circle path between the two endpoints.

The endpoints are set the way the notebook sets them: from the chosen
contour point, ``shelf_distance`` degrees inshore along the inward
normal and ``ocean_distance`` degrees offshore along the outward
normal.

Potential density (sigma2) contours are overlaid on the temperature
section, as in the notebook.
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


def section_endpoints(deptho_cube, contour_depth=1000.0, target_lon=-60.0,
                      shelf_distance=3.0, ocean_distance=5.0):
    """Endpoints of a section normal to the bathymetry contour.

    Returns ``((lon0, lat0), (lon1, lat1))`` running from the shelf to
    the open ocean through the contour point nearest ``target_lon``.
    """
    depth = cc.masked_data(deptho_cube)
    x_contour, y_contour = cc.extract_contour(
        np.ma.filled(depth, 0.0), contour_depth)
    if x_contour.size == 0:
        raise ValueError(
            f"no {contour_depth:.0f} m contour found in deptho")

    rows, cols = cc.contour_cells(x_contour, y_contour, depth.shape)
    lat = cc.lat_1d(deptho_cube)
    lon = cc.lon_1d(deptho_cube)
    contour_lat, contour_lon = cc.contour_indices_to_coordinates(
        rows, cols, lat, lon)

    # Longitudes may be 0-360 or -180-180; compare on the circle.
    offset = np.abs(((contour_lon - target_lon + 180.0) % 360.0) - 180.0)
    index = int(np.argmin(offset))
    point_lon = float(contour_lon[index])
    point_lat = float(contour_lat[index])
    logger.info("Contour point at %.2fE %.2fN", point_lon, point_lat)

    _, _, n_x, n_y = cc.topographic_unit_vectors(depth, deptho_cube)
    normal_x = float(np.ma.filled(n_x[rows[index], cols[index]], 0.0))
    normal_y = float(np.ma.filled(n_y[rows[index], cols[index]], -1.0))
    magnitude = np.hypot(normal_x, normal_y)
    if magnitude == 0.0:
        normal_x, normal_y, magnitude = 0.0, -1.0, 1.0
    normal_x, normal_y = normal_x / magnitude, normal_y / magnitude

    start = (point_lon - shelf_distance * normal_x,
             point_lat - shelf_distance * normal_y)
    end = (point_lon + ocean_distance * normal_x,
           point_lat + ocean_distance * normal_y)
    return start, end


# sample_section now lives in cosima_common so that this script and
# notebooks/10-Cross_slope_section.ipynb cannot drift apart: the two had
# different longitude handling, and only the notebook's was correct.
sample_section = cc.sample_section


def main(cfg):
    """Run the diagnostic."""
    contour_depth = float(cfg.get("contour_depth", 1000.0))
    target_lon = float(cfg.get("section_longitude", -60.0))
    shelf_distance = float(cfg.get("shelf_distance", 3.0))
    ocean_distance = float(cfg.get("ocean_distance", 5.0))
    n_points = int(cfg.get("section_points", 60))
    grouped = group_metadata(cfg["input_data"].values(), "dataset")

    for dataset, metadata in sorted(grouped.items()):
        by_variable = group_metadata(metadata, "short_name")
        missing = [n for n in ("thetao", "deptho") if n not in by_variable]
        if missing:
            logger.warning("%s: missing %s, skipping",
                           dataset, ", ".join(missing))
            continue

        deptho = cc.load_cube(by_variable["deptho"][0]["filename"], "deptho")
        thetao = cc.time_mean(
            cc.load_cube(by_variable["thetao"][0]["filename"], "thetao"))

        try:
            start, end = section_endpoints(
                deptho, contour_depth=contour_depth, target_lon=target_lon,
                shelf_distance=shelf_distance,
                ocean_distance=ocean_distance)
        except ValueError as error:
            logger.warning("%s: %s", dataset, error)
            continue

        temperature, lons, lats = sample_section(
            thetao, start, end, n_points=n_points)
        depth_name = cc.depth_coord_name(thetao)
        depth = thetao.coord(depth_name).points

        sigma_section = None
        if "so" in by_variable:
            so = cc.time_mean(
                cc.load_cube(by_variable["so"][0]["filename"], "so"))
            sigma2 = cc.potential_density_anomaly(thetao, so, p_ref=2000)
            sigma_cube = thetao.copy(np.ma.filled(sigma2, np.nan))
            sigma_section, _, _ = sample_section(
                sigma_cube, start, end, n_points=n_points)

        figure, axes = plt.subplots(figsize=(10, 6))
        mesh = axes.pcolormesh(lats, depth, temperature, cmap="RdYlBu_r",
                               shading="auto")
        if sigma_section is not None:
            contours = axes.contour(lats, depth, sigma_section, colors="k",
                                    linewidths=0.7, levels=10)
            axes.clabel(contours, inline=True, fontsize=7, fmt="%.2f")
        axes.invert_yaxis()
        axes.set_xlabel("Latitude along section")
        axes.set_ylabel("Depth (m)")
        axes.set_title(
            f"{dataset}: cross-slope section through the "
            f"{contour_depth:.0f} m isobath\n"
            f"{start[0]:.1f}E {start[1]:.1f}N to "
            f"{end[0]:.1f}E {end[1]:.1f}N", fontsize=11)
        figure.colorbar(mesh, ax=axes, shrink=0.85, pad=0.03,
                        label="Potential temperature ($^\\circ$C)")

        plot_path = get_plot_filename(
            f"cosima_cross_slope_section_{dataset}", cfg)
        figure.savefig(plot_path, bbox_inches="tight", dpi=150)
        plt.close(figure)
        logger.info("Wrote %s", plot_path)

        out = cc.make_cube(
            np.ma.filled(temperature, np.nan),
            [(depth, "depth", "m"), (lats, "latitude", "degrees_north")],
            "thetao", "degC", "potential temperature on a cross-slope "
                              "section")
        path = get_diagnostic_filename(
            f"cosima_cross_slope_section_{dataset}", cfg)
        iris.save(out, path)
        record = {
            "caption": f"Cross-slope hydrographic section of {dataset} "
                       f"normal to the {contour_depth:.0f} m isobath.",
            "statistics": ["mean"],
            "domains": ["shpolar"],
            "plot_types": ["vert"],
            "authors": ["beucher_romain"],
            "references": ["cosima_recipes"],
            "ancestors": [m["filename"] for m in metadata],
        }
        with ProvenanceLogger(cfg) as logger_:
            logger_.log(path, record)
            logger_.log(plot_path, record)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
