"""Relative vorticity (COSIMA Relative_Vorticity).

Translation of ``03-Advanced-Recipes/Relative_Vorticity.ipynb`` onto
CMORised data.  The vertical component of relative vorticity is

    zeta = dv/dx - du/dy

The notebook presents two methods: an ``xgcm`` calculation on MOM5's
staggered B-grid (recommended) and a naive one that converts degrees to
metres and differentiates with ``xarray`` (shown to be wrong, because
it ignores the grid staggering and the convergence of the meridians).

Neither applies unchanged here.  After the ``regrid`` preprocessor
``uo`` and ``vo`` are colocated at tracer points on a regular grid, so
the staggering the ``xgcm`` method exists to handle is already gone.
What remains of the notebook's warning is the metric term, and that is
kept: :func:`cosima_common.cell_lengths` shrinks the zonal cell length
by cos(lat), which is precisely the correction the naive method
dropped.

Vorticity is normalised by the local Coriolis parameter, so the plotted
Rossby number is comparable between latitudes.
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


def report_eddy_resolution(cube):
    """State plainly whether this grid can support mesoscale eddies.

    The COSIMA original runs on ACCESS-OM2-01, where the point of the
    figure is the eddy field.  Standard CMIP6 ocean grids are 1 degree
    and cannot resolve a deformation radius anywhere outside the deep
    tropics, so the same figure shows the shear of the resolved
    large-scale currents instead -- a real field, but not the one the
    original is about.  Saying which of the two you are looking at is
    the difference between a diagnostic and a decoration.
    """
    lat, ratio = cc.eddy_resolution_ratio(cube)
    permitting = ratio >= cc.EDDY_PERMITTING_RATIO
    resolving = ratio >= cc.EDDY_RESOLVING_RATIO

    if not permitting.any():
        logger.warning(
            "this grid resolves no deformation radius anywhere (max "
            "L_d/dx = %.2f, eddy-permitting needs %.1f). The vorticity "
            "below is the shear of the resolved mean currents, NOT a "
            "mesoscale eddy field. For eddies, use a HighResMIP member "
            "at 1/4 degree or finer.", ratio.max(),
            cc.EDDY_PERMITTING_RATIO)
    else:
        span = lat[permitting]
        logger.info(
            "eddy-permitting (L_d/dx >= %.1f) between %.0f and %.0f "
            "degrees; eddy-resolving over %.0f%% of latitudes",
            cc.EDDY_PERMITTING_RATIO, span.min(), span.max(),
            100 * resolving.mean())
    return lat, ratio


def vorticity_at_depth(uo, vo, depth=0.0):
    """Relative vorticity on the level nearest ``depth``."""
    u_cube = cc.time_mean(uo)
    v_cube = cc.time_mean(vo)
    depth_name = cc.depth_coord_name(u_cube)
    if depth_name is not None:
        levels = u_cube.coord(depth_name).points
        index = int(np.argmin(np.abs(levels - depth)))
        constraint = iris.Constraint(**{depth_name: levels[index]})
        u_cube = u_cube.extract(constraint)
        v_cube = v_cube.extract(constraint)
        logger.info("Using level %.1f m (requested %.1f m)",
                    levels[index], depth)

    u_field = cc.masked_data(u_cube)
    v_field = cc.masked_data(v_cube)
    return cc.relative_vorticity(u_field, v_field, u_cube), u_cube


def main(cfg):
    """Run the diagnostic."""
    depth = float(cfg.get("depth", 0.0))
    grouped = group_metadata(cfg["input_data"].values(), "dataset")

    for dataset, metadata in sorted(grouped.items()):
        by_variable = group_metadata(metadata, "short_name")
        if "uo" not in by_variable or "vo" not in by_variable:
            logger.warning("%s: needs both uo and vo, skipping", dataset)
            continue

        uo = cc.load_cube(by_variable["uo"][0]["filename"], "uo")
        vo = cc.load_cube(by_variable["vo"][0]["filename"], "vo")
        zeta, grid = vorticity_at_depth(uo, vo, depth=depth)
        _, ratio = report_eddy_resolution(grid)

        lat = cc.lat_1d(grid)
        fcor = cc.coriolis(lat)[:, None]
        with np.errstate(divide="ignore", invalid="ignore"):
            rossby = zeta / fcor
        rossby = np.ma.masked_where(
            np.broadcast_to(np.abs(lat)[:, None] < 2.0, rossby.shape), rossby)

        figure, axes_list = plt.subplots(1, 2, figsize=(16, 5))
        scale = np.ma.max(np.ma.abs(zeta)) or 1.0
        mesh = axes_list[0].pcolormesh(
            cc.lon_1d(grid), lat, zeta, cmap="RdBu_r", shading="auto",
            vmin=-scale, vmax=scale)
        axes_list[0].set_title(
            f"{dataset}: relative vorticity at {depth:.0f} m", fontsize=11)
        figure.colorbar(mesh, ax=axes_list[0], shrink=0.8, pad=0.03,
                        label="s$^{-1}$")

        mesh = axes_list[1].pcolormesh(
            cc.lon_1d(grid), lat, rossby, cmap="RdBu_r", shading="auto",
            vmin=-0.5, vmax=0.5)
        axes_list[1].set_title(
            f"{dataset}: Rossby number $\\zeta/f$", fontsize=11)
        figure.colorbar(mesh, ax=axes_list[1], shrink=0.8, pad=0.03,
                        label="dimensionless")
        for axes in axes_list:
            axes.set_xlabel("Longitude")
            axes.set_ylabel("Latitude")

        # The reader needs to know which field this is before reading it.
        if ratio.max() < cc.EDDY_PERMITTING_RATIO:
            figure.suptitle(
                f"{dataset}: grid resolves no deformation radius "
                f"(max $L_d/\\Delta x$ = {ratio.max():.2f}) -- this is "
                f"large-scale current shear, not mesoscale eddies",
                fontsize=10, color="0.25")
        else:
            figure.suptitle(
                f"{dataset}: eddy-permitting where "
                f"$L_d/\\Delta x \\geq$ "
                f"{cc.EDDY_PERMITTING_RATIO:.0f} "
                f"(max {ratio.max():.1f})", fontsize=10, color="0.25")

        plot_path = get_plot_filename(
            f"cosima_relative_vorticity_{dataset}", cfg)
        figure.savefig(plot_path, bbox_inches="tight", dpi=150)
        plt.close(figure)
        logger.info("Wrote %s", plot_path)

        coords = [(lat, "latitude", "degrees_north"),
                  (cc.lon_1d(grid), "longitude", "degrees_east")]
        cubes = iris.cube.CubeList([
            cc.make_cube(np.ma.filled(zeta, np.nan), coords, "rv", "s-1",
                         "ocean relative vorticity"),
            cc.make_cube(np.ma.filled(rossby, np.nan), coords, "rossby", "1",
                         "Rossby number"),
        ])
        path = get_diagnostic_filename(
            f"cosima_relative_vorticity_{dataset}", cfg)
        iris.save(cubes, path)
        record = {
            "caption": f"Relative vorticity and Rossby number of {dataset} "
                       f"at {depth:.0f} m.",
            "statistics": ["mean"],
            "domains": ["global"],
            "plot_types": ["map"],
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
