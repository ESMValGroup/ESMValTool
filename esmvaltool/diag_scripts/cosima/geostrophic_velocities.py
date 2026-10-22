"""Surface geostrophic velocities and EKE (COSIMA
Geostrophic_Velocities_from_Sea_Level).

Translation of ``03-Advanced-Recipes/Geostrophic_Velocities_from_Sea_Level.ipynb``
onto CMORised daily sea level:

    u_g = -(g/f) deta/dy       v_g =  (g/f) deta/dx

The notebook does the differencing with ``xgcm`` on MOM5's staggered
B-grid.  CMORised ``zos`` arrives on a regular latitude/longitude grid
after ``regrid``, so the derivatives here are ordinary central
differences with spherical metrics, evaluated at the tracer points
rather than at cell corners.  No interpolation between grid positions
is therefore needed, and the result sits on the same grid as ``zos``.

Because geostrophic balance fails as ``f`` approaches zero at the
equator, the diagnostic masks a band either side of it (5 degrees by
default).  The mask is applied to ``f`` before the division.

The notebook compares against the model's own surface velocities;
this diagnostic does the same when ``uo``/``vo`` are in the recipe.

It also computes the variance of the geostrophic velocity about its
own time mean, and this needs a caveat that the notebook does not have
to make.  CMIP6 publishes **no daily** ``zos``: the ``Oday`` table
holds only ``chlos``, ``omldamax``, ``phycos``, ``sos``, ``sossq``,
``tos`` and ``tossq``, and ``zos`` exists only in ``Omon``.  Monthly
sea level cannot resolve mesoscale eddies, so on monthly input this
quantity is *low-frequency sea-level-driven kinetic energy*, not eddy
kinetic energy, and the diagnostic labels it that way.  It only calls
the result EKE when the input really is sub-weekly, which in practice
means non-CMIP6 input such as native ACCESS-OM2 output.

That also settles what can be done with
``Eddy-Mean_Kinetic_Energy_Decomposition.ipynb``: it needs daily 3D
velocities, CMIP6 has neither those nor daily sea level, so no form of
mesoscale EKE is available from CMIP6 and the notebook has no faithful
translation.
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

#: Sampling interval below which the velocity variance may honestly be
#: called eddy kinetic energy.
EDDY_RESOLVING_DAYS = 7.0


def variance_labels(zos_cube):
    """Name the velocity variance according to the input frequency.

    Returns ``(var_name, long_name, plot_label)``.
    """
    step = cc.time_step_days(zos_cube)
    if step is not None and step <= EDDY_RESOLVING_DAYS:
        return ("eke", "surface eddy kinetic energy",
                "surface EKE")
    logger.warning(
        "sea level is sampled every %.0f days; the velocity variance is "
        "low-frequency sea-level-driven kinetic energy, not mesoscale EKE",
        step if step is not None else float("nan"))
    return ("ske_lowfreq",
            "low-frequency surface kinetic energy from sea level "
            "(not mesoscale EKE: input is not eddy resolving)",
            "low-frequency surface KE\n(monthly sea level: not eddy EKE)")


def geostrophic_series(zos_cube, min_lat=5.0):
    """Geostrophic velocity time series from a ``(time, lat, lon)`` cube."""
    data = cc.masked_data(zos_cube)
    if data.ndim == 2:
        data = data[None, ...]
    horizontal = next(zos_cube.slices(["latitude", "longitude"]))

    u_series, v_series = [], []
    for step in range(data.shape[0]):
        u_field, v_field = cc.geostrophic_velocities(
            data[step], horizontal, min_lat=min_lat)
        u_series.append(u_field)
        v_series.append(v_field)
    return np.ma.stack(u_series), np.ma.stack(v_series), horizontal


def _plot_speed(figure, position, speed, cube, title):
    """Filled speed map with the coastline."""
    import cartopy.crs as ccrs

    axes = figure.add_subplot(*position, projection=ccrs.PlateCarree())
    mesh = axes.pcolormesh(cc.lon_2d(cube), cc.lat_2d(cube), speed,
                           cmap="viridis", transform=ccrs.PlateCarree(),
                           shading="auto")
    axes.coastlines(resolution="110m", linewidth=0.5)
    axes.set_title(title, fontsize=11)
    figure.colorbar(mesh, ax=axes, shrink=0.7, pad=0.03,
                    label="m s$^{-1}$")
    return axes


def main(cfg):
    """Run the diagnostic."""
    min_lat = float(cfg.get("equator_mask_degrees", 5.0))
    grouped = group_metadata(cfg["input_data"].values(), "dataset")

    for dataset, metadata in sorted(grouped.items()):
        by_variable = group_metadata(metadata, "short_name")
        if "zos" not in by_variable:
            logger.warning("%s: no zos, skipping", dataset)
            continue

        zos = cc.load_cube(by_variable["zos"][0]["filename"], "zos")
        u_series, v_series, grid = geostrophic_series(zos, min_lat=min_lat)

        u_mean = np.ma.mean(u_series, axis=0)
        v_mean = np.ma.mean(v_series, axis=0)
        speed = np.ma.sqrt(u_mean**2 + v_mean**2)
        variance = cc.eddy_kinetic_energy(u_series, v_series)
        var_name, long_name, plot_label = variance_labels(zos)

        panels = [(speed, f"{dataset}: mean geostrophic speed")]
        if "uo" in by_variable and "vo" in by_variable:
            u_model = cc.surface_field(
                cc.time_mean(cc.load_cube(
                    by_variable["uo"][0]["filename"], "uo")))
            v_model = cc.surface_field(
                cc.time_mean(cc.load_cube(
                    by_variable["vo"][0]["filename"], "vo")))
            model_speed = np.ma.sqrt(
                cc.masked_data(u_model)**2
                + cc.masked_data(v_model)**2)
            panels.append(
                (model_speed,
                 f"{dataset}: modelled surface speed\n"
                 f"(geostrophic + ageostrophic)"))

        figure = plt.figure(figsize=(8 * len(panels) + 2, 5))
        for index, (field, title) in enumerate(panels):
            _plot_speed(figure, (1, len(panels) + 1, index + 1), field,
                        grid, title)
        axes = figure.add_subplot(1, len(panels) + 1, len(panels) + 1)
        mesh = axes.pcolormesh(
            cc.lon_1d(grid), cc.lat_1d(grid),
            np.ma.log10(np.ma.masked_less_equal(variance, 0)),
            cmap="magma", shading="auto")
        axes.set_title(f"{dataset}: {plot_label}", fontsize=11)
        axes.set_xlabel("Longitude")
        axes.set_ylabel("Latitude")
        figure.colorbar(mesh, ax=axes, shrink=0.7, pad=0.03,
                        label="log$_{10}$ (m$^2$ s$^{-2}$)")

        plot_path = get_plot_filename(
            f"cosima_geostrophic_velocities_{dataset}", cfg)
        figure.savefig(plot_path, bbox_inches="tight", dpi=150)
        plt.close(figure)
        logger.info("Wrote %s", plot_path)

        lat = cc.lat_1d(grid)
        lon = cc.lon_1d(grid)
        coords = [(lat, "latitude", "degrees_north"),
                  (lon, "longitude", "degrees_east")]
        cubes = iris.cube.CubeList([
            cc.make_cube(np.ma.filled(u_mean, np.nan), coords,
                         "ug", "m s-1", "surface geostrophic eastward "
                         "velocity"),
            cc.make_cube(np.ma.filled(v_mean, np.nan), coords,
                         "vg", "m s-1", "surface geostrophic northward "
                         "velocity"),
            cc.make_cube(np.ma.filled(variance, np.nan), coords,
                         var_name, "m2 s-2", long_name),
        ])
        path = get_diagnostic_filename(
            f"cosima_geostrophic_velocities_{dataset}", cfg)
        iris.save(cubes, path)
        record = {
            "caption": f"Surface geostrophic velocities of {dataset} derived "
                       f"from sea level, with the {long_name}.",
            "statistics": ["mean", "var"],
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
