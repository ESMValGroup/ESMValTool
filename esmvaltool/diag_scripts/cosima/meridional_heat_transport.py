"""Ocean meridional heat transport (COSIMA Meridional_heat_transport).

Translation of ``03-Advanced-Recipes/Meridional_heat_transport.ipynb``
onto CMORised CMIP6 data.  The notebook offers two methods and this
diagnostic keeps both:

Method 1 -- online diagnostic.  MOM5 integrates the resolved advective
heat flux online into ``temp_yflux_adv_int_z``.  The CMIP6 equivalent is
``hfbasin``, the basin-integrated northward ocean heat transport, which
is already zonally integrated, so the diagnostic only converts W to PW
and time-averages.  Fewer models publish ``hfbasin`` than ``hfds``.

Method 2 -- surface flux integration.  Where ``hfbasin`` is missing, the
transport is inferred by integrating the net surface heat flux
``hfds`` (plus the frazil term ``hfsifrazil`` where available) from the
south pole northwards.  This assumes a steady state, exactly as the
notebook warns, and the result is centred so that the implied transport
vanishes at both ends.

Both methods are plotted together when both are available.
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

WATTS_TO_PW = 1.0e-15


def _basin_name(cube):
    """Readable basin label from an ``hfbasin`` region coordinate."""
    for name in ("region", "basin", "sector"):
        if cube.coords(name):
            points = cube.coord(name).points
            if points.size:
                value = points[0]
                return value.decode() if isinstance(value, bytes) else str(
                    value)
    return "global_ocean"


def transport_from_hfbasin(cube):
    """Method 1: read the online transport straight off ``hfbasin``."""
    cube = cc.time_mean(cube)
    data = cc.masked_data(cube) * WATTS_TO_PW
    return cc.lat_1d(cube), np.squeeze(data)


def transport_from_surface_flux(flux_cube, area, n_bins=None):
    """Method 2: integrate the surface heat flux from the south pole.

    ``MHF(y) = int_{-90}^{y} (Q * dA)``, accumulated over latitude bins
    that span the whole grid, then centred by ``(MHF[0] - MHF[-1]) / 2``
    as the notebook does.

    The bins come from :func:`cosima_common.latitude_bin_edges`, not
    from the cube's latitude coordinate.  On a tripolar grid the
    latter is 2D and any single column of it can stop at the fold, so
    edges taken from it drop the Arctic without saying so.
    """
    flux_cube = cc.time_mean(flux_cube)
    flux = cc.masked_data(flux_cube)
    lat_field = cc.lat_2d(flux_cube)

    edges = cc.latitude_bin_edges(flux_cube, n_bins=n_bins)
    total_flux = flux * area

    # Nothing may fall outside the bins: a discarded cell is heat that
    # silently never enters the integral.
    binned_cells = ((lat_field >= edges[0]) & (lat_field <= edges[-1]))
    covered = float(np.ma.mean(
        np.ma.masked_where(np.ma.getmaskarray(flux[..., :, :]),
                           binned_cells).astype(float)))
    if covered < 0.999:
        logger.warning(
            "only %.1f%% of ocean cells fall within the latitude bins; the "
            "remaining flux is not in the transport", 100 * covered)

    binned = cc.latitude_bin_sum(total_flux, lat_field, edges)
    transport = np.cumsum(binned)
    transport = transport + (transport[0] - transport[-1]) / 2.0
    return cc.bin_centres(edges), transport * WATTS_TO_PW


def _save(cube_data, lat, var_name, long_name, basename, cfg, ancestors,
          caption):
    """Write a 1D latitude profile and log its provenance."""
    cube = cc.make_cube(
        np.ma.filled(cube_data, np.nan), [(lat, "latitude", "degrees_north")],
        var_name, "PW", long_name)
    path = get_diagnostic_filename(basename, cfg)
    iris.save(cube, path)
    record = {
        "caption": caption,
        "statistics": ["mean"],
        "domains": ["global"],
        "plot_types": ["line"],
        "authors": ["beucher_romain"],
        "references": ["cosima_recipes"],
        "ancestors": ancestors,
    }
    with ProvenanceLogger(cfg) as logger_:
        logger_.log(path, record)
    return path


def main(cfg):
    """Run the diagnostic."""
    grouped = group_metadata(cfg["input_data"].values(), "dataset")
    figure, axes = plt.subplots(figsize=(9, 6))
    plotted = 0

    for index, (dataset, metadata) in enumerate(sorted(grouped.items())):
        by_variable = group_metadata(metadata, "short_name")
        style = cc.style_for(dataset, index)
        ancestors = [m["filename"] for m in metadata]

        if "hfbasin" in by_variable:
            meta = by_variable["hfbasin"][0]
            cube = cc.load_cube(meta["filename"], "hfbasin")
            lat, transport = transport_from_hfbasin(cube)
            if transport.ndim > 1:
                # Several basins; the global one is what the notebook plots.
                logger.info("%s: hfbasin has %d basins, using the first",
                            dataset, transport.shape[0])
                transport = transport[0]
            axes.plot(lat, transport, color=style["color"], linestyle="-",
                      linewidth=style["thick"],
                      label=f"{dataset} (hfbasin)")
            _save(transport, lat, "hfbasin",
                  f"northward ocean heat transport ({_basin_name(cube)})",
                  f"cosima_mht_hfbasin_{dataset}", cfg, ancestors,
                  f"Meridional heat transport of {dataset} from the online "
                  f"hfbasin diagnostic.")
            plotted += 1

        if "hfds" in by_variable:
            meta = by_variable["hfds"][0]
            cube = cc.load_cube(meta["filename"], "hfds")
            # This diagnostic stays on the model's own grid, which for
            # almost every CMIP6 ocean is tripolar, so the areas cannot
            # come from a representative 1D latitude. cell_area falls
            # back through areacello -> the 2D corner bounds (which
            # CMOR mandates) -> 1D coordinates, so a missing areacello
            # is no longer fatal.
            try:
                area = cc.cell_area(cube)
            except ValueError as error:
                logger.error("%s: %s", dataset, error)
                continue
            if "hfsifrazil" in by_variable:
                frazil = cc.load_cube(
                    by_variable["hfsifrazil"][0]["filename"], "hfsifrazil")
                cube = cc.time_mean(cube).copy(
                    cc.masked_data(cc.time_mean(cube))
                    + cc.masked_data(cc.time_mean(frazil)))
                logger.info("%s: added hfsifrazil to hfds", dataset)
            lat, transport = transport_from_surface_flux(cube, area)
            axes.plot(lat, transport, color=style["color"], linestyle="--",
                      linewidth=style["thick"],
                      label=f"{dataset} (surface flux)")
            _save(transport, lat, "hfds_mht",
                  "northward ocean heat transport inferred from hfds",
                  f"cosima_mht_surfaceflux_{dataset}", cfg, ancestors,
                  f"Meridional heat transport of {dataset} inferred by "
                  f"integrating the net surface heat flux.")
            plotted += 1

    if not plotted:
        logger.warning("Neither hfbasin nor hfds was available")
        return

    axes.axhline(0.0, color="black", linewidth=1)
    axes.set_xlim(-90, 90)
    axes.set_xticks(np.arange(-90, 91, 30))
    axes.set_xlabel("Latitude")
    axes.set_ylabel("Heat transport (PW)")
    axes.set_title("Global ocean meridional heat transport")
    axes.legend(frameon=False, fontsize=9)
    axes.grid(alpha=0.3)

    plot_path = get_plot_filename("cosima_meridional_heat_transport", cfg)
    figure.savefig(plot_path, bbox_inches="tight", dpi=150)
    plt.close(figure)
    with ProvenanceLogger(cfg) as logger_:
        logger_.log(plot_path, {
            "caption": "Model meridional heat transports from archived hfbasin and surface fluxes.",
            "statistics": ["mean"],
            "domains": ["global"],
            "plot_types": ["line"],
            "authors": ["beucher_romain"],
            "references": ["cosima_recipes"],
            "ancestors": [m["filename"] for m in cfg["input_data"].values()],
        })
    logger.info("Wrote %s", plot_path)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
