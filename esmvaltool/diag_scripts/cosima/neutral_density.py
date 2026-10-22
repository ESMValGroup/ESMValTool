"""Neutral density (COSIMA Neutral_density).

Translation of ``03-Advanced-Recipes/Neutral_density.ipynb`` onto
CMORised data.

The notebook computes neutral density gamma_n with ``pygamma_n``, a
wrapper around the Jackett & McDougall (1997) algorithm.  ``pygamma_n``
is not on conda-forge and is not an ESMValTool dependency, so this
diagnostic uses it when it happens to be importable and otherwise falls
back to :func:`cosima_common.patched_potential_density`: potential
density referenced to 0, 1000, 2000, 3000 and 4000 dbar, with each
patch offset to make the composite continuous, after Reid (1994).

The fallback is an approximation and the diagnostic says so in the
plot title, the netCDF ``long_name`` and the provenance caption.  It
follows neutral surfaces well within each pressure patch and has a
discontinuous slope at the joins, so it is fine for identifying water
masses and wrong for quantitative neutral-surface work.

The notebook's own warning carries over regardless of the method:
gamma_n is referenced to a 1980s climatology and should not be applied
to future projections.
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

FUTURE_EXPERIMENTS = ("ssp", "rcp", "abrupt", "1pct")


def _warn_about_future(metadata, dataset):
    """gamma_n is tied to a 1980s climatology; say so for future runs."""
    for meta in metadata:
        exp = str(meta.get("exp", "")).lower()
        if any(exp.startswith(prefix) for prefix in FUTURE_EXPERIMENTS):
            logger.warning(
                "%s: experiment %s is not a present-day climate; neutral "
                "density is referenced to a 1980s climatology and should "
                "not be used here", dataset, meta["exp"])
            return


def main(cfg):
    """Run the diagnostic."""
    grouped = group_metadata(cfg["input_data"].values(), "dataset")
    section_lon = cfg.get("section_longitude")

    for dataset, metadata in sorted(grouped.items()):
        by_variable = group_metadata(metadata, "short_name")
        if "thetao" not in by_variable or "so" not in by_variable:
            logger.warning("%s: needs both thetao and so, skipping", dataset)
            continue
        _warn_about_future(metadata, dataset)

        thetao = cc.time_mean(
            cc.load_cube(by_variable["thetao"][0]["filename"], "thetao"))
        so = cc.time_mean(
            cc.load_cube(by_variable["so"][0]["filename"], "so"))
        gamma, method = cc.neutral_density(thetao, so)
        logger.info("%s: computed %s", dataset, method)

        depth_name = cc.depth_coord_name(thetao)
        depth = thetao.coord(depth_name).points
        lat = cc.lat_1d(thetao)
        lon = cc.lon_1d(thetao)

        if section_lon is None:
            section = np.ma.mean(gamma, axis=-1)
            section_label = "zonal mean"
        else:
            index = int(np.argmin(np.abs(lon - float(section_lon))))
            section = gamma[:, :, index]
            section_label = f"{lon[index]:.1f}$^\\circ$E"

        figure, axes = plt.subplots(figsize=(10, 6))
        mesh = axes.pcolormesh(lat, depth, section, cmap="viridis",
                               shading="auto")
        contours = axes.contour(lat, depth, section, colors="k",
                                linewidths=0.6, levels=14)
        axes.clabel(contours, inline=True, fontsize=7, fmt="%.2f")
        axes.invert_yaxis()
        axes.set_xlabel("Latitude")
        axes.set_ylabel("Depth (m)")
        axes.set_title(f"{dataset}: {method}\nsection at {section_label}",
                       fontsize=11)
        figure.colorbar(mesh, ax=axes, shrink=0.85, pad=0.03,
                        label="kg m$^{-3}$")

        plot_path = get_plot_filename(f"cosima_neutral_density_{dataset}",
                                      cfg)
        figure.savefig(plot_path, bbox_inches="tight", dpi=150)
        plt.close(figure)
        logger.info("Wrote %s", plot_path)

        long_name = (f"neutral density: {method}"
                     if method.startswith("gamma_n") else method)
        out = cc.make_cube(
            np.ma.filled(section, np.nan),
            [(depth, "depth", "m"), (lat, "latitude", "degrees_north")],
            "gamma_n", "kg m-3", long_name)
        path = get_diagnostic_filename(f"cosima_neutral_density_{dataset}",
                                       cfg)
        iris.save(out, path)
        record = {
            "caption": f"Neutral density section of {dataset}, computed as "
                       f"{method}.",
            "statistics": ["mean"],
            "domains": ["global"],
            "plot_types": ["zonal"],
            "authors": ["beucher_romain"],
            "references": ["cosima_recipes", "jackett97jpo"],
            "ancestors": [m["filename"] for m in metadata],
        }
        with ProvenanceLogger(cfg) as logger_:
            logger_.log(path, record)
            logger_.log(plot_path, record)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
