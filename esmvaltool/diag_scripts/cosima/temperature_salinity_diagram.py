"""Volume-weighted T-S diagram (COSIMA Temperature_Salinity_Diagram).

Translation of ``03-Advanced-Recipes/Temperature_Salinity_Diagram.ipynb``
onto CMORised data.  The notebook weights the joint temperature-salinity
histogram by ``volcello`` using ``xhistogram``; here the weighting is
done with :func:`cosima_common.volume_weighted_ts_histogram`, which is
``numpy.histogram2d`` with volume weights and therefore conserves the
total volume exactly.

Axes are TEOS-10 quantities, as in the notebook: practical salinity
``so`` is converted to absolute salinity and potential temperature
``thetao`` to conservative temperature, because mixing CMOR's potential
temperature with absolute salinity on a T-S plot misplaces the water
masses relative to the density contours drawn over them.

Potential density (sigma2) contours are overlaid so the classical
Southern Ocean water masses can be read off directly.
"""

import logging
import os

import dask
import iris
import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors  # noqa: E402
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

def _cell_volume(thetao, volcello_cube=None):
    """Cell volume in m3, from ``volcello`` or from area * thickness."""
    if volcello_cube is not None:
        return cc.masked_data(volcello_cube)

    depth_name = cc.depth_coord_name(thetao)
    if depth_name is None:
        raise ValueError("cannot build a volume without a depth coordinate")
    depth = thetao.coord(depth_name)
    if not depth.has_bounds():
        depth.guess_bounds()
    thickness = np.diff(depth.bounds, axis=1)[:, 0]
    area = cc.cell_area(thetao)
    volume = thickness[:, None, None] * area[None, ...]
    logger.info("volcello not supplied; using layer thickness * areacello")
    return np.ma.masked_array(volume)


def _density_contours(axes, temperature_bins, salinity_bins):
    """Overlay sigma2 contours on the T-S plane."""
    if cc.gsw is None:
        return
    sal_grid, temp_grid = np.meshgrid(
        cc.bin_centres(salinity_bins), cc.bin_centres(temperature_bins))
    sigma = cc.gsw.sigma2(sal_grid, temp_grid)
    contours = axes.contour(sal_grid, temp_grid, sigma, colors="0.4",
                            linewidths=0.6, levels=12)
    axes.clabel(contours, inline=True, fontsize=7, fmt="%.1f")


def main(cfg):
    """Run the diagnostic."""
    temperature_width = float(
        cfg.get("temperature_bin_width",
                cc.DEFAULT_TEMPERATURE_BIN_WIDTH))
    salinity_width = float(
        cfg.get("salinity_bin_width",
                cc.DEFAULT_SALINITY_BIN_WIDTH))
    grouped = group_metadata(cfg["input_data"].values(), "dataset")

    datasets = sorted(grouped.items())
    if not datasets:
        logger.warning("No input data")
        return

    figure, axes_list = plt.subplots(
        1, len(datasets), figsize=(6 * len(datasets), 5.5), squeeze=False)

    for axes, (dataset, metadata) in zip(axes_list[0], datasets):
        by_variable = group_metadata(metadata, "short_name")
        if "thetao" not in by_variable or "so" not in by_variable:
            logger.warning("%s: needs both thetao and so, skipping", dataset)
            continue

        thetao = cc.load_cube(by_variable["thetao"][0]["filename"], "thetao")
        so = cc.load_cube(by_variable["so"][0]["filename"], "so")
        volcello = None
        if "volcello" in by_variable:
            volcello = cc.load_cube(
                by_variable["volcello"][0]["filename"], "volcello")

        # Construct one graph for both fields.  Iterating through time
        # slices and accessing ``cube.data`` separately would submit 240
        # small graphs one after another, leaving most workers idle.
        thetao_data, so_data = dask.compute(thetao.core_data(),
                                             so.core_data())
        thetao.data = thetao_data
        so.data = so_data
        volume = _cell_volume(thetao, volcello)

        abs_sal, cons_temp, _ = cc.teos10_fields(thetao, so)
        n_times = thetao.coord("time").shape[0] if thetao.coords("time") else 1
        logger.info("%s: histogramming %d time step(s) in one array", dataset,
                    n_times)

        # Bins from the data, as the original builds them, so the axes
        # frame the water masses instead of a fixed rectangle.
        temperature_bins = cc.ts_bin_edges(
            cons_temp.ravel(),
            temperature_width)
        salinity_bins = cc.ts_bin_edges(
            abs_sal.ravel(),
            salinity_width)

        histogram = cc.volume_weighted_ts_histogram(
            cons_temp, abs_sal, volume, temperature_bins, salinity_bins)
        histogram /= n_times
        cc.report_aliasing(histogram)
        plotted = np.ma.masked_less_equal(histogram, 0.0)

        mesh = axes.pcolormesh(
            cc.bin_centres(salinity_bins), cc.bin_centres(temperature_bins),
            plotted, cmap="viridis", shading="auto",
            norm=mcolors.LogNorm(
                vmin=max(np.ma.min(plotted), 1.0) if plotted.count() else 1.0,
                vmax=max(np.ma.max(plotted), 10.0) if plotted.count() else 10.0
            ))
        _density_contours(axes, temperature_bins, salinity_bins)
        axes.set_xlabel("Absolute salinity (g kg$^{-1}$)")
        axes.set_ylabel("Conservative temperature ($^\\circ$C)")
        axes.set_title(f"{dataset}\nvolume-weighted T-S", fontsize=11)
        figure.colorbar(mesh, ax=axes, shrink=0.8, pad=0.03,
                        label="Volume (m$^3$)")

        out = cc.make_cube(
            histogram,
            [(cc.bin_centres(temperature_bins),
              "sea_water_conservative_temperature", "degC"),
             (cc.bin_centres(salinity_bins),
              "sea_water_absolute_salinity", "g kg-1")],
            "ts_volume", "m3", "volume-weighted temperature-salinity "
                               "distribution")
        path = get_diagnostic_filename(
            f"cosima_ts_diagram_{dataset}", cfg)
        iris.save(out, path)
        record = {
            "caption": f"Volume-weighted temperature-salinity distribution "
                       f"of {dataset} over the region selected by the "
                       f"recipe.",
            "statistics": ["other"],
            "domains": ["global"],
            "plot_types": ["scatter"],
            "authors": ["beucher_romain"],
            "references": ["cosima_recipes"],
            "ancestors": [m["filename"] for m in metadata],
        }
        with ProvenanceLogger(cfg) as logger_:
            logger_.log(path, record)

    plot_path = get_plot_filename("cosima_temperature_salinity_diagram", cfg)
    figure.savefig(plot_path, bbox_inches="tight", dpi=150)
    plt.close(figure)
    with ProvenanceLogger(cfg) as logger_:
        logger_.log(plot_path, {
            **record,
            "caption": "Volume-weighted temperature-salinity distributions by model.",
            "ancestors": [m["filename"] for m in cfg["input_data"].values()],
        })
    logger.info("Wrote %s", plot_path)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
