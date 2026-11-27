"""Meridional overturning circulation (COSIMA Overturning_Circulation).

Translation of ``03-Advanced-Recipes/Overturning_Circulation.ipynb``,
merged with ``Transformation_from_Depth_to_Potential_Density.ipynb``
as the inventory recommends, since the two notebooks are the depth- and
density-space halves of one diagnostic.

Depth space
    CMIP6 publishes the streamfunction directly as ``msftyz`` (or
    ``msftmz``), so the diagnostic only selects the basin, converts to
    Sv and time-averages.

Density space
    MOM bins the meridional transport into density layers online, into
    ``ty_trans_rho``.  CMIP6 has no such variable, so this rebins the
    z-coordinate ``vmo`` using sigma2 computed from monthly ``thetao``
    and ``so``.  Rebinning monthly means cannot recover the correlation
    between the transport and the density within a month, so the eddy
    contribution is underestimated relative to an online binning.  This
    is a known methodological difference, not an equivalent
    calculation, and the plot is labelled accordingly.

Both are integrated from the sea floor upwards, so the streamfunction
vanishes at the bottom.
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

#: Default sigma2 bins.  Deliberately wider than the range worth
#: plotting: a cell whose density falls outside the bins contributes
#: nothing, which offsets that latitude's streamfunction by a constant
#: at *every* density and shows up as a vertical stripe across the
#: section.  Brackish water (river plumes, the Baltic, sea-ice melt)
#: reaches sigma2 below 20, so the bins start well below anything the
#: open ocean occupies and :data:`DEFAULT_PLOT_RANGE` narrows the view
#: afterwards.
DEFAULT_SIGMA2_BINS = np.arange(15.0, 38.51, 0.05)

#: Density range actually shown, as (densest, lightest).  The COSIMA
#: original does the same thing, plotting only part of its rho2 axis.
DEFAULT_PLOT_RANGE = (37.8, 29.5)

#: Latitudes north of this are excluded: on a tripolar grid the local
#: y-direction stops being northward, so a zonal sum of vmo is wrong
#: there.  See COSIMA/cosima-recipes#510.
MAX_LATITUDE = 65.0


def streamfunction_depth_space(cube):
    """Depth-space streamfunction in Sv from ``msftyz``/``msftmz``."""
    cube = cc.time_mean(cube)
    data = cc.masked_data(cube)
    while data.ndim > 2:
        # Collapse any leading basin dimension by taking the first
        # (global) basin, which is what the notebook plots.
        data = data[0]
    if str(cube.units) in ("kg s-1", "kg/s"):
        data = (data / cc.RHO0) / 1.0e6
    return data


def streamfunction_density_space(vmo_cube, sigma2, bins):
    """Density-space streamfunction in Sv from ``vmo`` and sigma2.

    ``vmo`` is northward mass transport in kg s-1 on ``(depth, lat,
    lon)``.  The transport is binned by the density of the same cell,
    summed zonally, converted to Sv and integrated from the densest
    class upwards.
    """
    cc.check_transport_grid(vmo_cube, "vmo")

    transport = cc.masked_data(vmo_cube)
    if transport.ndim == 4:  # a time axis survived the preprocessor
        transport = np.ma.mean(transport, axis=0)

    report_bin_coverage(sigma2, bins, transport)

    # Bin over depth, keeping (lat, lon); then sum zonally.
    binned = cc.rebin_to_density(transport, sigma2, bins, depth_axis=0)
    zonal = np.ma.sum(np.ma.masked_invalid(binned), axis=-1)
    zonal = (zonal / cc.RHO0) / 1.0e6
    psi = cc.cumulative_from_bottom(zonal, axis=0)
    report_closure(psi)
    return psi


def report_bin_coverage(sigma2, bins, transport):
    """Warn when transport falls outside the density bins.

    Such a cell contributes nothing to any bin, so the latitude it sits
    at ends up offset by a constant at every density.  On the section
    that is a vertical stripe, which is easy to mistake for a feature.
    """
    inside = (sigma2 >= bins[0]) & (sigma2 <= bins[-1])
    weight = np.ma.abs(transport)
    total = float(np.ma.sum(weight))
    if total <= 0:
        return 1.0
    kept = float(np.ma.sum(np.ma.where(inside, weight, 0.0))) / total
    if kept < 0.999:
        logger.warning(
            "%.2f%% of the transport falls outside the sigma2 bins "
            "%.2f to %.2f (density spans %.2f to %.2f); those cells are "
            "dropped and will appear as vertical stripes. Widen "
            "sigma2_bins.", 100 * (1 - kept), bins[0], bins[-1],
            float(np.ma.min(sigma2)), float(np.ma.max(sigma2)))
    else:
        logger.info("%.3f%% of the transport is inside the sigma2 bins",
                    100 * kept)
    return kept


def report_closure(psi):
    """Report how close the streamfunction comes to closing.

    Integrating from the densest class up, psi is zero at the sea floor
    by construction; it should also return to ~0 in the lightest
    classes, because the net transport across a latitude is close to
    zero.  A large residual means transport is missing.
    """
    residual = float(np.ma.max(np.ma.abs(psi[0])))
    message = "streamfunction closure residual at the lightest class: %.2f Sv"
    if residual > 5.0:
        logger.warning(
            message + " -- the net transport does not balance, so some "
            "transport is missing (check the bin coverage above)", residual)
    else:
        logger.info(message, residual)
    return residual


def _plot_section(axes, x, y, field, title, ylabel, invert=True,
                  max_value=25.0, ylim=None):
    """Filled contour plot of a streamfunction section."""
    levels = np.hstack((np.arange(-max_value, 0, 2),
                        np.flip(-np.arange(-max_value, 0, 2))))
    mesh = axes.contourf(x, y, field, levels=levels, cmap="RdBu_r",
                         extend="both")
    axes.contour(x, y, field, levels=levels, colors="k", linewidths=0.25)
    if invert:
        axes.invert_yaxis()
    if ylim is not None:
        # The bins run wider than this on purpose; see
        # DEFAULT_SIGMA2_BINS. Narrow the view, not the binning.
        axes.set_ylim(*ylim)
    axes.set_xlabel("Latitude")
    axes.set_ylabel(ylabel)
    axes.set_title(title, fontsize=11)
    axes.set_xlim(-75, MAX_LATITUDE)
    return mesh


def _provenance(cfg, path, caption, ancestors):
    record = {
        "caption": caption,
        "statistics": ["mean"],
        "domains": ["global"],
        "plot_types": ["zonal"],
        "authors": ["beucher_romain"],
        "references": ["cosima_recipes"],
        "ancestors": ancestors,
    }
    with ProvenanceLogger(cfg) as logger_:
        logger_.log(path, record)


def main(cfg):
    """Run the diagnostic."""
    bins = np.asarray(cfg.get("sigma2_bins", DEFAULT_SIGMA2_BINS), dtype=float)
    max_psi = float(cfg.get("max_streamfunction", 25.0))
    plot_range = tuple(cfg.get("density_plot_range", DEFAULT_PLOT_RANGE))
    grouped = group_metadata(cfg["input_data"].values(), "dataset")

    for index, (dataset, metadata) in enumerate(sorted(grouped.items())):
        by_variable = group_metadata(metadata, "short_name")
        ancestors = [m["filename"] for m in metadata]
        panels = []

        # ---- depth space -------------------------------------------------
        for name in ("msftyz", "msftmz"):
            if name not in by_variable:
                continue
            cube = cc.load_cube(by_variable[name][0]["filename"], name)
            psi = streamfunction_depth_space(cube)
            depth_name = cc.depth_coord_name(cube)
            depth = cube.coord(depth_name).points if depth_name else None
            lat = cc.row_latitude(cube)
            panels.append(("depth", lat, depth, psi,
                           f"{dataset}: depth space ({name})",
                           "Depth (m)"))

            out = cc.make_cube(
                np.ma.filled(psi, np.nan),
                [(depth, "depth", "m"), (lat, "latitude", "degrees_north")],
                "msftyz", "Sv", "meridional overturning streamfunction")
            path = get_diagnostic_filename(
                f"cosima_moc_depth_{dataset}", cfg)
            iris.save(out, path)
            _provenance(cfg, path,
                        f"Depth-space meridional overturning streamfunction "
                        f"of {dataset}.", ancestors)
            break

        # ---- density space -----------------------------------------------
        if all(name in by_variable for name in ("vmo", "thetao", "so")):
            vmo = cc.load_cube(by_variable["vmo"][0]["filename"], "vmo")
            thetao = cc.load_cube(by_variable["thetao"][0]["filename"],
                                  "thetao")
            so = cc.load_cube(by_variable["so"][0]["filename"], "so")
            sigma2 = cc.potential_density_anomaly(thetao, so, p_ref=2000)
            psi_rho = streamfunction_density_space(vmo, sigma2, bins)
            lat = cc.row_latitude(vmo)
            centres = cc.bin_centres(bins)
            panels.append(("density", lat, centres, psi_rho,
                           f"{dataset}: density space (offline rebin)",
                           r"$\sigma_2$ (kg m$^{-3}$)"))

            out = cc.make_cube(
                np.ma.filled(psi_rho, np.nan),
                [(centres, "sea_water_sigma_theta", "kg m-3"),
                 (lat, "latitude", "degrees_north")],
                "msftyzrho", "Sv",
                "meridional overturning streamfunction in sigma2 space")
            path = get_diagnostic_filename(
                f"cosima_moc_density_{dataset}", cfg)
            iris.save(out, path)
            _provenance(
                cfg, path,
                f"Density-space meridional overturning streamfunction of "
                f"{dataset}, obtained by rebinning period-mean vmo with "
                f"period-mean sigma2. Seasonal transport-density covariance "
                f"and the online eddy contribution are omitted.", ancestors)
        elif "vmo" in by_variable:
            logger.warning("%s: vmo present but thetao/so missing, "
                           "skipping the density-space streamfunction",
                           dataset)

        if not panels:
            logger.warning("%s: nothing to plot", dataset)
            continue

        figure, axes_list = plt.subplots(
            len(panels), 1, figsize=(10, 5 * len(panels)), squeeze=False)
        for axes, (kind, lat, yaxis, field, title, ylabel) in zip(
                axes_list[:, 0], panels):
            mesh = _plot_section(
                axes, lat, yaxis, field, title, ylabel,
                max_value=max_psi,
                ylim=plot_range if kind == "density" else None)
            figure.colorbar(mesh, ax=axes, shrink=0.85, pad=0.03,
                            label="Overturning (Sv)")
            if kind == "density":
                axes.text(0.01, -0.16,
                          "Offline rebin of period means: seasonal "
                          "transport-density covariance is omitted.",
                          transform=axes.transAxes, fontsize=8, color="0.3")
        figure.suptitle(
            f"Meridional overturning circulation: {dataset}\n"
            f"(latitudes north of {MAX_LATITUDE:.0f}N excluded: "
            f"tripolar grid)", fontsize=12)
        plot_path = get_plot_filename(f"cosima_moc_{dataset}", cfg)
        figure.savefig(plot_path, bbox_inches="tight", dpi=150)
        plt.close(figure)
        _provenance(
            cfg, plot_path,
            f"Depth- and density-space meridional overturning sections for "
            f"{dataset}; density space rebins period means.", ancestors)
        logger.info("Wrote %s", plot_path)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
