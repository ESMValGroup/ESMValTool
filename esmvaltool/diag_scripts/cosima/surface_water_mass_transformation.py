"""Surface water-mass transformation (COSIMA
Surface_Water_Mass_Transformation).

Translation of ``03-Advanced-Recipes/Surface_Water_Mass_Transformation.ipynb``
onto CMORised CMIP6 data, following Newsom et al. (2016) and Abernathey
et al. (2016).  Surface water-mass transformation is the volume flux
into a density class driven by surface buoyancy forcing:

    Omega(sigma) = d/dsigma int_{sigma' < sigma}
                       (dsigma/dtheta * theta_dot
                        + dsigma/dS * S_dot) dx dy

which separates into a heat and a salt contribution:

    heat transformation  =  alpha * Q / c_p
    salt transformation  =  beta * (SA * wfo - 1000 * sfdsi)

Both are accumulated into sigma1 bins using the surface density of the
same cell and the same month, then normalised by the total number of
days and by the bin width, and finally sign-flipped so that positive
means conversion towards denser water.

Variable mapping from the notebook's MOM6 diagnostics:

===================== ==========================================
notebook              CMIP6
===================== ==========================================
``hfds``              ``hfds`` (Omon)
``wfo``               ``wfo`` (Omon)
``salt_flux``         ``sfdsi`` (Omon), downward salt flux
``thetao``, ``so``    ``tos``, ``sos`` (Omon) surface fields
``areacello``         ``areacello`` (Ofx)
===================== ==========================================

``sfdsi`` is the sea-ice salt flux only; models that also restore
salinity fold that term elsewhere, so the salt contribution is a lower
bound where restoring matters.
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

#: sigma1 bins; the notebook's np.arange(31, 33.5, 0.02).
DEFAULT_SIGMA1_BINS = np.arange(31.0, 33.5, 0.02)

SECONDS_PER_DAY = 86400.0


def _days_per_month(cube):
    """Length of each month in days, in the model's own calendar.

    See :func:`cosima_common.time_interval_days`: this must not go
    through a hard-coded reference unit, both because ``cf_units``
    refuses to convert between calendars and because the month lengths
    themselves are calendar-dependent.
    """
    return cc.time_interval_days(cube)


def surface_transformation(tos, sos, hfds, wfo, sfdsi=None,
                           temperature_is_conservative=False):
    """Heat and salt transformation fields, in kg m-2 per bin-normalising.

    Returns ``(heat, salt, sigma1, days)``, each of the first three
    shaped ``(time, lat, lon)``, with ``days`` the month lengths.
    """
    salinity = cc.masked_data(sos)
    temperature = cc.masked_data(tos)
    lat = np.broadcast_to(cc.lat_2d(sos)[None, ...], salinity.shape)
    lon = np.broadcast_to(cc.lon_2d(sos)[None, ...], salinity.shape)
    pressure = np.zeros_like(salinity)  # surface fields

    abs_sal = cc.absolute_salinity(salinity, pressure, lon, lat)
    if temperature_is_conservative:
        cons_temp = temperature
    else:
        cons_temp = cc.conservative_temperature(abs_sal, temperature)

    alpha = cc.thermal_expansion(abs_sal, cons_temp, pressure)
    beta = cc.haline_contraction(abs_sal, cons_temp, pressure)
    with cc.quiet_gsw():
        sigma1 = np.ma.masked_invalid(
            cc.gsw.sigma1(abs_sal, cons_temp))

    days = _days_per_month(tos)
    days_field = days[:, None, None]

    heat_flux = cc.masked_data(hfds)
    water_flux = cc.masked_data(wfo)
    # wfo is kg m-2 s-1 of water; multiplying by SA turns it into the
    # equivalent salt flux in g m-2 s-1.
    salt_from_water = abs_sal * water_flux
    if sfdsi is not None:
        # sfdsi is kg m-2 s-1 of salt; beta is in kg/g, so convert.
        salt_flux = cc.masked_data(sfdsi) * 1000.0
    else:
        salt_flux = np.ma.zeros(salt_from_water.shape)

    heat = alpha * heat_flux * days_field
    salt = beta * (salt_from_water - salt_flux) * days_field
    return heat, salt, sigma1, days


def report_bin_coverage(sigma1, bins):
    """Log how much of the domain falls inside the density bins.

    Cells outside the bin range contribute nothing, silently.  The
    default bins are tuned for the Southern Ocean south of 60S; point
    the recipe at a wider region and most of the water leaves the range
    without any sign that the answer is now a partial one.
    """
    inside = (sigma1 >= bins[0]) & (sigma1 <= bins[-1])
    fraction = float(np.ma.mean(np.ma.filled(inside, False).astype(float)))
    message = ("%.0f%% of cells fall inside the sigma1 bin range "
               "%.2f to %.2f (data spans %.2f to %.2f)")
    arguments = (100 * fraction, bins[0], bins[-1],
                 float(np.ma.min(sigma1)), float(np.ma.max(sigma1)))
    if fraction < 0.9:
        logger.warning(
            message + "; the transformation only accounts for that "
            "fraction of the region. Widen sigma1_bins or narrow the "
            "region.", *arguments)
    else:
        logger.info(message, *arguments)
    return fraction


def bin_transformation(heat, salt, sigma1, days, bins):
    """Bin the transformation fields by surface density.

    Returns ``(heat_binned, salt_binned)`` of shape ``(bin, lat, lon)``
    in m s-1, normalised by the total number of days and the bin width,
    with the sign flipped so that positive is densification.
    """
    heat_binned = cc.bin_sum(heat, sigma1, bins, axis=0)
    salt_binned = cc.bin_sum(salt, sigma1, bins, axis=0)

    total_days = float(np.sum(days))
    widths = np.diff(np.asarray(bins, dtype=float))[:, None, None]

    heat_binned = heat_binned / (cc.CP_OCEAN * total_days) / widths
    salt_binned = salt_binned / total_days / widths

    heat_binned = np.ma.masked_equal(heat_binned, 0.0)
    salt_binned = np.ma.masked_equal(salt_binned, 0.0)
    return -heat_binned, -salt_binned


def integrate_to_sverdrups(binned, area):
    """Horizontally integrate a binned transformation rate into Sv."""
    weighted = binned * np.broadcast_to(area, binned.shape)
    return np.ma.sum(np.ma.sum(weighted, axis=-1), axis=-1) / 1.0e6


def main(cfg):
    """Run the diagnostic."""
    bins = np.asarray(cfg.get("sigma1_bins", DEFAULT_SIGMA1_BINS),
                      dtype=float)
    centres = cc.bin_centres(bins)
    grouped = group_metadata(cfg["input_data"].values(), "dataset")

    figure, axes = plt.subplots(figsize=(7, 7))
    plotted = 0

    for index, (dataset, metadata) in enumerate(sorted(grouped.items())):
        by_variable = group_metadata(metadata, "short_name")
        required = ("tos", "sos", "hfds", "wfo")
        missing = [name for name in required if name not in by_variable]
        if missing:
            logger.warning("%s: missing %s, skipping",
                           dataset, ", ".join(missing))
            continue

        cubes = {name: cc.load_cube(by_variable[name][0]["filename"], name)
                 for name in required}
        sfdsi = None
        if "sfdsi" in by_variable:
            sfdsi = cc.load_cube(by_variable["sfdsi"][0]["filename"], "sfdsi")
        else:
            logger.warning("%s: no sfdsi, the salt contribution omits the "
                           "sea-ice salt flux", dataset)

        heat, salt, sigma1, days = surface_transformation(
            cubes["tos"], cubes["sos"], cubes["hfds"], cubes["wfo"], sfdsi)
        report_bin_coverage(sigma1, bins)
        heat_binned, salt_binned = bin_transformation(
            heat, salt, sigma1, days, bins)

        area = cc.cell_area(cubes["tos"])
        swmt_heat = integrate_to_sverdrups(heat_binned, area)
        swmt_salt = integrate_to_sverdrups(salt_binned, area)
        swmt_net = swmt_heat + swmt_salt

        style = cc.style_for(dataset, index)
        axes.plot(swmt_net, centres, color=style["color"], linestyle="-",
                  linewidth=2.0, label=f"{dataset}: net")
        axes.plot(swmt_heat, centres, color=style["color"], linestyle="--",
                  linewidth=1.3, label=f"{dataset}: heat")
        axes.plot(swmt_salt, centres, color=style["color"], linestyle=":",
                  linewidth=1.3, label=f"{dataset}: salt")
        plotted += 1

        out = iris.cube.CubeList([
            cc.make_cube(np.ma.filled(values, np.nan),
                         [(centres, "sea_water_sigma_theta", "kg m-3")],
                         var_name, "Sv", long_name)
            for values, var_name, long_name in (
                (swmt_net, "swmt", "surface water-mass transformation"),
                (swmt_heat, "swmt_heat",
                 "surface water-mass transformation from heat fluxes"),
                (swmt_salt, "swmt_salt",
                 "surface water-mass transformation from salt fluxes"),
            )])
        path = get_diagnostic_filename(f"cosima_swmt_{dataset}", cfg)
        iris.save(out, path)
        record = {
            "caption": f"Surface water-mass transformation of {dataset} "
                       f"south of the region selected by the recipe, "
                       f"partitioned into heat and salt contributions.",
            "statistics": ["mean"],
            "domains": ["shpolar"],
            "plot_types": ["line"],
            "authors": ["beucher_romain"],
            "references": ["cosima_recipes", "newsom16jclim"],
            "ancestors": [m["filename"] for m in metadata],
        }
        with ProvenanceLogger(cfg) as logger_:
            logger_.log(path, record)

    if not plotted:
        logger.warning("No dataset had the full set of surface fluxes")
        plt.close(figure)
        return

    axes.axvline(0.0, color="black", linewidth=0.6)
    axes.set_ylim(centres.max(), centres.min())
    axes.set_ylabel(r"$\sigma_1$ (kg m$^{-3}$)", fontsize=12)
    axes.set_xlabel("Surface water-mass transformation (Sv)", fontsize=12)
    axes.set_title("Surface water-mass transformation\n"
                   "positive = densification", fontsize=12)
    axes.legend(loc="lower left", fontsize=8, frameon=False)
    axes.grid(alpha=0.3)

    plot_path = get_plot_filename("cosima_surface_water_mass_transformation",
                                  cfg)
    figure.savefig(plot_path, bbox_inches="tight", dpi=150)
    plt.close(figure)
    with ProvenanceLogger(cfg) as logger_:
        logger_.log(plot_path, {
            **record,
            "caption": "Surface water-mass transformation by model and flux contribution.",
            "ancestors": [m["filename"] for m in cfg["input_data"].values()],
        })
    logger.info("Wrote %s", plot_path)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
