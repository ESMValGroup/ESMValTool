"""Shared machinery for russell18jgr figures 6a and 6b.

Density-layer based volume (6a) and heat (6b) transport across 30S,
with layers defined as in Talley (2003, J. Phys. Oceanogr. 33,
530-560).  Port of the common parts of russell18jgr-fig6a.ncl and
russell18jgr-fig6b.ncl.
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

# Layer bounds: each condition is (sigma_variable, lower, upper) and all
# conditions of a layer are combined with AND.  sigma variables are
# 's0', 's2', 's4' (potential density anomaly referenced to 0, 2000 and
# 4000 dbar).  Bounds of None are unbounded.
MAIN_LAYERS = [
    [("s0", None, 26.10)],
    [("s0", 26.10, 26.40)],
    [("s0", 26.40, 26.90)],
    [("s0", 26.90, 27.10)],
    [("s0", 27.10, 27.40)],
    [("s0", 27.40, None), ("s2", None, 36.80)],
    [("s2", 36.80, None), ("s4", None, 45.80)],
    [("s4", 45.80, 45.86)],
    [("s4", 45.86, 45.92)],
    [("s4", 45.92, 46.00)],
    [("s4", 46.00, None)],
]

SUB_LAYERS = [
    # layer 1
    [("s0", None, 24.900)],
    [("s0", 24.900, 25.300)],
    [("s0", 25.300, 25.700)],
    [("s0", 25.700, 26.100)],
    # layer 2
    [("s0", 26.100, 26.175)],
    [("s0", 26.175, 26.250)],
    [("s0", 26.250, 26.325)],
    [("s0", 26.325, 26.400)],
    # layer 3
    [("s0", 26.400, 26.525)],
    [("s0", 26.525, 26.650)],
    [("s0", 26.650, 26.775)],
    [("s0", 26.775, 26.900)],
    # layer 4
    [("s0", 26.900, 26.950)],
    [("s0", 26.950, 27.000)],
    [("s0", 27.000, 27.050)],
    [("s0", 27.050, 27.100)],
    # layer 5
    [("s0", 27.100, 27.175)],
    [("s0", 27.175, 27.250)],
    [("s0", 27.250, 27.325)],
    [("s0", 27.325, 27.400)],
    # layer 6
    [("s0", 27.400, 27.500)],
    [("s0", 27.500, None), ("s2", None, 36.700)],
    [("s2", 36.700, 36.750)],
    [("s2", 36.750, 36.800)],
    # layer 7
    [("s2", 36.800, 36.850)],
    [("s2", 36.850, 36.900)],
    [("s2", 36.900, 36.950)],
    [("s2", 36.950, None), ("s4", None, 45.800)],
    # layer 8
    [("s4", 45.800, 45.815)],
    [("s4", 45.815, 45.830)],
    [("s4", 45.830, 45.845)],
    [("s4", 45.845, 45.860)],
    # layer 9
    [("s4", 45.860, 45.875)],
    [("s4", 45.875, 45.890)],
    [("s4", 45.890, 45.905)],
    [("s4", 45.905, 45.920)],
    # layer 10
    [("s4", 45.920, 45.940)],
    [("s4", 45.940, 45.960)],
    [("s4", 45.960, 45.980)],
    [("s4", 45.980, 46.000)],
    # layer 11
    [("s4", 46.000, 46.050)],
    [("s4", 46.050, 46.100)],
    [("s4", 46.100, 46.150)],
    [("s4", 46.150, None)],
]

YAXIS_LABELS = [
    "Net",
    "Surface",
    r"26.10$\sigma_0$",
    r"26.40$\sigma_0$",
    r"26.90$\sigma_0$",
    r"27.10$\sigma_0$",
    r"27.40$\sigma_0$",
    r"36.8$\sigma_2$",
    r"45.8$\sigma_4$",
    r"45.86$\sigma_4$",
    r"45.92$\sigma_4$",
    r"46.0$\sigma_4$",
    "Bottom",
]


def _interp_to_lat(cube, exact_lat):
    """Linear interpolation of a (lev, lat, lon) field to one latitude."""
    data = np.ma.masked_invalid(cube.data)
    lat = rc.lat_1d(cube)
    if lat[0] > lat[-1]:  # ensure ascending
        lat = lat[::-1]
        data = data[:, ::-1, :]
    k = int(np.searchsorted(lat, exact_lat)) - 1
    k = int(np.clip(k, 0, len(lat) - 2))
    weight = (exact_lat - lat[k]) / (lat[k + 1] - lat[k])
    weight = float(np.clip(weight, 0.0, 1.0))
    return (1.0 - weight) * data[:, k, :] + weight * data[:, k + 1, :]


def _interp_to_lon(field, src_lon, dst_lon):
    """Cyclic linear interpolation of a (lev, lon) field in longitude."""
    src_lon = np.mod(np.asarray(src_lon, dtype=float), 360.0)
    dst_lon = np.mod(np.asarray(dst_lon, dtype=float), 360.0)
    if not (np.all(np.diff(src_lon) > 0) and len(src_lon) > 1):
        return field  # not monotonic: use as-is (as the NCL script did)
    out = np.ma.masked_all((field.shape[0], len(dst_lon)))
    src_ext = np.concatenate(
        [src_lon[-1:] - 360.0, src_lon, src_lon[:1] + 360.0]
    )
    for k in range(field.shape[0]):
        row = np.ma.filled(field[k], np.nan)
        row_ext = np.concatenate([row[-1:], row, row[:1]])
        out[k] = np.ma.masked_invalid(np.interp(dst_lon, src_ext, row_ext))
    return out


def prepare_section(thetao_cube, so_cube, vo_cube, volcello_cube):
    """Compute everything needed at the 30S section of the vo grid.

    Returns a dict with theta (degC), sigma0/2/4, vo, cross-section
    area (m^2), all with shape (lev, lon), plus exact_lat.
    """
    theta = rc.time_mean(thetao_cube)
    salt = rc.time_mean(so_cube)
    vo = rc.time_mean(vo_cube)

    vo_lat = rc.lat_1d(vo)
    vo_lon = rc.lon_1d(vo)
    a = rc.closest_index(-30.0, vo_lat)
    exact_lat = float(vo_lat[a])

    theta_sec = _interp_to_lat(theta, exact_lat)
    so_sec = _interp_to_lat(salt, exact_lat)
    theta_sec = _interp_to_lon(theta_sec, rc.lon_1d(theta), vo_lon)
    so_sec = _interp_to_lon(so_sec, rc.lon_1d(salt), vo_lon)

    if np.ma.max(theta_sec) > 250.0:  # make sure temperature is in C
        theta_sec = theta_sec - 273.15

    sigma = {}
    for name, depth in (("s0", 0.0), ("s2", 1977.0), ("s4", 3948.0)):
        rho = rc.rho_mwjf(theta_sec, so_sec, depth)
        sigma[name] = 1000.0 * (rho - 1.0)

    # north-south cross-section area of each cell at the section
    volcello = np.ma.masked_invalid(volcello_cube.data)
    vol_lev = volcello_cube.coord(axis="Z").points
    dlat = abs(vo_lat[a] - vo_lat[a + 1]) * rc.EARTH_RADIUS * rc.DEG2RAD
    area = volcello / dlat
    if abs(vol_lev[0]) > abs(vol_lev[1]):
        area_sec = area[::-1, a, :]
    else:
        area_sec = area[:, a, :]

    vo_sec = np.ma.masked_invalid(vo.data)[:, a, :]

    return {
        "theta": theta_sec,
        "sigma": sigma,
        "vo": vo_sec,
        "area": area_sec,
        "exact_lat": exact_lat,
    }


def layer_mask(sigma, conditions):
    """Boolean mask for one (sub)layer definition."""
    mask = np.ones(sigma["s0"].shape, dtype=bool)
    for var, lower, upper in conditions:
        values = sigma[var]
        if lower is not None:
            mask &= np.ma.filled(values >= lower, False)
        if upper is not None:
            mask &= np.ma.filled(values < upper, False)
    return mask


def layer_sums(sigma, transport):
    """Transport summed per main layer (12 incl. net) and sublayer (44)."""
    main = np.zeros(12)
    for i, conditions in enumerate(MAIN_LAYERS):
        masked = np.ma.masked_where(~layer_mask(sigma, conditions), transport)
        value = masked.sum()
        main[i + 1] = 0.0 if np.ma.is_masked(value) else float(value)
    main[0] = main[1:].sum()

    sub = np.zeros(44)
    for i, conditions in enumerate(SUB_LAYERS):
        masked = np.ma.masked_where(~layer_mask(sigma, conditions), transport)
        value = masked.sum()
        sub[i] = 0.0 if np.ma.is_masked(value) else float(value)
    return main, sub


def plot_layer_bars(
    main, sub, talley, dataset, meta, exact_lat, xlimit, xstep, unit, net_label
):
    """One bar-chart figure (blue main bars, red sublayers, magenta
    Talley reference).
    """
    fig, axes = plt.subplots(figsize=(7, 9))
    y_edges = np.arange(-1.0, 12.0)  # bottom edges of the 12 main bars

    for i in range(12):
        axes.barh(
            y_edges[i] + 0.5,
            main[i],
            height=1.0,
            color="#00008b",
            edgecolor="#00008b",
            zorder=2,
        )
    for i in range(44):
        layer = i // 4
        offset = (i % 4) * 0.25
        ypos = layer + offset + 0.125
        axes.barh(
            ypos,
            sub[i],
            height=0.25,
            color="red",
            edgecolor="black",
            linewidth=0.3,
            zorder=3,
        )
    for i in range(12):
        axes.plot(
            [talley[i], talley[i]],
            [y_edges[i], y_edges[i] + 1.0],
            color="magenta",
            linewidth=2.0,
            zorder=4,
        )
    for i in range(12):
        axes.text(
            xlimit * 0.98,
            y_edges[i] + 0.5,
            f"{main[i]:4.3f}",
            fontsize=7,
            va="center",
            ha="right",
            zorder=5,
        )

    axes.axvline(0, color="black", linewidth=0.8)
    axes.set_xlim(-xlimit, xlimit)
    axes.set_xticks(np.arange(-xlimit, xlimit + xstep / 2, xstep))
    axes.set_ylim(11, -1)
    axes.set_yticks(np.arange(-1, 12))
    axes.set_yticklabels(YAXIS_LABELS, fontsize=8)
    axes.tick_params(axis="x", labelsize=7)
    axes.set_title(dataset, fontsize=12)
    axes.text(
        0.0,
        1.02,
        f"({meta['start_year']} - {meta['end_year']}) at "
        f"({abs(exact_lat):4.2f}S)",
        transform=axes.transAxes,
        fontsize=8,
    )
    axes.text(
        1.0,
        1.02,
        f"{net_label} = {main[0]:4.2f}{unit}",
        ha="right",
        transform=axes.transAxes,
        fontsize=8,
    )
    fig.tight_layout()
    return fig


def run_fig6(cfg, mode):
    """Run figure 6a (mode='volume') or 6b (mode='heat').

    The ESMValTool imports live here rather than at module level so
    that the notebooks can reuse the functions above with only
    ESMValCore installed.
    """
    from esmvaltool.diag_scripts.shared import (
        ProvenanceLogger,
        get_diagnostic_filename,
        get_plot_filename,
    )

    metadata = list(cfg["input_data"].values())
    vo_data = rc.select_by_short_name(metadata, "vo")
    so_data = rc.select_by_short_name(metadata, "so")
    thetao_data = rc.select_by_short_name(metadata, "thetao")
    vol_data = rc.select_by_short_name(metadata, "volcello")

    if mode == "volume":
        talley = [
            0.0,
            -7.34,
            -1.9,
            9.71,
            2.37,
            -5.86,
            -10.02,
            -11.11,
            -2.84,
            9.96,
            16.49,
            0.51,
        ]
        xlimit, xstep, unit = 20.0, 4.0, "Sv"
        net_label = "Net Transport out of Southern ocean"
        var_name = "transport_per_layer"
        tag, caption = "6a", "Russell et al 2018 figure 6 part a"
    else:
        talley = [
            -0.91,
            -0.89,
            -0.13,
            0.43,
            0.04,
            -0.16,
            -0.14,
            -0.12,
            -0.03,
            0.05,
            0.05,
            0.00,
        ]
        xlimit, xstep, unit = 1.4, 0.2, "PW"
        net_label = "Net energy out of Southern ocean"
        var_name = "energy_transport_per_layer"
        tag, caption = "6b", "Russell et al 2018 figure 6b"

    with ProvenanceLogger(cfg) as prov:
        for vo_meta in vo_data:
            dataset = vo_meta["dataset"]
            thetao_meta = rc.match_by_dataset(thetao_data, dataset)
            so_meta = rc.match_by_dataset(so_data, dataset)
            vol_meta = rc.match_by_dataset(vol_data, dataset)
            if None in (thetao_meta, so_meta, vol_meta):
                raise ValueError(
                    f"thetao/so/volcello for {dataset} not all found; "
                    "please keep the four variable groups consistent."
                )

            section = prepare_section(
                rc.load_cube(thetao_meta["filename"], "thetao"),
                rc.load_cube(so_meta["filename"], "so"),
                rc.load_cube(vo_meta["filename"], "vo"),
                rc.load_cube(vol_meta["filename"], "volcello"),
            )

            if mode == "volume":
                transport = section["vo"] * section["area"] / 1.0e6
            else:
                # 4.2 kJ/(kg K) specific heat, 1035 kg/m^3 density,
                # 1e12 converts to PW
                transport = (
                    section["vo"]
                    * section["area"]
                    * section["theta"]
                    * (1035.0 * 4.2)
                    / 1.0e12
                )

            main, sub = layer_sums(section["sigma"], transport)
            fig = plot_layer_bars(
                main,
                sub,
                talley,
                dataset,
                vo_meta,
                section["exact_lat"],
                xlimit,
                xstep,
                unit,
                net_label,
            )
            plot_file = get_plot_filename(
                f"russell18jgr-fig{tag}_{dataset}_"
                f"{vo_meta['start_year']}-{vo_meta['end_year']}",
                cfg,
            )
            fig.savefig(plot_file, bbox_inches="tight", dpi=200)
            plt.close(fig)
            logger.info("Wrote %s", plot_file)

            out = iris.cube.Cube(
                np.concatenate([main, sub]),
                var_name=var_name,
                long_name=(
                    "transport in main layers (blue bars, i=0-11) and "
                    "sub layers (red bars, i=12-55) at "
                    f"{abs(section['exact_lat']):.2f}S"
                ),
            )
            nc_name = get_diagnostic_filename(
                f"russell18jgr-figure{tag}_{dataset}_"
                f"{vo_meta['start_year']}-{vo_meta['end_year']}",
                cfg,
            )
            iris.save(out, nc_name)

            record = rc.provenance_record(
                caption,
                [
                    vo_meta["filename"],
                    thetao_meta["filename"],
                    so_meta["filename"],
                ],
                plot_types=("bar", "vert"),
            )
            prov.log(plot_file, record)
            prov.log(nc_name, record)
