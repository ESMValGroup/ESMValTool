"""Zonal velocity through Drake Passage (russell18jgr figure 4).

Python port of russell18jgr-fig4.ncl.

- Takes the time average of uo and extracts the longitude column
  closest to 69W (291E).
- If a volcello file is available, computes the total transport through
  the passage: volcello divided by the east-west grid distance gives
  the cross-section area of each cell, transport per cell is
  uo * area / 1e6 (Sv), summed between 76S and 48S with the NCL
  1/cos(lat) compensation.
- Plots the velocity section (depth vs latitude) per dataset, paged
  max_vert x max_hori, with the net transport in the panel title.
"""

import logging
import os
import sys

import iris
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    run_diagnostic,
)

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import russell_common as rc

logger = logging.getLogger(os.path.basename(__file__))

# Colour table of the NCL script (RGB 0-256)
COLORS = np.array([
    [15.0, 69, 168], [36, 118, 205], [57, 162, 245], [96, 190, 250],
    [131, 212, 253], [146, 230, 253], [161, 241, 255], [188, 246, 255],
    [205, 226, 229], [234, 231, 211], [251, 246, 190], [255, 232, 154],
    [252, 224, 97], [254, 173, 26], [251, 136, 10], [238, 91, 12],
    [209, 49, 7], [178, 0, 0]]) / 256.0


def drake_passage_section(cube):
    """Section at the longitude closest to 69W: (section, lat, lev, a)."""
    lon = rc.lon_1d(cube)
    if np.max(lon) > 291.0:
        a = rc.closest_index(291.0, lon)
    else:
        a = rc.closest_index(-69.0, lon)
    lat = rc.lat_1d(cube)
    lev = cube.coord(axis="Z").points
    section = np.ma.masked_invalid(cube.data[:, :, a])
    return section, lat, lev, a, lon


def total_transport(section, vol_cube, lat, lon, a, lev):
    """Net transport (Sv) through the passage from volcello."""
    volcello = np.ma.masked_invalid(vol_cube.data)
    vol_lev = vol_cube.coord(axis="Z").points
    dlon = abs(lon[a] - lon[a + 1]) * rc.EARTH_RADIUS * rc.DEG2RAD
    area = volcello / dlon
    if abs(vol_lev[0]) > abs(vol_lev[1]):
        area_2d = area[::-1, :, a]
    else:
        area_2d = area[:, :, a]
    if area_2d.shape != section.shape:
        vol_lat = rc.lat_1d(vol_cube)
        interp = np.ma.masked_all(section.shape)
        filled = area_2d.filled(np.nan)
        for k in range(min(area_2d.shape[0], section.shape[0])):
            interp[k] = np.ma.masked_invalid(
                np.interp(lat, vol_lat, filled[k]))
        area_2d = interp
    transport_per_cell = section * area_2d / 1.0e6  # m3/s -> Sv
    per_lat = transport_per_cell.sum(axis=0)
    per_lat = per_lat / np.cos(lat * rc.DEG2RAD)
    b1 = rc.closest_index(-76.0, lat)
    b2 = rc.closest_index(-48.0, lat)
    return float(per_lat[b1:b2 + 1].sum())


def main(cfg):
    """Run the diagnostic."""
    metadata = list(cfg["input_data"].values())
    uo_data = rc.select_by_short_name(metadata, "uo")
    vol_data = rc.select_by_short_name(metadata, "volcello")
    years = rc.year_range_str(uo_data)

    nvert = int(cfg.get("max_vert", 1))
    nhori = int(cfg.get("max_hori", 1))
    per_page = nvert * nhori
    factor = float(cfg.get("unitCorrectionalFactor", 1.0))
    units_str = f" ({cfg.get('new_units', 'm/s')})"

    levels = np.arange(-40.0, 40.1, 5.0)
    cmap = ListedColormap(COLORS)
    norm = BoundaryNorm(levels, ncolors=cmap.N, extend="both")

    ancestors, nc_files, plot_files = [], [], []
    pages = [uo_data[i:i + per_page]
             for i in range(0, len(uo_data), per_page)]
    for ipage, page in enumerate(pages):
        fig, axs = plt.subplots(nvert, nhori,
                                figsize=(8 * nhori, 5 * nvert),
                                squeeze=False)
        for iax in range(per_page):
            axes = axs.flat[iax]
            if iax >= len(page):
                axes.set_visible(False)
                continue
            meta = page[iax]
            cube = rc.time_mean(
                rc.load_cube(meta["filename"], meta["short_name"]))
            section, lat, lev, a, lon = drake_passage_section(cube)

            vol_meta = rc.match_by_dataset(vol_data, meta["dataset"])
            transport = None
            if vol_meta is not None:
                try:
                    vol_cube = rc.load_cube(vol_meta["filename"],
                                            "volcello")
                    transport = total_transport(
                        section, vol_cube, lat, lon, a, lev)
                except Exception:
                    logger.warning(
                        "Transport calculation failed for %s",
                        meta["dataset"], exc_info=True)
            else:
                logger.warning(
                    "volcello file for %s not found in the recipe; "
                    "skipping the transport calculation.",
                    meta["dataset"])

            exact_lon = lon[a]
            exact_lon = (360.0 - exact_lon if exact_lon > 200.0
                         else -exact_lon)

            b1 = rc.closest_index(-76.0, lat)
            b2 = rc.closest_index(-48.0, lat)
            plot_section = section[:, b1:b2 + 1] * factor
            plot_lat = lat[b1:b2 + 1]

            filled = axes.contourf(plot_lat, lev, plot_section,
                                   levels=levels, cmap=cmap, norm=norm,
                                   extend="both")
            axes.contour(plot_lat, lev, plot_section, levels=levels,
                         colors="black", linewidths=0.3)
            axes.set_facecolor("dimgrey")
            axes.set_xlim(-74, -50)
            axes.set_xticks(np.arange(-74, -49, 2))
            axes.invert_yaxis()
            axes.set_xlabel("Latitude")
            axes.set_ylabel("Depth (m)")
            axes.set_title(
                f"Section velocity of {meta['dataset']}{units_str}",
                fontsize=11)
            left = (f"Net transport : {transport:4.1f}Sv"
                    if transport is not None else " no volcello file ")
            axes.text(0.0, 1.02, left, transform=axes.transAxes,
                      fontsize=8)
            axes.text(1.0, 1.02,
                      f"Drake passage ({exact_lon:4.2f} W)",
                      ha="right", transform=axes.transAxes, fontsize=8)
            fig.colorbar(filled, ax=axes, orientation="vertical",
                         shrink=0.9)

            nc_name = get_diagnostic_filename(
                f"russell18jgr_fig4_uo_{meta['dataset']}_"
                f"{meta['start_year']}-{meta['end_year']}", cfg)
            iris.save(cube, nc_name)
            nc_files.append(nc_name)
            ancestors.append(meta["filename"])
        suffix = f"_page{ipage + 1}" if len(pages) > 1 else ""
        fig.tight_layout()
        plot_file = get_plot_filename(
            f"Russell18jgr-fig4_{years}{suffix}", cfg)
        fig.savefig(plot_file, bbox_inches="tight", dpi=200)
        plt.close(fig)
        plot_files.append(plot_file)
        logger.info("Wrote %s", plot_file)

    record = rc.provenance_record(
        "Russell et al 2018 figure 4", ancestors)
    with ProvenanceLogger(cfg) as prov:
        for filename in nc_files + plot_files:
            prov.log(filename, record)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
