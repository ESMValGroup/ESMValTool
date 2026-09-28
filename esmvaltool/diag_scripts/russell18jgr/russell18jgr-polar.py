"""Polar contour plot with land masked (russell18jgr figures 1, 7a, 8a).

Python port of russell18jgr-polar.ncl.

- Uses the original grid (no regridding).
- Expects time-averaged (and for tauu/fgco2 land-masked) input from the
  ESMValTool preprocessor.
- Plots a Southern-Hemisphere polar contour map per dataset and panels
  them max_vert x max_hori per page.

Diagnostic script options (same as the NCL version):
  max_lat               : equatorward plot limit (e.g. -30.)
  max_vert, max_hori    : panel layout per page
  grid_min/max/step     : contour levels
  colors                : list of RGB triples (0-255) used as palette
  colormap              : name of a colormap (NCL 'BlWhRe' -> 'bwr')
  labelBar_end_type     : 'ExcludeOuterBoxes' or triangle ends
  unitCorrectionalFactor: multiplicative unit conversion
  new_units             : unit string shown in the title
"""

import logging
import os
import sys

import cartopy.crs as ccrs
import iris
import matplotlib

matplotlib.use("Agg")
import matplotlib.path as mpath
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


def colormap_from_cfg(cfg, nlevels):
    """Build the colormap requested in the recipe."""
    if "colors" in cfg:
        rgb = np.array(cfg["colors"], dtype=float) / 256.0
        return ListedColormap(np.clip(rgb, 0, 1))
    name = cfg.get("colormap", "nrl_sirkes_nowhite")
    ncl_to_mpl = {"BlWhRe": "bwr", "nrl_sirkes_nowhite": "RdYlBu_r"}
    cmap = plt.get_cmap(ncl_to_mpl.get(name, name))
    return ListedColormap([cmap(i / max(nlevels, 1))
                           for i in range(nlevels + 1)])


def circular_boundary():
    """Circle in axes coordinates to clip the polar plot."""
    theta = np.linspace(0, 2 * np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    return mpath.Path(verts * radius + center)


def plot_dataset(axes, cube, levels, cmap, extend, cfg, title, subtitle,
                 units_str, var_name):
    """Draw one polar contour panel."""
    lat = rc.coord_1d_or_2d(cube, "latitude")
    lon = rc.coord_1d_or_2d(cube, "longitude")
    if lat.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon, lat)
    else:
        lat2d, lon2d = lat, lon

    max_lat = float(cfg.get("max_lat", 0.0))
    data = np.ma.masked_invalid(cube.data)
    # The cubes are global but the map shows the Southern Ocean only.
    # Without this, contours of northern hemisphere features (the
    # Baltic, the Great Lakes, the Kara Sea) are computed and can place
    # labels on the plot.
    data = np.ma.masked_where(lat2d > max_lat, data)
    logger.info(
        "%s: %s range %.4g to %.4g south of %g degrees latitude",
        title, var_name, data.min(), data.max(), max_lat)

    axes.set_extent([-180, 180, -90, max_lat], ccrs.PlateCarree())
    axes.set_boundary(circular_boundary(), transform=axes.transAxes)

    norm = BoundaryNorm(levels, ncolors=cmap.N, extend=extend)
    filled = axes.contourf(
        lon2d, lat2d, data, levels=levels, cmap=cmap, norm=norm,
        extend=extend, transform=ccrs.PlateCarree())
    # A line at every level is unreadable where the field is steep, so
    # draw at most ~8 of them
    stride = max(1, len(levels) // 8)
    axes.contour(
        lon2d, lat2d, data, levels=levels[::stride], colors="black",
        linewidths=0.3, transform=ccrs.PlateCarree())
    try:
        axes.coastlines(linewidth=0.4)
        axes.add_feature(
            __import__("cartopy.feature", fromlist=["LAND"]).LAND,
            facecolor=(0.5, 0.5, 0.5), zorder=2)
    except Exception:  # noqa: BLE001 - offline cartopy data
        logger.warning("Could not draw coastlines/land feature")
    grid_color = cfg.get("grid_color", "green")
    gridlines = axes.gridlines(color=grid_color, linewidth=0.5,
                               linestyle="-")
    gridlines.ylocator = plt.MultipleLocator(10)
    # title on top, then the right and left strings below it, so that
    # they never overlap (the NCL layout)
    axes.set_title(title, fontsize=11, pad=30)
    axes.text(1.0, 1.06, subtitle, fontsize=8,
              transform=axes.transAxes, ha="right")
    axes.text(0.0, 1.01, f"{var_name}{units_str}", fontsize=8,
              transform=axes.transAxes, ha="left")
    return filled


def main(cfg):
    """Run the diagnostic."""
    input_data = list(cfg["input_data"].values())
    var0 = input_data[0]["short_name"]
    years = rc.year_range_str(input_data)

    grid_min = float(cfg.get("grid_min", 0.0))
    grid_max = float(cfg.get("grid_max", 1.0))
    grid_step = float(cfg.get("grid_step", 0.1))
    nsteps = int(round((grid_max - grid_min) / grid_step)) + 1
    levels = np.linspace(grid_min, grid_max, nsteps)
    extend = ("neither" if cfg.get("labelBar_end_type")
              == "ExcludeOuterBoxes" else "both")
    cmap = colormap_from_cfg(cfg, nsteps)

    nvert = int(cfg.get("max_vert", 1))
    nhori = int(cfg.get("max_hori", 1))
    per_page = nvert * nhori

    factor = float(cfg.get("unitCorrectionalFactor", 1.0))
    ancestors, nc_files = [], []

    pages = [input_data[i:i + per_page]
             for i in range(0, len(input_data), per_page)]
    plot_files = []
    for ipage, page in enumerate(pages):
        fig = plt.figure(figsize=(5.5 * nhori, 5.5 * nvert))
        for ipanel, meta in enumerate(page):
            cube = rc.load_cube(meta["filename"], meta["short_name"])
            cube = rc.time_mean(cube)
            # ph is a full-depth field in CMIP5 and figure 8 is the
            # surface: take the surface here too, so that the saved
            # netCDF matches what is plotted
            cube = rc.surface_field(cube)
            if factor != 1.0:
                cube.data = cube.data * factor
            units_str = f" ({cfg.get('new_units', cube.units)})"
            axes = fig.add_subplot(nvert, nhori, ipanel + 1,
                                   projection=ccrs.SouthPolarStereo())
            subtitle = (f"annual mean {meta['start_year']} - "
                        f"{meta['end_year']}")
            filled = plot_dataset(
                axes, cube, levels, cmap, extend, cfg,
                meta["dataset"], subtitle, units_str, meta["short_name"])
            cbar = fig.colorbar(filled, ax=axes, orientation="vertical",
                                shrink=0.85, ticks=levels)
            cbar.ax.tick_params(labelsize=6)
            nc_name = get_diagnostic_filename(
                f"russell18jgr_polar_{meta['short_name']}_"
                f"{meta['dataset']}_{meta['start_year']}-"
                f"{meta['end_year']}", cfg)
            iris.save(cube, nc_name)
            nc_files.append(nc_name)
            ancestors.append(meta["filename"])
        suffix = f"_page{ipage + 1}" if len(pages) > 1 else ""
        plot_file = get_plot_filename(
            f"Russell_polar-contour_{var0}_{years}{suffix}", cfg)
        fig.savefig(plot_file, bbox_inches="tight", dpi=200)
        plt.close(fig)
        plot_files.append(plot_file)
        logger.info("Wrote %s", plot_file)

    record = rc.provenance_record(
        f"Russell et al 2018 polar plot {var0}", ancestors)
    with ProvenanceLogger(cfg) as prov:
        for filename in nc_files + plot_files:
            prov.log(filename, record)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
