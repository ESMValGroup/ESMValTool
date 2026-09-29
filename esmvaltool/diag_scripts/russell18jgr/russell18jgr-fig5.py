"""Sea-ice max/min extent polar plot (russell18jgr figure 5).

Python port of russell18jgr-fig5.ncl.

- Builds the monthly climatology of sic.
- Marks cells with September (max) coverage > 15% in blue and cells
  with March (min) coverage > 15% in red (March wins where both).
- Plots one SH polar panel per dataset, paged max_vert x max_hori.
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


def circular_boundary():
    """Circle in axes coordinates to clip the polar plot."""
    theta = np.linspace(0, 2 * np.pi, 100)
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    return mpath.Path(verts * 0.5 + [0.5, 0.5])


def ice_extent_field(cube):
    """Return the 1 (September) / 2 (March) masked field."""
    clim = rc.monthly_climatology(cube)
    sic = rc.sea_ice_percent(clim)
    field = np.ma.masked_all(sic[0].shape)
    # the edge of full coverage is defined by 15% areal coverage;
    # March is applied last so it wins where both months qualify
    field[np.ma.filled(sic[8] > 15.0, False)] = 1.0
    field[np.ma.filled(sic[2] > 15.0, False)] = 2.0
    if field.count() == 0:
        logger.warning(
            "No cell exceeds 15%% sea ice concentration - the panel will "
            "be empty (sic range %s to %s)",
            sic.min(),
            sic.max(),
        )
    return field, clim


def main(cfg):
    """Run the diagnostic."""
    input_data = [
        m
        for m in cfg["input_data"].values()
        if m["short_name"] in ("sic", "siconc")
    ]
    var0 = input_data[0]["short_name"]
    years = rc.year_range_str(input_data)

    nvert = int(cfg.get("max_vert", 1))
    nhori = int(cfg.get("max_hori", 1))
    per_page = nvert * nhori
    max_lat = float(cfg.get("max_lat", 0.0))

    cmap = ListedColormap(["#00008b", "#ff0000"])  # blue4, red
    norm = BoundaryNorm([0.5, 1.5, 2.5], ncolors=2)

    ancestors, nc_files, plot_files = [], [], []
    pages = [
        input_data[i : i + per_page]
        for i in range(0, len(input_data), per_page)
    ]
    for ipage, page in enumerate(pages):
        fig = plt.figure(figsize=(5.5 * nhori, 5.5 * nvert))
        for ipanel, meta in enumerate(page):
            cube = rc.load_cube(meta["filename"], meta["short_name"])
            field, clim = ice_extent_field(cube)
            lat = rc.coord_1d_or_2d(clim, "latitude")
            lon = rc.coord_1d_or_2d(clim, "longitude")
            if lat.ndim == 1:
                lon2d, lat2d = np.meshgrid(lon, lat)
            else:
                lat2d, lon2d = lat, lon

            axes = fig.add_subplot(
                nvert, nhori, ipanel + 1, projection=ccrs.SouthPolarStereo()
            )
            axes.set_extent([-180, 180, -90, max_lat], ccrs.PlateCarree())
            axes.set_boundary(circular_boundary(), transform=axes.transAxes)
            axes.pcolormesh(
                lon2d,
                lat2d,
                field,
                cmap=cmap,
                norm=norm,
                transform=ccrs.PlateCarree(),
                shading="auto",
            )
            try:
                import cartopy.feature as cfeature

                axes.add_feature(
                    cfeature.LAND, facecolor=(0.5, 0.5, 0.5), zorder=2
                )
                axes.coastlines(linewidth=0.4, zorder=3)
            except Exception:  # noqa: BLE001
                logger.warning("Could not draw coastlines/land feature")
            axes.gridlines(color="green", linewidth=0.5)
            axes.set_title(meta["dataset"], fontsize=10)
            axes.text(
                0.5,
                1.06,
                "(Blue - September & Red - March)",
                fontsize=8,
                ha="center",
                transform=axes.transAxes,
            )
            axes.text(
                0.0,
                1.12,
                "Southern Ocean Max Min Sea ice extent",
                fontsize=8,
                ha="left",
                transform=axes.transAxes,
            )
            axes.text(
                1.0,
                1.12,
                f"annual mean {meta['start_year']} - {meta['end_year']}",
                fontsize=7,
                ha="right",
                transform=axes.transAxes,
            )

            out_cube = iris.cube.Cube(
                field,
                var_name=var0,
                long_name="sea ice extent (1 September / 2 March)",
            )
            for coord in ("latitude", "longitude"):
                dims = clim.coord_dims(clim.coord(coord))
                out_cube.add_aux_coord(
                    clim.coord(coord), tuple(d - 1 for d in dims)
                )
            nc_name = get_diagnostic_filename(
                f"russell18jgr_fig5_{var0}_{meta['dataset']}_"
                f"{meta['start_year']}-{meta['end_year']}",
                cfg,
            )
            iris.save(out_cube, nc_name)
            nc_files.append(nc_name)
            ancestors.append(meta["filename"])
        suffix = f"_page{ipage + 1}" if len(pages) > 1 else ""
        plot_file = get_plot_filename(
            f"russell18jgr-fig5_{var0}_{years}{suffix}", cfg
        )
        fig.savefig(plot_file, bbox_inches="tight", dpi=200)
        plt.close(fig)
        plot_files.append(plot_file)
        logger.info("Wrote %s", plot_file)

    record = rc.provenance_record(
        "Russell et al 2018 figure 5 -polar", ancestors
    )
    with ProvenanceLogger(cfg) as prov:
        for filename in nc_files + plot_files:
            prov.log(filename, record)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
