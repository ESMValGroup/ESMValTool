"""Cross-contour transport (COSIMA Cross-contour_transport).

Translation of ``03-Advanced-Recipes/Cross-contour_transport.ipynb``
onto CMORised data: the volume transport across a closed bathymetry
contour, for example onto the Antarctic continental shelf across the
1000 m isobath.

The upstream notebook carries a standing warning that its contour and
mask construction can produce contours that are not closed
(COSIMA/cosima-recipes#669, and the non-monotonic indexing of #383).
That is the part this translation changes rather than copies:
:func:`cosima_common.contour_cells` inserts an intermediate cell at
every diagonal step, so the cell sequence is adjacently connected by
construction, and :func:`cosima_common.cross_contour_masks` then
assigns exactly one signed face crossing per step.  The unit checks in
``validation/test_cosima_common.py`` assert both properties, including
that a closed contour has balanced crossings in each direction.

The transport is ``sum(mask_u * umo + mask_v * vmo)`` over the contour,
per depth level, converted from mass to volume transport with the
Boussinesq reference density.  Positive is outward, away from the
shallow side.
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


def contour_from_bathymetry(deptho_cube, contour_depth):
    """Adjacently-connected contour cells at ``contour_depth``."""
    depth = cc.masked_data(deptho_cube)
    x_contour, y_contour = cc.extract_contour(
        np.ma.filled(depth, 0.0), contour_depth)
    if x_contour.size == 0:
        raise ValueError(f"no {contour_depth:.0f} m contour in deptho")
    return cc.contour_cells(x_contour, y_contour, depth.shape)


def cross_contour_transport(umo, vmo, rows, cols):
    """Volume transport per depth level across the contour, in Sv.

    Returns ``(transport, depth)``; ``transport`` is positive outward.
    """
    # umo/vmo are kg s-1 through a cell face, not a flux per unit area.
    # Interpolating them and summing scales the total with the cell
    # count, so a regridded transport is simply wrong.
    cc.check_transport_grid(umo, "umo")
    cc.check_transport_grid(vmo, "vmo")

    u_data = cc.masked_data(umo)
    v_data = cc.masked_data(vmo)
    if u_data.ndim == 4:
        u_data = np.ma.mean(u_data, axis=0)
        v_data = np.ma.mean(v_data, axis=0)

    # ``umo`` and ``vmo`` are native C-grid faces.  Geographic extraction
    # can retain one extra unused edge row from one component, so construct
    # the contour on the shared interior index region.  The helper permits
    # exactly that stagger and rejects any larger/domain-changing mismatch.
    shape = cc.common_face_shape(u_data.shape[-2:], v_data.shape[-2:])
    mask_u, mask_v = cc.cross_contour_masks(rows, cols, shape)
    kind, message = cc.describe_contour(rows, cols, mask_u, mask_v, shape)
    if kind == "open":
        raise ValueError(message)
    logger.info(message)

    face_rows, face_cols, components, signs = cc.cross_contour_faces(
        rows, cols, shape)
    transport_faces = np.ma.masked_all((u_data.shape[0], len(face_rows)))
    for face, (row, col, component, sign) in enumerate(
            zip(face_rows, face_cols, components, signs)):
        values = u_data[:, row, col] if component == "u" else \
            v_data[:, row, col]
        transport_faces[:, face] = sign * values / cc.RHO0 / 1.0e6
    transport = np.ma.sum(transport_faces, axis=-1)

    depth_name = cc.depth_coord_name(umo)
    depth = umo.coord(depth_name).points if depth_name else np.array([0.0])
    return transport, depth


def main(cfg):
    """Run the diagnostic."""
    contour_depth = float(cfg.get("contour_depth", 1000.0))
    grouped = group_metadata(cfg["input_data"].values(), "dataset")

    figure, axes_list = plt.subplots(1, 2, figsize=(14, 5.5))
    plotted = 0

    for index, (dataset, metadata) in enumerate(sorted(grouped.items())):
        by_variable = group_metadata(metadata, "short_name")
        missing = [n for n in ("umo", "vmo", "deptho")
                   if n not in by_variable]
        if missing:
            logger.warning("%s: missing %s, skipping",
                           dataset, ", ".join(missing))
            continue

        deptho = cc.load_cube(by_variable["deptho"][0]["filename"], "deptho")
        umo = cc.load_cube(by_variable["umo"][0]["filename"], "umo")
        vmo = cc.load_cube(by_variable["vmo"][0]["filename"], "vmo")

        # A native C grid can retain one more edge row for one velocity
        # component after geographic extraction. A deptho contour then has
        # the wrong index domain for the shared transport faces. Use deptho
        # when it fits; otherwise derive bottom depth on vmo's native grid.
        # This preserves the extensive transports rather than interpolating.
        try:
            common_shape = cc.common_face_shape(umo.shape[-2:],
                                                vmo.shape[-2:])
            contour_grid_cube = (umo if umo.shape[-2:] == common_shape
                                 else vmo)
            bathymetry, bathymetry_source = cc.bathymetry_on_grid(
                deptho, contour_grid_cube)
        except ValueError as error:
            logger.warning("%s: %s", dataset, error)
            continue
        try:
            x_contour, y_contour = cc.extract_contour(
                np.ma.filled(bathymetry, 0.0), contour_depth)
            if x_contour.size == 0:
                raise ValueError(f"no {contour_depth:.0f} m contour")
            rows, cols = cc.contour_cells(
                x_contour, y_contour, bathymetry.shape)
        except ValueError as error:
            logger.warning("%s: %s", dataset, error)
            continue
        logger.info("%s: contour has %d cells (%s)", dataset, rows.size,
                    bathymetry_source)
        transport, depth = cross_contour_transport(umo, vmo, rows, cols)

        style = cc.style_for(dataset, index)
        axes_list[1].plot(transport, depth, color=style["color"],
                          linestyle=style["dash"], linewidth=2.0,
                          label=f"{dataset} (net "
                                f"{np.ma.sum(transport):+.2f} Sv)")
        plotted += 1

        if plotted == 1:
            # Plot against the coordinate system on which the contour was
            # calculated. It is deptho when grids agree, otherwise vmo.
            bathymetry_cube = (deptho if bathymetry_source == "deptho"
                               else contour_grid_cube)
            # 2D coordinates, since the cube is on its native
            # (curvilinear) grid: a representative 1D latitude would
            # put the contour in the wrong place on the map.
            # Native tripolar coordinate centres are non-monotonic at the
            # fold.  ``pcolormesh`` would infer spurious cell corners there;
            # draw the actual centres instead of inventing edges.
            map_lat, map_lon = (cc.lat_2d(bathymetry_cube),
                                cc.lon_2d(bathymetry_cube))
            valid = ~np.ma.getmaskarray(bathymetry)
            axes_list[0].scatter(
                map_lon[valid], map_lat[valid], c=bathymetry[valid],
                s=4, marker="s", linewidths=0, cmap="Blues",
                rasterized=True)
            contour_lat, contour_lon = cc.contour_indices_to_coordinates(
                rows, cols, cc.lat_2d(bathymetry_cube),
                cc.lon_2d(bathymetry_cube))
            axes_list[0].plot(contour_lon, contour_lat, color="red",
                              linewidth=1.0)
            axes_list[0].set_title(
                f"{dataset}: {contour_depth:.0f} m contour "
                f"({rows.size} cells; {bathymetry_source})", fontsize=11)
            axes_list[0].set_xlabel("Longitude")
            axes_list[0].set_ylabel("Latitude")

        out = cc.make_cube(
            np.ma.filled(transport, np.nan), [(depth, "depth", "m")],
            "cross_contour_transport", "Sv",
            f"volume transport across the {contour_depth:.0f} m isobath")
        path = get_diagnostic_filename(
            f"cosima_cross_contour_transport_{dataset}", cfg)
        iris.save(out, path)
        record = {
            "caption": f"Volume transport of {dataset} across the "
                       f"{contour_depth:.0f} m isobath, by depth. Positive "
                       f"is outward, away from the shallow side.",
            "statistics": ["mean"],
            "domains": ["shpolar"],
            "plot_types": ["vert"],
            "authors": ["beucher_romain"],
            "references": ["cosima_recipes"],
            "ancestors": [m["filename"] for m in metadata],
        }
        with ProvenanceLogger(cfg) as logger_:
            logger_.log(path, record)

    if not plotted:
        logger.warning("No dataset had umo, vmo and deptho")
        plt.close(figure)
        return

    axes_list[1].axvline(0.0, color="black", linewidth=0.6)
    axes_list[1].invert_yaxis()
    axes_list[1].set_xlabel("Cross-contour transport (Sv), positive outward")
    axes_list[1].set_ylabel("Depth (m)")
    axes_list[1].set_title("Transport across the isobath", fontsize=11)
    axes_list[1].legend(fontsize=8, frameon=False)
    axes_list[1].grid(alpha=0.3)

    plot_path = get_plot_filename("cosima_cross_contour_transport", cfg)
    figure.savefig(plot_path, bbox_inches="tight", dpi=150)
    plt.close(figure)
    with ProvenanceLogger(cfg) as logger_:
        logger_.log(plot_path, {
            **record,
            "caption": "Bathymetry contour and cross-contour transport by model.",
            "ancestors": [m["filename"] for m in cfg["input_data"].values()],
        })
    logger.info("Wrote %s", plot_path)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
