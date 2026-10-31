"""Along-isobath averaging (COSIMA Along_Isobath_Average).

Translation of ``03-Advanced-Recipes/Along_Isobath_Average.ipynb`` onto
CMORised data.  Antarctic shelf analyses want an average taken along
bathymetric contours rather than along latitude circles, because the
shelf break does not follow a latitude circle.  The notebook bins a
scalar field by bottom depth with ``xhistogram``, weighting by cell
thickness and area; here the binning is
:func:`cosima_common.bin_sum`, which is the same operation without the
xarray dependency.

The notebook's two caveats carry over unchanged and are worth
repeating, because they are properties of the method rather than of the
implementation:

1. deep troughs cutting into the shelf are averaged together with the
   open ocean at the same depth, and offshore seamounts with the slope;
   the region mask has to be chosen with that in mind;
2. MOM's ``ht``/``hu`` bottom depths do not line up with the model's
   own vertical levels, so the notebook masks with the bathymetry but
   bins by the bottom depth implied by the 3D field.  The same is done
   here: the mask comes from ``deptho`` and the binning coordinate is
   the depth of the deepest wet cell of the field itself.

The result is a depth-binned profile of the scalar, over the region the
recipe selects.
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

DEFAULT_ISOBATH_BINS = np.arange(0.0, 4001.0, 100.0)


def column_bottom_depth(cube):
    """Depth of the deepest wet cell of each column.

    This is the notebook's ``ht_coarse``: the bathymetry as the 3D field
    itself sees it, discretised onto the model's own levels, rather than
    the continuous ``deptho``.
    """
    depth_name = cc.depth_coord_name(cube)
    if depth_name is None:
        raise ValueError("a 3D field is needed to derive the bottom depth")
    levels = cube.coord(depth_name).points
    marker = cube.copy(
        np.ma.masked_array(
            np.broadcast_to(levels[:, None, None],
                            cube.shape[-3:]).astype(float).copy(),
            mask=np.ma.getmaskarray(
                cc.masked_data(cube))))
    return cc.bottom_cell(marker)


def isobath_average(field, bottom_depth, bins, area, mask=None):
    """Area-weighted mean of ``field`` in each bottom-depth bin."""
    if mask is not None:
        field = np.ma.masked_where(~mask, field)
        bottom_depth = np.ma.masked_where(~mask, bottom_depth)
    weights = np.ma.masked_where(np.ma.getmaskarray(field), area)
    weighted = cc.bin_sum(field * weights, bottom_depth, bins)
    totals = cc.bin_sum(weights, bottom_depth, bins)
    with np.errstate(divide="ignore", invalid="ignore"):
        profile = np.where(totals > 0, weighted / totals, np.nan)
    return np.ma.masked_invalid(profile)


def main(cfg):
    """Run the diagnostic."""
    bins = np.asarray(cfg.get("isobath_bins", DEFAULT_ISOBATH_BINS),
                      dtype=float)
    mask_depth = float(cfg.get("mask_depth", 2000.0))
    variables = cfg.get("variables", ["thetao", "so"])
    grouped = group_metadata(cfg["input_data"].values(), "dataset")
    centres = cc.bin_centres(bins)

    figure, axes_list = plt.subplots(
        1, max(len(variables), 1), figsize=(6 * max(len(variables), 1), 5.5),
        squeeze=False)
    plotted = 0

    for index, (dataset, metadata) in enumerate(sorted(grouped.items())):
        by_variable = group_metadata(metadata, "short_name")
        if "deptho" not in by_variable:
            logger.warning("%s: no deptho, skipping", dataset)
            continue
        deptho = cc.load_cube(by_variable["deptho"][0]["filename"], "deptho")
        bathymetry = cc.masked_data(deptho)
        region_mask = bathymetry < mask_depth
        style = cc.style_for(dataset, index)

        for axes, name in zip(axes_list[0], variables):
            if name not in by_variable:
                logger.warning("%s: %s not available", dataset, name)
                continue
            cube = cc.time_mean(
                cc.load_cube(by_variable[name][0]["filename"], name))
            try:
                bottom_depth = column_bottom_depth(cube)
            except ValueError as error:
                logger.warning("%s/%s: %s", dataset, name, error)
                continue
            surface = cc.surface_field(cube)
            field = cc.masked_data(surface)
            area = cc.cell_area(surface)

            profile = isobath_average(field, bottom_depth, bins, area,
                                      mask=region_mask)
            axes.plot(profile, centres, color=style["color"],
                      linestyle=style["dash"], linewidth=2.0, marker="o",
                      markersize=3, label=dataset)
            axes.set_title(f"{name} along isobaths", fontsize=11)
            axes.set_xlabel(f"{name} ({cube.units})")
            plotted += 1

            out = cc.make_cube(
                np.ma.filled(profile, np.nan),
                [(centres, "depth", "m")], name, str(cube.units),
                f"{name} averaged along isobaths")
            path = get_diagnostic_filename(
                f"cosima_along_isobath_{name}_{dataset}", cfg)
            iris.save(out, path)
            record = {
                "caption": f"{name} of {dataset} averaged along isobaths "
                           f"over the region shallower than "
                           f"{mask_depth:.0f} m.",
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
        logger.warning("Nothing to average along isobaths")
        plt.close(figure)
        return

    for axes in axes_list[0]:
        axes.invert_yaxis()
        axes.set_ylabel("Bottom depth (m)")
        axes.grid(alpha=0.3)
        axes.legend(fontsize=8, frameon=False)

    plot_path = get_plot_filename("cosima_along_isobath_average", cfg)
    figure.savefig(plot_path, bbox_inches="tight", dpi=150)
    plt.close(figure)
    with ProvenanceLogger(cfg) as logger_:
        logger_.log(plot_path, {
            **record,
            "caption": "Model hydrographic profiles grouped by bottom-depth bins.",
            "ancestors": [m["filename"] for m in cfg["input_data"].values()],
        })
    logger.info("Wrote %s", plot_path)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
