"""Sea-ice area, extent and volume against observations (COSIMA
Sea_Ice_Area_Concentration_Volume_with_Obs).

Translation of
``03-Advanced-Recipes/Sea_Ice_Area_Concentration_Volume_with_Obs.ipynb``
onto CMORised data.  The notebook compares ACCESS-OM2 CICE output to
the NSIDC climate data record pulled live from an ERDDAP server; here
the observations come through ESMValTool's own CMORisers
(NSIDC-G02202-sh, OSI-450 for ``siconc``, PIOMAS for ``sivol``) and are
therefore just further datasets in the recipe.

Three things are computed per dataset and hemisphere:

- **area**, the concentration-weighted total, sum(siconc * areacello);
- **extent**, the area of every cell at or above 15% concentration;
- **volume**, sum(sivol * areacello), where ``sivol`` is the CMOR
  volume per unit area.

The seasonal cycle of each is plotted, with observations drawn in black.

The percent-versus-fraction trap is worth naming: CMOR ``siconc`` is a
percentage, but several models publish a fraction, and deciding between
them from ``max(siconc) < 5`` fails whenever a single fill value
survives unmasked.  :func:`cosima_common.fraction_from_percent` masks
fill values first and trusts the units when they are meaningful.
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

EXTENT_THRESHOLD = 0.15
M2_TO_MILLION_KM2 = 1.0e-12
M3_TO_THOUSAND_KM3 = 1.0e-12


def _month_numbers(cube):
    """Month number of each time step."""
    coord = cube.coord("time")
    return np.array([cell.point.month for cell in coord.cells()])


def seasonal_cycle(cube, area, hemisphere="south", weighted=True,
                   as_volume=False):
    """Monthly climatology of sea-ice area, extent or volume.

    Returns a 12-element array; months with no data are masked.
    """
    lat = cc.lat_2d(cube)
    in_hemisphere = lat < 0 if hemisphere == "south" else lat > 0

    if as_volume:
        values = cc.masked_data(cube)
    else:
        values = cc.fraction_from_percent(cube)
    if values.ndim == 2:
        values = values[None, ...]

    masked_area = np.ma.masked_where(~in_hemisphere, np.asarray(area))
    months = _month_numbers(cube)

    totals = np.ma.masked_all(12)
    for month in range(1, 13):
        steps = np.where(months == month)[0]
        if steps.size == 0:
            continue
        per_step = [
            cc.seaice_area(values[step], masked_area,
                           threshold=None if as_volume else EXTENT_THRESHOLD,
                           weighted=weighted)
            for step in steps]
        totals[month - 1] = np.mean(per_step)
    return totals


def main(cfg):
    """Run the diagnostic."""
    hemisphere = cfg.get("hemisphere", "south")
    grouped = group_metadata(cfg["input_data"].values(), "dataset")

    figure, axes_list = plt.subplots(1, 3, figsize=(17, 5))
    months = np.arange(1, 13)
    results = {}
    have_data = False

    for index, (dataset, metadata) in enumerate(sorted(grouped.items())):
        by_variable = group_metadata(metadata, "short_name")
        is_obs = any(m.get("project", "").startswith(("OBS", "obs4"))
                     for m in metadata)
        style = cc.style_for(dataset, index)
        colour = "black" if is_obs else style["color"]
        dashes = "-" if not is_obs else "--"
        width = 2.2 if is_obs else 1.6

        if "siconc" in by_variable:
            cube = cc.load_cube(by_variable["siconc"][0]["filename"],
                                "siconc")
            area = cc.cell_area(cube)
            ice_area = seasonal_cycle(cube, area, hemisphere,
                                      weighted=True) * M2_TO_MILLION_KM2
            extent = seasonal_cycle(cube, area, hemisphere,
                                    weighted=False) * M2_TO_MILLION_KM2
            axes_list[0].plot(months, ice_area, color=colour,
                              linestyle=dashes, linewidth=width,
                              label=dataset)
            axes_list[1].plot(months, extent, color=colour,
                              linestyle=dashes, linewidth=width,
                              label=dataset)
            results[(dataset, "siarea")] = ice_area
            results[(dataset, "siextent")] = extent
            have_data = True

        if "sivol" in by_variable:
            cube = cc.load_cube(by_variable["sivol"][0]["filename"], "sivol")
            area = cc.cell_area(cube)
            volume = seasonal_cycle(cube, area, hemisphere,
                                    as_volume=True) * M3_TO_THOUSAND_KM3
            axes_list[2].plot(months, volume, color=colour,
                              linestyle=dashes, linewidth=width,
                              label=dataset)
            results[(dataset, "sivolume")] = volume
            have_data = True

        if not any(key[0] == dataset for key in results):
            logger.warning("%s: neither siconc nor sivol available", dataset)

    if not have_data:
        logger.warning("No sea-ice data to plot")
        plt.close(figure)
        return

    titles = [
        f"Sea-ice area ({hemisphere}ern hemisphere)",
        f"Sea-ice extent (>= {EXTENT_THRESHOLD:.0%})",
        "Sea-ice volume",
    ]
    labels = ["10$^6$ km$^2$", "10$^6$ km$^2$", "10$^3$ km$^3$"]
    for axes, title, label in zip(axes_list, titles, labels):
        axes.set_xticks(months)
        axes.set_xlabel("Month")
        axes.set_ylabel(label)
        axes.set_title(title, fontsize=11)
        axes.grid(alpha=0.3)
        axes.legend(fontsize=8, frameon=False)

    plot_path = get_plot_filename(
        f"cosima_seaice_area_volume_{hemisphere}", cfg)
    figure.savefig(plot_path, bbox_inches="tight", dpi=150)
    plt.close(figure)
    logger.info("Wrote %s", plot_path)

    cubes = iris.cube.CubeList([
        cc.make_cube(np.ma.filled(values, np.nan),
                     [(months, "time", "months")], var_name,
                     "1e6 km2" if var_name != "sivolume" else "1e3 km3",
                     f"{dataset} seasonal cycle of {var_name}")
        for (dataset, var_name), values in sorted(results.items())])
    path = get_diagnostic_filename(
        f"cosima_seaice_area_volume_{hemisphere}", cfg)
    iris.save(cubes, path)
    record = {
        "caption": f"Seasonal cycle of {hemisphere}ern hemisphere sea-ice "
                   f"area, extent and volume, model against observations.",
        "statistics": ["clim"],
        "domains": ["shpolar" if hemisphere == "south" else "nhpolar"],
        "plot_types": ["seas"],
        "authors": ["beucher_romain"],
        "references": ["cosima_recipes"],
        "ancestors": [m["filename"] for m in cfg["input_data"].values()],
    }
    with ProvenanceLogger(cfg) as logger_:
        logger_.log(path, record)
        logger_.log(plot_path, record)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
