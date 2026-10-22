"""Heave / along-isopycnal decomposition (COSIMA
Heaving_Water_mass_transformation_decomposition).

Translation of
``03-Advanced-Recipes/Heaving_Water_mass_transformation_decomposition.ipynb``
onto CMORised data.  A temperature anomaly between two states has two
causes that the raw difference conflates:

- *along-isopycnal* change, the water on a given density surface
  actually warming or cooling, which is water-mass transformation;
- *heave*, the density surfaces moving up or down while the water on
  them is unchanged.

The decomposition follows the notebook exactly:

1. remap temperature from depth onto density coordinates for both
   states, ``T(z) -> T(sigma)``;
2. difference in density space; that is the along-isopycnal part;
3. remap that difference back onto depth coordinates;
4. the residual, total anomaly minus the remapped along-isopycnal part,
   is the heave.

The notebook contrasts a control and a perturbation experiment.  Here
the two states are two ``exp`` facets of the same dataset (for example
``historical`` against ``piControl``), which the recipe supplies as two
variable groups.

``xgcm``'s ``transform`` does the remapping upstream; here it is a
per-column linear interpolation onto the density bins, which is the
same operation for a monotonic density profile and is the honest
fallback where the profile inverts (those columns are masked).
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

#: The notebook's two-part bins, coarse in light water and fine in dense.
DEFAULT_DENSITY_BINS = np.concatenate(
    [np.arange(20.0, 25.0, 0.1), np.arange(25.0, 38.0, 0.02)])


def remap_to_density(field, density, bins):
    """Interpolate a ``(depth, ...)`` field onto density bins.

    Columns whose density profile is not monotonically increasing are
    masked rather than silently reordered; a non-monotonic profile means
    the column is statically unstable and the remapping is undefined.
    """
    field = np.ma.masked_invalid(np.ma.asarray(field, dtype=float))
    density = np.ma.masked_invalid(np.ma.asarray(density, dtype=float))
    centres = cc.bin_centres(bins)

    flat_field = field.reshape(field.shape[0], -1)
    flat_density = density.reshape(density.shape[0], -1)
    output = np.ma.masked_all((centres.size, flat_field.shape[1]))

    for column in range(flat_field.shape[1]):
        values = flat_field[:, column]
        coordinate = flat_density[:, column]
        valid = ~np.ma.getmaskarray(values) & ~np.ma.getmaskarray(coordinate)
        if valid.sum() < 2:
            continue
        rho = np.ma.filled(coordinate, np.nan)[valid]
        data = np.ma.filled(values, np.nan)[valid]
        if not np.all(np.diff(rho) > 0):
            continue
        inside = (centres >= rho[0]) & (centres <= rho[-1])
        if not inside.any():
            continue
        output[inside, column] = np.interp(centres[inside], rho, data)
    return output.reshape((centres.size,) + field.shape[1:])


def remap_to_depth(field_rho, density, bins):
    """Interpolate a density-coordinate field back onto depth levels."""
    field_rho = np.ma.masked_invalid(np.ma.asarray(field_rho, dtype=float))
    density = np.ma.masked_invalid(np.ma.asarray(density, dtype=float))
    centres = cc.bin_centres(bins)

    flat_field = field_rho.reshape(field_rho.shape[0], -1)
    flat_density = density.reshape(density.shape[0], -1)
    output = np.ma.masked_all(density.shape).reshape(density.shape[0], -1)

    for column in range(flat_density.shape[1]):
        target = flat_density[:, column]
        values = flat_field[:, column]
        valid_source = ~np.ma.getmaskarray(values)
        valid_target = ~np.ma.getmaskarray(target)
        if valid_source.sum() < 2 or not valid_target.any():
            continue
        rho = centres[valid_source]
        data = np.ma.filled(values, np.nan)[valid_source]
        wanted = np.ma.filled(target, np.nan)
        inside = valid_target & (wanted >= rho[0]) & (wanted <= rho[-1])
        if not inside.any():
            continue
        output[inside, column] = np.interp(wanted[inside], rho, data)
    return output.reshape(density.shape)


def decompose(thetao_a, so_a, thetao_b, so_b, bins, p_ref=1000):
    """Total, along-isopycnal and heave temperature anomalies.

    ``*_b`` is the reference state; the anomaly is ``b - a`` reversed,
    that is ``perturbation - control``.
    """
    sigma_a = cc.potential_density_anomaly(thetao_a, so_a, p_ref=p_ref)
    sigma_b = cc.potential_density_anomaly(thetao_b, so_b, p_ref=p_ref)

    temp_a = cc.masked_data(thetao_a)
    temp_b = cc.masked_data(thetao_b)

    total = temp_b - temp_a

    temp_a_rho = remap_to_density(temp_a, sigma_a, bins)
    temp_b_rho = remap_to_density(temp_b, sigma_b, bins)
    along_rho = temp_b_rho - temp_a_rho

    # Bring the along-isopycnal anomaly back onto the control's depths.
    along = remap_to_depth(along_rho, sigma_a, bins)
    heave = total - along
    return total, along, heave


def main(cfg):
    """Run the diagnostic."""
    bins = np.asarray(cfg.get("density_bins", DEFAULT_DENSITY_BINS),
                      dtype=float)
    p_ref = int(cfg.get("reference_pressure", 1000))
    reference_exp = cfg.get("reference_experiment", "piControl")
    grouped = group_metadata(cfg["input_data"].values(), "dataset")

    for dataset, metadata in sorted(grouped.items()):
        by_experiment = group_metadata(metadata, "exp")
        if reference_exp not in by_experiment or len(by_experiment) < 2:
            logger.warning(
                "%s: need two experiments including %s, found %s",
                dataset, reference_exp, ", ".join(sorted(by_experiment)))
            continue
        perturbation_exp = next(
            name for name in sorted(by_experiment) if name != reference_exp)

        cubes = {}
        for exp in (reference_exp, perturbation_exp):
            by_variable = group_metadata(by_experiment[exp], "short_name")
            if "thetao" not in by_variable or "so" not in by_variable:
                cubes = {}
                break
            cubes[exp] = (
                cc.time_mean(cc.load_cube(
                    by_variable["thetao"][0]["filename"], "thetao")),
                cc.time_mean(cc.load_cube(
                    by_variable["so"][0]["filename"], "so")))
        if not cubes:
            logger.warning("%s: thetao and so needed for both experiments",
                           dataset)
            continue

        total, along, heave = decompose(
            *cubes[reference_exp], *cubes[perturbation_exp], bins,
            p_ref=p_ref)

        reference_cube = cubes[reference_exp][0]
        depth_name = cc.depth_coord_name(reference_cube)
        depth = reference_cube.coord(depth_name).points
        lat = cc.lat_1d(reference_cube)

        # Zonal mean sections, which is how the notebook presents it.
        # The three fields must be averaged over the same columns, or
        # the decomposition stops closing in the mean: `along` is masked
        # wherever the remapping failed, and averaging `total` over the
        # wider set of columns leaves a residual that looks like heave
        # but is a bookkeeping artefact.
        common_mask = (np.ma.getmaskarray(total)
                       | np.ma.getmaskarray(along)
                       | np.ma.getmaskarray(heave))
        sections = [np.ma.mean(np.ma.masked_where(common_mask, field),
                               axis=-1)
                    for field in (total, along, heave)]
        titles = ["Total anomaly", "Along-isopycnal (water-mass)", "Heave"]

        figure, axes_list = plt.subplots(1, 3, figsize=(18, 5))
        scale = max(
            float(np.ma.max(np.ma.abs(section))) for section in sections
            if section.count()) or 1.0
        for axes, section, title in zip(axes_list, sections, titles):
            mesh = axes.pcolormesh(lat, depth, section, cmap="RdBu_r",
                                   shading="auto", vmin=-scale, vmax=scale)
            axes.invert_yaxis()
            axes.set_title(title, fontsize=11)
            axes.set_xlabel("Latitude")
            axes.set_ylabel("Depth (m)")
            figure.colorbar(mesh, ax=axes, shrink=0.8, pad=0.03,
                            label="$^\\circ$C")
        figure.suptitle(
            f"{dataset}: {perturbation_exp} - {reference_exp} temperature "
            f"anomaly decomposition (sigma{p_ref // 1000})", fontsize=12)

        plot_path = get_plot_filename(
            f"cosima_heaving_decomposition_{dataset}", cfg)
        figure.savefig(plot_path, bbox_inches="tight", dpi=150)
        plt.close(figure)
        logger.info("Wrote %s", plot_path)

        coords = [(depth, "depth", "m"), (lat, "latitude", "degrees_north")]
        out = iris.cube.CubeList([
            cc.make_cube(np.ma.filled(section, np.nan), coords, var_name,
                         "degC", long_name)
            for section, var_name, long_name in zip(
                sections, ("thetao_total", "thetao_along", "thetao_heave"),
                ("total temperature anomaly",
                 "along-isopycnal temperature anomaly",
                 "temperature anomaly due to isopycnal heave"))])
        path = get_diagnostic_filename(
            f"cosima_heaving_decomposition_{dataset}", cfg)
        iris.save(out, path)
        record = {
            "caption": f"Decomposition of the {dataset} "
                       f"{perturbation_exp} minus {reference_exp} "
                       f"temperature anomaly into along-isopycnal and heave "
                       f"components.",
            "statistics": ["anomaly"],
            "domains": ["global"],
            "plot_types": ["zonal"],
            "authors": ["beucher_romain"],
            "references": ["cosima_recipes"],
            "ancestors": [m["filename"] for m in metadata],
        }
        with ProvenanceLogger(cfg) as logger_:
            logger_.log(path, record)
            logger_.log(plot_path, record)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
