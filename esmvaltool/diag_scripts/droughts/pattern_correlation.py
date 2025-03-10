"""Overview plots for pattern correlations of multiple variables.
=================================================================

This diagnostic calculates the pattern correlation between pairs of datasets
over several datasets and plots them for all variables and datasets in one
overview plot. In each figure the correlations are grouped by variable and
project. Only projects listed in the recipe are used everything else is added
with individual markers and labels if `extra_datasets` is True.
A single reference dataset needs to be specified in the recipe.

The diagnostic requires 2d cubes as input. They can be prepared using the
preprocessor or or by any ancestor diagnostic, that provide a metadata.yml.
The additional `group_by` option allows to apply this diagnostic to each entry
of a meta key seperatly including extra facets. For example this can be used to
plot different aggregations of the same variable, like mean and trend in a
single run.

For example:
.. code-block:: yaml
    script: droughtindex/pattern_correlation.py
    reference: ERA5
    group_by: diffmap_metric  # first, last, diff, total
    ancestors: [obs/diffmaps, models/diffmaps]


Configuration options in recipe
-------------------------------
reference: str, required
    Dataset name used to correlate all other datasets against.
group_by: str, optional (default: None)
    Plot figures for each entry of the given key.
project_plot_kwargs: dict, optional (default: {})
    Kwargs passed to the patches for specific projects the plot function.
plot_kwargs: dict, optional (default: {})
    Kwargs passed to the plot function.
projects: list, optional (default: ["CMIP6"])
    List of projects to include as individual colors in the plot.
extra_datasets: bool, optional (default: True)
    Add datasets not belonging to any of the projects as individual markers.
relative_change: bool, optional (default: False)
    Creates an additional plot with global relative changes fo each variable
    in all datasets.
labels: dict, optional
    Mapping of variable names to custom xtick labels. Falls back to variable
    name if not present.
"""

# from esmvalcore import preprocessor as pp
import logging
from collections import defaultdict
from pprint import PrettyPrinter

import iris
import iris.analysis
import iris.analysis.cartography
import matplotlib.pyplot as plt
from iris.analysis import MEAN
from iris.analysis.stats import pearsonr

from esmvaltool.diag_scripts import shared
from esmvaltool.diag_scripts.droughtindex import utils as ut

logger = logging.getLogger(__file__)
p = PrettyPrinter(indent=4)


def pattern_correlation(cube1, cube2, centered=False, weighted=False):
    """Calculate pattern correlation between two 2d-cubes.
    uses pearsonr from scipy.stats to calculate the correlation coefficient
    along latitude and longitude coordinates. Returns area weighted
    coefficient.
    weighted applies only to the centering (weighted mean subtraction) and not
    to the correlation itself.
    """
    weights = None
    if weighted:
        weights = iris.analysis.cartography.area_weights(cube1)
    if centered:
        dims = ["latitude", "longitude"]
        cube1 -= cube1.collapsed(dims, MEAN, weights=weights)
        cube2 -= cube2.collapsed(dims, MEAN, weights=weights)
    result = pearsonr(cube1, cube2, corr_coords=["latitude", "longitude"])
    return result.data


def ds_markers(extra_labels):
    datasets = set()
    for labels in extra_labels.values():
        datasets.update(labels)
    markers = {
        "ERA5": "o",
        "ERA5-Land": "s",
        "CRU": "x",
    }
    more_markers = iter("v^<>1sp*hH+xDd|_")
    for dataset in datasets:
        if dataset not in markers:
            markers[dataset] = next(more_markers)
    return markers


def process(metas, cfg, key=None):
    results = defaultdict(dict)
    # reference = shared.select_metadata(metas, dataset=cfg["reference"])
    cfg["variables"] = shared.group_metadata(metas, "short_name").keys()
    extra_labels = defaultdict(list)
    # print(
    #     [
    #         m["short_name"]
    #         for m in shared.group_metadata(metas, "dataset")["ERA5"]
    #     ]
    # )
    for var, var_metas in shared.group_metadata(metas, "short_name").items():
        logger.info("Processing %s (%s datasets)", var, len(var_metas))
        reference = ut.select_single_metadata(
            var_metas, dataset=cfg["reference"], short_name=var
        )
        ref_cube = iris.load_cube(reference["filename"])
        results[var] = defaultdict(list)
        for ds_meta in var_metas:
            if ds_meta["dataset"] in ["MMM", cfg["reference"]]:
                continue
            ds_meta["project"] = ds_meta.get("project", "unknown")
            # TODO: this need to be fixed in diag_pet.R:
            if var == "evspsblpot" and ds_meta["project"] == "unknown":
                ds_meta["project"] = "CMIP6"
            cube = iris.load_cube(ds_meta["filename"])
            # print(ds_meta['filename'])
            # print(cube)
            ds_result = float(pattern_correlation(ref_cube, cube))
            if ds_meta["project"] in cfg.get("projects", ["CMIP6"]):
                results[var][ds_meta["project"]].append(ds_result)
            elif cfg.get("extra_datasets", True):
                results[var]["extra"].append(ds_result)
                extra_labels[var].append(ds_meta["dataset"])
            # print("RSULTS")
            # print(ds_result)
    title = None
    if cfg.get("plot_title", False):
        title = f"Pattern Correlation with {cfg['reference']} ({key})"
    plot(
        cfg,
        results,
        f"pattern_correlation_{key}",
        extra_labels=extra_labels,
        ylims=(0, 1),
        title=title,
    )


def process_relative_change(metas, cfg):
    results = defaultdict(dict)
    cfg["variables"] = shared.group_metadata(metas, "short_name").keys()
    extra_labels = defaultdict(list)
    for var, var_metas in shared.group_metadata(metas, "short_name").items():
        logger.info("Processing %s (%s datasets)", var, len(var_metas))
        results[var] = defaultdict(list)
        for ds_meta in var_metas:
            ds_meta["project"] = ds_meta.get("project", "unknown")
            # TODO: this need to be fixed in diag_pet.R:
            if var == "evspsblpot" and ds_meta["project"] == "unknown":
                ds_meta["project"] = "CMIP6"
            cube = iris.load_cube(ds_meta["filename"])
            iris.analysis.maths.abs(cube, in_place=True)
            rel = cube.collapsed(["latitude", "longitude"], iris.analysis.MEAN)
            print(rel)
            ds_result = float(rel.data)
            if ds_meta["project"] in cfg.get("projects", ["CMIP6"]):
                results[var][ds_meta["project"]].append(ds_result)
            elif cfg.get("extra_datasets", True):
                results[var]["extra"].append(ds_result)
                extra_labels[var].append(ds_meta["dataset"])
    title = None
    if cfg.get("plot_title", False):
        title = "Relative Change over 10 Years"
    plot(
        cfg,
        results,
        "relative_change",
        extra_labels=extra_labels,
        ylims=(0, 10),
        title=title,
    )


def plot(
    cfg,
    results,
    fname,
    extra_labels=None,
    title="Pattern Correlation",
    ylims=None,
):
    if "order" in cfg:
        if sorted(cfg["order"]) != sorted(results.keys()):
            logger.warning(
                "Order does not match result keys: %s, %s",
                cfg["order"],
                results.keys(),
            )
        else:
            results = {key: results[key] for key in cfg["order"]}

    fig, axes = plt.subplots(figsize=(4, 3))
    plot_kwargs = dict(marker="_", linestyle="", markersize="5")
    plot_kwargs.update(cfg.get("plot_kwargs", {}))
    var_count = len(results.keys())
    # TODO: Start with simple case for 2 groups (add multiple groups later)
    projects = cfg.get("projects", ["CMIP6"])
    # proj_count = len(groups)
    # shift = 0.8 / group_count

    axes.set_title(title)
    axes.plot([], **plot_kwargs)
    axes.set_xticks(range(var_count))
    labels = results.keys()
    if "labels" in cfg:
        labels = [cfg["labels"].get(lab, lab) for lab in labels]
    axes.set_xticklabels(labels, rotation=90, ha="right")
    axes.set_xlim(-0.5, var_count - 0.5)
    if ylims is not None:
        axes.set_ylim(*ylims)
    columns = len(projects)
    if cfg.get("extra_datasets", True):
        columns += 1
    # single iteration yet
    for j, project in enumerate(projects):
        # NOTE: dict needs to be complete (each pair of var/proj required)
        # TODO: make results a xarray dataset instead?
        y_pos = []
        x_pos = []
        shift = -0.2 + 0.4 * j / columns
        for i, var in enumerate(results.keys()):
            y_pos.extend(results[var][project])
            x_pos.extend([i + shift] * len(results[var][project]))
        label = project
        axes.plot(x_pos, y_pos, **plot_kwargs, label=label)
    # Add extra datasets
    if cfg.get("extra_datasets", True):
        markers = ds_markers(extra_labels)
        scatters = defaultdict(list)
        print("plotting extra data")
        for i, var in enumerate(results.keys()):
            # TODO add label here manually?
            y_pos = results[var]["extra"]
            print(var)
            print(y_pos)
            for j, value in enumerate(y_pos):
                scatters[extra_labels[var][j]].append((i + 0.2, value))
        for ds, values in scatters.items():
            label = ds
            if label.startswith("CDS"):
                label = "CDS-SM"
            x_pos, y_pos = zip(*values)
            axes.scatter(
                x_pos, y_pos, marker=markers[ds], color="gray", label=label
            )
    axes.legend()
    if cfg.get("plot_properties"):
        axes.set(**cfg["plot_properties"])
    # fig.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    # fig.tight_layout()
    axes.grid(axis="y")
    shared.save_figure(fname, {}, cfg, figure=fig, bbox_inches="tight")


def main(cfg):
    metas = cfg["input_data"].values()
    # print(shared.group_metadata(metas, "short_name")["evspsblpot"])
    if not cfg.get("group_by"):
        process(metas, cfg)
        if cfg.get("relative_change", False):
            process_relative_change(metas, cfg)
    else:
        grouped = shared.group_metadata(metas, cfg["group_by"])
        groups = cfg.get("groups", list(grouped.keys()))
        for key, group in grouped.items():
            if key is None or key not in groups:
                continue
            if cfg.get("groups", None) is not None:
                if key not in cfg["groups"]:
                    continue
            process(group, cfg, key=key)
            if key == "percent" and cfg.get("relative_change", False):
                process_relative_change(group, cfg)


if __name__ == "__main__":
    with shared.run_diagnostic() as cfg:
        main(cfg)
