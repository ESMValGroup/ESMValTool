#!/usr/bin/env python
"""Plot carbon cycle--related multi-model mean time series.

Description
-----------
Plot multi-model mean time series in different panels based on different
(customizable) model properties (e.g., project, experiment, etc.). Moreover,
the lines can be adjusted based on another property. Additional facets can be
added to the datasets with a ``*.yml`` file (see e.g.,
``diag_scripts/carbon_cycle/carbon_model_characteristics.yml``). Reference
datasets can be specified with ``reference_datasets``.

Author
------
Manuel Schlund (DLR)

Projects
--------
Eval4CMIP
4C

Configuration options in recipe
-------------------------------
additional_facets_file: str, optional
    Additional facets that are added to the input datasets. See
    ``esmvaltool/diag_scripts/carbon_cycle/carbon_model_characteristics.yml``
    for an example. Accepts absolute paths or paths relative to
    ``esmvaltool/diag_scripts/carbon_cycle/``.
aliases: dict, optional
    Aliases used for labels, e.g., ``'project=CMIP6': CMIP6`` or
    ``'project=CMIP5_exp=historical': CMIP5 historical``. Can be elements
    distinguished by ``facets_to_separate_lines``,
    ``facets_to_separate_panels``, reference datasets, units, var_names,
    ``'All'``, etc.
facets_to_separate_lines: list of str, optional (default: False)
    Facets used to create the individual lines in the plots with different line
    properties. These properties can be customized with ``plot_kwargs``.  Each
    unique combination of elements will get individual line properties.  At
    most 3 different lines are supported.  If ``False``, do not use different
    line properties. All facets listed must be given for every input dataset
    (exception: datasets that are listed as ``reference_datasets``).
facets_to_separate_panels: list of str, optional (default: ['project'])
    Facets used to create the individual panels. Each unique combination of
    elements will get an individual panel. If ``False``, only a single panel is
    plotted. All facets listed must be given for every input dataset
    (exception: datasets that are listed as ``reference_datasets``).
gridspec_kwargs: dict, optional
    Keyword arguments for :func:`matplotlib.gridspec.GridSpec`.
input_units: str, optional
    Units in which the input data is given. This is useful if the metadata of
    the input data does not reflect the true units, e.g., when preprocessors
    like ``area_statistics`` with ``operator: sum`` are used.
legend_kwargs: dict, optional
    Optional keyword arguments of :func:`matplotlib.pyplot.legend`.
legend: bool or str, optional (default: 'all')
    Defines which elements are shows in the legend. If ``False``, do not show
    legends. If ``True`` or ``'all'``, show every single element that is
    plotted in the different panels and lines, i.e., every single combination
    of the unique elements from ``facets_to_separate_panels`` and
    ``facets_to_separate_lines``. If ``'panels'``, only show elements
    distinguished by ``facets_to_separate_panels``. If ``'lines'``, only show
    elements distinguished by ``facets_to_separate_lines``. If
    ``'panels+lines'``, show the elements from the options ``'panels'`` and
    ``'lines'``. In contrast to ``'All'``, this does not show all possible
    combinations.  Reference datasets are always shown.
n_columns: int, optional (default: 2)
    Number of columns used for the panels.
panel_text_position: list of float, optional (default: [0.02, 0.88])
    xy-position (in axes coordinates) for the panel heading.
plot_hline: float, optional
    If given, plot horizontal line at the given y-position.
plot_kwargs: dict of dict, optional
    Custom plot kwargs for :func:`matplotlib.pyplot.plot` for the different
    elements to plot in the different panels and lines. Must inlcude the unique
    elements identified by ``facets_to_separate_panels`` and
    ``facets_to_separate_lines``  as keys (individual facets separated by
    ``'_'``). Might also include plot kwargs for the reference datasets listed
    in ``reference_datasets`` and/or an entry ``default``.
plot_units: str, optional
    Units used to plot the data.
plot_vline: float, optional
    If given, plot vertical line at the given x-position.
pyplot_functions: dict, optional
    Optional calls of functions of :mod:`matplotlib.pyplot` (keys) with
    arbitrary arguments (values) used for plot appearance. Examples: ``title``,
    ``suptitle``, ``xlabel``, ``xlim``, etc.
reference_datasets: list of str, optional
    Reference datasets that are plotted in every panel.
savefig_kwargs: dict, optional
    Keyword arguments for :func:`matplotlib.pyplot.savefig`.
seaborn_settings: dict, optional
    Options for :func:`seaborn.set` (affects all plots).
separate_summary_plot: bool, optional (default: False)
    If ``True``, create extra summary plot (i.e., extra file) of all panels. If
    ``False``, use first panel for summary plot.
figure_kwargs: dict, optional
    Keyword arguments for :func:`matplotlib.pyplot.figure`.
uncertainty_kwargs: bool or dict, optional
    Controls plotting of uncertainties (standard deviations). If ``False``, do
    not plot uncertainties.  If ``True`` or :obj:`dict`, plot uncertainties.
    The :obj:`dict` defines optional keyword arguments for
    :func:`matplotlib.pyplot.fill_between`.
x_coord: str, optional (default: 'time')
    Name of the x-coordinate used for the plots.

"""

import logging
import warnings
from copy import deepcopy
from pathlib import Path

import iris
import iris.cube
import iris.plot
import iris.quickplot
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import yaml
from cf_units import Unit
from esmvalcore.preprocessor import multi_model_statistics
from matplotlib.lines import Line2D

# import iris.plot as iplt
from esmvaltool.diag_scripts.shared import (
    get_plot_filename,
    group_metadata,
    run_diagnostic,
    select_metadata,
)

logger = logging.getLogger(Path(__file__).stem)


def _add_facets(cfg, cfg_option, new_facet_name, input_data):
    """Add facets to input data."""
    # No separation necessary
    if cfg[cfg_option] is False:
        for dataset in input_data:
            dataset[new_facet_name] = None
        return input_data

    # Separation necessary
    for dataset in input_data:
        if dataset["ref"]:
            continue
        new_facets = []
        for facet in cfg[cfg_option]:
            if facet not in dataset:
                raise ValueError(
                    f"Facet '{facet}' in {cfg_option} not given in "
                    f"non-reference dataset {dataset['dataset']} ("
                    f"{dataset['project']})"
                )
            new_facets.append(f"{facet}={str(dataset[facet])}")
        dataset[new_facet_name] = "_".join(new_facets)

    return input_data


def _create_legend(cfg, mpl_obj, input_data):
    """Create legend with desired content."""
    # TODO: option for panels+lines
    if cfg["legend"] is False:
        return None
    ref_datasets = select_metadata(input_data, ref=True)
    non_ref_datasets = select_metadata(input_data, ref=False)
    n_lines = len(group_metadata(non_ref_datasets, "line"))

    handles = []
    labels = []

    # Reference datasets
    for dataset in ref_datasets:
        ref_name = dataset["dataset"]
        plot_kwargs = get_plot_kwargs(cfg, "default")
        plot_kwargs.update(get_plot_kwargs(cfg, ref_name))
        handles.append(Line2D([0.0], [0.0], **plot_kwargs))
        labels.append(alias(cfg, ref_name))

    # Panel information
    if cfg["legend"] in ["all", "panels"]:
        grouped_data = group_metadata(non_ref_datasets, "panel")
        for panel_name, datasets in grouped_data.items():
            if panel_name is None:
                panel_name = "All"
            plot_kwargs = get_plot_kwargs(cfg, "default")
            plot_kwargs.update(get_plot_kwargs(cfg, panel_name))
            handles.append(Line2D([0.0], [0.0], **plot_kwargs))
            labels.append(alias(cfg, panel_name))

            # All information
            if cfg["legend"] == "panels":
                continue
            if n_lines > 1:
                subgrouped_data = group_metadata(datasets, "line", sort=True)
                for line_name in subgrouped_data:
                    sub_label = panel_name + "_" + line_name
                    sub_plot_kwargs = dict(plot_kwargs)
                    sub_plot_kwargs.update(get_plot_kwargs(cfg, line_name))
                    sub_plot_kwargs.update(get_plot_kwargs(cfg, sub_label))
                    handles.append(Line2D([0.0], [0.0], **sub_plot_kwargs))
                    labels.append(alias(cfg, sub_label))

    # Only line information in legend
    elif cfg["legend"] == "lines" and n_lines > 1:
        plot_kwargs = get_plot_kwargs(cfg, "default")
        plot_kwargs.update(get_plot_kwargs(cfg, "All"))
        handles.append(Line2D([0.0], [0.0], **plot_kwargs))
        labels.append(alias(cfg, "All"))

        grouped_data = group_metadata(non_ref_datasets, "line")
        for line_name in grouped_data:
            plot_kwargs = get_plot_kwargs(cfg, "default")
            plot_kwargs.update(get_plot_kwargs(cfg, line_name))
            handles.append(Line2D([0.0], [0.0], **plot_kwargs))
            labels.append(alias(cfg, line_name))

    # Unsupported legend option
    else:
        raise ValueError(f"Option '{cfg['legend']}' for legend not supported")

    # Create legend
    legend = mpl_obj.legend(handles, labels, **cfg["legend_kwargs"])
    return legend


def _plot_mm_cube(cfg, axes, datasets, uncertainty=True, **plot_kwargs):
    """Plot single multi-model cube."""
    cubes = iris.cube.CubeList([d["cube"] for d in datasets])

    # Calculate necessary statistics
    stats = ["mean"]
    if uncertainty and cfg["uncertainty_kwargs"] is not False:
        stats.append("std_dev")
    mm_cubes = multi_model_statistics(cubes, "full", stats)

    # Plot multi-model mean
    mmm_cube = mm_cubes["mean"]
    mmm_line = iris.plot.plot(
        mmm_cube.coord(cfg["x_coord"]), mmm_cube, axes=axes, **plot_kwargs
    )[0]
    x_data = mmm_line.get_xdata()

    # Plot uncertainties if desired
    if "std_dev" not in mm_cubes:
        return
    std_cube = mm_cubes["std_dev"]
    uncertainty_kwargs = dict(plot_kwargs)
    uncertainty_kwargs["alpha"] = 0.25
    uncertainty_kwargs["zorder"] = 1.1
    if isinstance(cfg["uncertainty_kwargs"], dict):
        uncertainty_kwargs.update(cfg["uncertainty_kwargs"])
    axes.fill_between(
        x_data,
        mmm_cube.data - std_cube.data,
        mmm_cube.data + std_cube.data,
        **uncertainty_kwargs,
    )


def _plot_ref_datasets(cfg, axes, input_data):
    """Plot reference datasets."""
    ref_datasets = select_metadata(input_data, ref=True)
    logger.info(ref_datasets)
    for dataset in ref_datasets:
        ref_name = dataset["dataset"]
        ref_cube = dataset["cube"]
        plot_kwargs = get_plot_kwargs(cfg, "default")
        plot_kwargs.update(get_plot_kwargs(cfg, ref_name))
        plot_kwargs["zorder"] = 2.2
        logger.info(dataset["dataset"])
        logger.info(ref_cube.data)
        iris.plot.plot(
            ref_cube.coord(cfg["x_coord"]), ref_cube, axes=axes, **plot_kwargs
        )


def _plot_summary_panel(cfg, axes, input_data, n_lines):
    """Plot summary panel in given axes."""
    logger.info("Plotting summary panel")
    non_ref_datasets = select_metadata(input_data, ref=False)
    grouped_data = group_metadata(non_ref_datasets, "panel")

    _plot_ref_datasets(cfg, axes, input_data)

    # Plot all panels iteratively
    for panel_name, datasets in grouped_data.items():
        plot_kwargs = get_plot_kwargs(cfg, "default")
        plot_kwargs.update(get_plot_kwargs(cfg, panel_name))
        plot_kwargs["zorder"] = 2.1
        _plot_mm_cube(cfg, axes, datasets, uncertainty=False, **plot_kwargs)
        plot_kwargs.pop("zorder", None)

        # Plot individual lines distinguished by facets_to_separate_lines if
        # necessary
        if n_lines > 1:
            subgrouped_data = group_metadata(datasets, "line", sort=True)
            for line_name, sub_datasets in subgrouped_data.items():
                sub_label = panel_name + "_" + line_name
                sub_plot_kwargs = dict(plot_kwargs)
                sub_plot_kwargs.update(get_plot_kwargs(cfg, line_name))
                sub_plot_kwargs.update(get_plot_kwargs(cfg, sub_label))
                _plot_mm_cube(
                    cfg,
                    axes,
                    sub_datasets,
                    uncertainty=False,
                    **sub_plot_kwargs,
                )
    plot_h_and_v_lines(cfg, axes, color="gray", linestyle="--", linewidth=1)
    axes.set_xlabel(alias(cfg, cfg["x_coord"]))
    axes.set_ylabel(
        f"{alias(cfg, input_data[0]['short_name'])} "
        f"[{alias(cfg, cfg['plot_units'])}]"
    )
    process_pyplot_functions(cfg)


def _unify_time_coord(cube):
    """Unify time coordinate of cube."""
    if not cube.coords("time", dim_coords=True):
        return
    time_coord = cube.coord("time")
    dates_points = time_coord.units.num2date(time_coord.points)
    dates_bounds = time_coord.units.num2date(time_coord.bounds)
    new_units = Unit("days since 1850-01-01 00:00:00")
    new_time_coord = iris.coords.DimCoord(
        new_units.date2num(dates_points),
        bounds=new_units.date2num(dates_bounds),
        var_name="time",
        standard_name="time",
        long_name="time",
        units=new_units,
    )
    coord_dims = cube.coord_dims("time")
    cube.remove_coord("time")
    cube.add_dim_coord(new_time_coord, coord_dims)


def add_additional_facets(input_data, path_to_facets):
    """Add additional facets to the input data."""
    # Get absolute path
    path_to_facets = Path(path_to_facets)
    if not path_to_facets.is_absolute():
        path_to_facets = (Path(__file__).parent / path_to_facets).resolve()

    # Read facets file and add information to input datasets
    with open(path_to_facets) as infile:
        facets_dict = yaml.safe_load(infile)
    logger.info("Read additional facets from %s", path_to_facets)
    input_data = deepcopy(input_data)
    grouped_data = group_metadata(input_data, "dataset")
    for dataset_name, datasets in grouped_data.items():
        if dataset_name not in facets_dict:
            continue
        for dataset in datasets:
            dataset.update(facets_dict[dataset_name])

    return input_data


def mark_reference_datasets(cfg, input_data):
    """Mark reference datasets."""
    input_data = deepcopy(input_data)
    for dataset in input_data:
        if dataset["dataset"] in cfg["reference_datasets"]:
            dataset["ref"] = True
        else:
            dataset["ref"] = False
    return input_data


def add_facets_for_plotting(cfg, input_data):
    """Add additional facets to the input data."""
    input_data = deepcopy(input_data)
    input_data = _add_facets(
        cfg, "facets_to_separate_panels", "panel", input_data
    )
    input_data = _add_facets(
        cfg, "facets_to_separate_lines", "line", input_data
    )
    non_ref_datasets = select_metadata(input_data, ref=False)
    lines = list(group_metadata(non_ref_datasets, "line"))
    if len(lines) > 3:
        raise ValueError(
            f"At most three different lines distinguished by "
            f"facets_to_separate_lines are currently supported got "
            f"{len(lines):d}: {lines}"
        )
    return input_data


def get_default_cfg(cfg):
    """Get default values for configuration :obj:`dict`."""
    cfg = deepcopy(cfg)

    # Kwargs
    cfg.setdefault("aliases", {})
    cfg.setdefault("figure_kwargs", {})
    cfg.setdefault("gridspec_kwargs", {})
    cfg.setdefault("legend_kwargs", {})
    cfg.setdefault("plot_kwargs", {})
    cfg.setdefault("pyplot_functions", {})
    cfg.setdefault("savefig_kwargs", {})
    cfg.setdefault("uncertainty_kwargs", False)

    # Others
    cfg.setdefault("facets_to_separate_lines", False)
    cfg.setdefault("facets_to_separate_panels", ["project"])
    cfg.setdefault("legend", "all")
    cfg.setdefault("n_columns", 2)
    cfg.setdefault("panel_text_position", [0.02, 0.88])
    cfg.setdefault("reference_datasets", [])
    cfg.setdefault("separate_summary_plot", False)
    cfg.setdefault("x_coord", "time")

    # Verification of options
    if cfg["legend"] is True:
        cfg["legend"] = "all"
    elif cfg["legend"] is False:
        pass
    else:
        allowed_legends = ["all", "panels", "lines"]
        if cfg["legend"] not in allowed_legends:
            raise ValueError(
                f"Got invalid str option '{cfg['legend']}' for legend, "
                f"expected bool or one of {allowed_legends}"
            )

    return cfg


def alias(cfg, key):
    """Get alias for key."""
    return cfg["aliases"].get(key, key)


def load_cubes(cfg, datasets):
    """Load all cubes and add as 'cube' facet."""
    datasets = deepcopy(datasets)
    for dataset in datasets:
        cube = iris.load_cube(dataset["filename"])
        if dataset["dataset"] == "GCP":
            cube = cube.collapsed(
                ["longitude", "latitude"], iris.analysis.MEAN
            )
            cube.data = cube.data * 148300000000000.0  # multiply by area
            iris.coord_categorisation.add_year(cube, "time")

        for coord in cube.coords(dim_coords=False):
            if coord.name() == cfg["x_coord"]:
                continue
            cube.remove_coord(coord)
        _unify_time_coord(cube)
        if "input_units" in cfg:
            cube.units = cfg["input_units"]
        if "plot_units" in cfg:
            cube.convert_units(cfg["plot_units"])
        dataset["cube"] = cube
        if dataset["dataset"] == "GCP":
            logger.info(dataset["dataset"])
            logger.info(cube.data)
    return datasets


def process_pyplot_functions(cfg):
    """Process arbitrary pyplot functions."""
    for func, arg in cfg["pyplot_functions"].items():
        getattr(plt, func)(arg)


def get_plot_kwargs(cfg, key):
    """Get plot_kwargs for a given panel."""
    all_plot_kwargs = cfg["plot_kwargs"]
    plot_kwargs = {}
    if key in all_plot_kwargs:
        plot_kwargs.update(all_plot_kwargs[key])
    if "linestyle" in plot_kwargs:
        if isinstance(plot_kwargs["linestyle"], list):
            plot_kwargs["linestyle"] = tuple(
                [
                    plot_kwargs["linestyle"][0],
                    tuple(plot_kwargs["linestyle"][1]),
                ]
            )
    return plot_kwargs


def plot_h_and_v_lines(cfg, axes, **line_kwargs):
    """Plot horizontal and vertical lines if desired."""
    if "plot_hline" in cfg:
        axes.axhline(cfg["plot_hline"], zorder=1.9, **line_kwargs)
    if "plot_vline" in cfg:
        axes.axvline(cfg["plot_vline"], zorder=1.9, **line_kwargs)


def create_time_series_plot(cfg, input_data):
    """Create time series plot."""
    non_ref_datasets = select_metadata(input_data, ref=False)
    n_panels = len(group_metadata(non_ref_datasets, "panel"))
    if n_panels == 1:
        plot_single_panel(cfg, input_data)
    else:
        plot_multiple_panels(cfg, input_data, n_panels)


def plot_single_panel(cfg, input_data):
    """Plot single-panel time series plot."""
    logger.info("Plotting single panel")
    (_, axes) = plt.subplots(**cfg["figure_kwargs"])

    # Plot reference datasets
    logger.info("Plotting reference datasets")
    _plot_ref_datasets(cfg, axes, input_data)

    # Plot line for all non-reference datasets
    non_ref_datasets = select_metadata(input_data, ref=False)
    grouped_data = group_metadata(non_ref_datasets, "line", sort=True)
    n_lines = len(grouped_data)
    logger.info("Plotting line containing all non-reference datasets")
    plot_kwargs = get_plot_kwargs(cfg, "default")
    plot_kwargs.update(get_plot_kwargs(cfg, "All"))
    plot_kwargs["zorder"] = 2.1
    _plot_mm_cube(cfg, axes, non_ref_datasets, **plot_kwargs)
    plot_kwargs.pop("zorder", None)

    # Plot lines for datasets distinguished by facets_to_separate_lines if
    # necessary
    if n_lines > 1:
        for line_name, datasets in grouped_data.items():
            logger.info("Plotting line %s", line_name)
            plot_kwargs = get_plot_kwargs(cfg, "default")
            plot_kwargs.update(get_plot_kwargs(cfg, line_name))
            _plot_mm_cube(cfg, axes, datasets, **plot_kwargs)

    # Plot appearance
    plot_h_and_v_lines(cfg, axes, color="gray", linestyle="--", linewidth=1)
    axes.set_xlabel(alias(cfg, cfg["x_coord"]))
    axes.set_ylabel(
        f"{alias(cfg, input_data[0]['short_name'])} "
        f"[{alias(cfg, cfg['plot_units'])}]"
    )
    process_pyplot_functions(cfg)
    _create_legend(cfg, axes, input_data)

    # Save plot
    plot_path = get_plot_filename("mmm_time_series", cfg)
    plt.savefig(plot_path, **cfg["savefig_kwargs"])
    logger.info("Wrote %s", plot_path)
    plt.close()


def plot_multiple_panels(cfg, input_data, n_panels):
    """Plot multi-panel time series plot."""
    ref_datasets = select_metadata(input_data, ref=True)
    non_ref_datasets = select_metadata(input_data, ref=False)
    grouped_data = group_metadata(non_ref_datasets, "panel")
    n_lines = len(group_metadata(non_ref_datasets, "line"))
    n_cols = cfg["n_columns"]

    # Plot summary panel in separate plot
    if cfg["separate_summary_plot"]:
        logger.info("Plotting %d panels and a separate summary plot", n_panels)
        n_rows = int(np.ceil(n_panels / n_cols))
        row_offset = 0

        (_, ax_sum) = plt.subplots(**cfg["figure_kwargs"])
        _plot_summary_panel(cfg, ax_sum, input_data, n_lines)
        _create_legend(cfg, ax_sum, input_data)
        sum_plot_path = get_plot_filename("summary_time_series", cfg)
        plt.savefig(sum_plot_path, **cfg["savefig_kwargs"])
        logger.info("Wrote %s", sum_plot_path)

        # Create new figure for multi-panel plot
        fig = plt.figure(**cfg["figure_kwargs"])
        gridspec = fig.add_gridspec(
            nrows=n_rows, ncols=n_cols, **cfg["gridspec_kwargs"]
        )

    # Plot summary panel in top panel (total number of panels = number of
    # subpanels [=n_panels] + 1
    else:
        logger.info(
            "Plotting %d panels (including summary panel at the top)",
            n_panels + 1,
        )
        n_rows = int(np.ceil(n_panels / n_cols)) + 1
        row_offset = 1

        fig = plt.figure(**cfg["figure_kwargs"])
        gridspec = fig.add_gridspec(
            nrows=n_rows, ncols=n_cols, **cfg["gridspec_kwargs"]
        )
        ax_sum = fig.add_subplot(gridspec[0, :])
        _plot_summary_panel(cfg, ax_sum, input_data, n_lines)
        ax_sum.set_xlabel("")  # This is set later (only for the bottom axes)

    # In order to use sharex and sharey with a gridspec, a reference axes needs
    # to be specified. To define this, use the first iteration in the following
    # loop
    ref_axes = None

    # Plot other panels
    for panel_idx, (panel_name, datasets) in enumerate(grouped_data.items()):
        logger.info("Plotting panel %s", panel_name)
        row_idx = panel_idx // n_cols + row_offset
        col_idx = panel_idx % n_cols
        if ref_axes is None:
            axes = fig.add_subplot(gridspec[row_idx, col_idx])
            ref_axes = axes
        else:
            axes = fig.add_subplot(
                gridspec[row_idx, col_idx], sharey=ref_axes, sharex=ref_axes
            )

        # Plot reference datasets
        logger.info("Plotting reference datasets")
        _plot_ref_datasets(cfg, axes, input_data)

        # Plot all datasets of specific panel
        plot_kwargs = get_plot_kwargs(cfg, "default")
        plot_kwargs.update(get_plot_kwargs(cfg, panel_name))
        plot_kwargs["zorder"] = 2.1
        _plot_mm_cube(cfg, axes, datasets, **plot_kwargs)
        plot_kwargs.pop("zorder", None)

        # Plot individual lines distinguished by facets_to_separate_lines if
        # necessary
        if n_lines > 1:
            subgrouped_data = group_metadata(datasets, "line", sort=True)
            for line_name, sub_datasets in subgrouped_data.items():
                sub_label = panel_name + "_" + line_name
                sub_plot_kwargs = dict(plot_kwargs)
                sub_plot_kwargs.update(get_plot_kwargs(cfg, line_name))
                sub_plot_kwargs.update(get_plot_kwargs(cfg, sub_label))
                _plot_mm_cube(cfg, axes, sub_datasets, **sub_plot_kwargs)

        # Plot appearance
        # axes.set_xlim(cfg['x_lim'][0], cfg['x_lim'][1]) #nopeeeeeeeeee
        axes.text(
            *cfg["panel_text_position"],
            alias(cfg, panel_name),
            transform=axes.transAxes,
        )
        plot_h_and_v_lines(
            cfg, axes, color="gray", linestyle="--", linewidth=1
        )
        axes.set_ylabel(
            f"{alias(cfg, input_data[0]['short_name'])} "
            f"[{alias(cfg, cfg['plot_units'])}]"
        )
        axes.set_xlabel(alias(cfg, cfg["x_coord"]))
        process_pyplot_functions(cfg)

        # Remove axis labels and tick labels for panels not located at the left
        # or bottom edge
        if row_idx != n_rows - 1:
            axes.set_xlabel("")
            # axes.set_xticklabels([]) does not work for time axis
            plt.setp(axes.get_xticklabels(), visible=False)
        if col_idx != 0:
            axes.set_ylabel("")
            plt.setp(axes.get_yticklabels(), visible=False)

    # Add legend to figure based on summary plot
    _create_legend(cfg, fig, input_data)

    # Save plot
    plot_path = get_plot_filename("mmm_time_series", cfg)
    plt.savefig(plot_path, **cfg["savefig_kwargs"])
    logger.info("Wrote %s", plot_path)
    plt.close()


def main(cfg):
    """Run the diagnostic."""
    warnings.filterwarnings(
        "ignore",
        message="Collapsing a non-contiguous coordinate",
        category=UserWarning,
        module="iris",
    )
    sns.set(**cfg.get("seaborn_settings", {}))
    cfg = get_default_cfg(cfg)

    # Read data, add additional facets, and load cubes
    input_data = list(cfg["input_data"].values())
    if "additional_facets_file" in cfg:
        input_data = add_additional_facets(
            input_data, cfg["additional_facets_file"]
        )
    input_data = mark_reference_datasets(cfg, input_data)
    input_data = add_facets_for_plotting(cfg, input_data)
    input_data = load_cubes(cfg, input_data)

    # Create time series plot
    create_time_series_plot(cfg, input_data)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
