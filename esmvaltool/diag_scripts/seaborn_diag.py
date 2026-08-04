#!/usr/bin/env python
"""Create arbitrary Seaborn plots.

Description
-----------
This diagnostic provides a high-level interface to Seaborn. For this, the input
data is arranged into a single :class:`pandas.DataFrame`, which is then used as
input for the Seaborn function defined by the option `seaborn_func`.

Caveats
-------
All datasets of a given variable must have the same units (e.g., it is not
allowed to use datasets with units `K` and datasets with units `°C` for the
variable `tas`).

Author
------
Manuel Schlund (DLR, Germany)

Configuration options in recipe
-------------------------------
add_ancillary_variables: bool, optional (default: False)
    Add :meth:`~iris.cube.Cube.ancillary_variables` to the main data frame.
    Note that this will assume that ancillary variables are identical across
    cubes within a group (see option `groupby_facet`). This equality is not
    checked!
add_aux_coords: bool, optional (default: False)
    Add :attr:`~iris.cube.Cube.aux_coords` to the main data frame. Note that
    this will assume that auxiliary coordinates are identical across cubes
    within a group (see option `groupby_facet`). This equality is not checked!
add_cell_measures: bool, optional (default: False)
    Add :meth:`~iris.cube.Cube.cell_measures` to the main data frame. Note that
    this will assume that cell measures are identical across cubes within a
    group (see option `groupby_facet`). This equality is not checked!
data_frame_ops: dict, optional
    Perform additional operations on the main data frame. Allowed operations
    are :meth:`pandas.DataFrame.query` (dict key `query`) and
    :meth:`pandas.DataFrame.eval` (dict key `eval`). Operations are defined by
    strings (dict values). Examples: ``{'query': 'latitude > 80', 'eval':
    'longitude = longitude - 180.0'}``.
dropna_kwargs: dict, optional
    Optional keyword arguments for :meth:`pandas.DataFrame.dropna` to drop
    missing values in the input data. If not given, do not drop NaNs. Note:
    NaNs are dropped after potential `data_frame_ops`.
facets_as_columns: list of str, optional
    Facets that will be added as a columns to the main data frame. Values for
    these facets must be identical across all datasets within a group (see
    option `groupby_facet`).
groupby_facet: str, optional (default: 'alias')
    Facet which is used to group input datasets when creating the main data
    frame. All datasets within a group are expected to have the same index
    after calling :func:`iris.pandas.as_data_frame` on them. These datasets
    within a group will then get merged (combined along axis 1, i.e., columns)
    into a single data frame per group. Finally, the data frames for all groups
    are concatenated (combined along axis 0, i.e., rows) into one main data
    frame. `groupby_facet` is also added as a column to this main data frame.
legend_title: str, optional (default: None)
    Title for legend. If ``None``, Seaborn will determine the legend title (if
    possible).
plot_filename: str, optional
    Filename for the final plot. By default, uses 'seaborn_(`seaborn_func`)'.
plot_object_methods: dict, optional
    Execute methods of the object returned by the plotting function
    (`seaborn_func`). This object will either be a
    :class:`matplotlib.axes.Axes` (e.g., :func:`~seaborn.scatterplot`,
    :func:`~seaborn.lineplot`), a :class:`seaborn.FacetGrid` (e.g.,
    :func:`~seaborn.relplot`, :func:`~seaborn.displot`), a
    :class:`seaborn.JointGrid` (e.g., :func:`~seaborn.jointplot`), or a
    :class:`seaborn.PairGrid` (e.g., :func:`~seaborn.pairplot`). Dictionary
    keys are method names, dictionary values function arguments (use a
    :obj:`dict` to specify keyword arguments). Example (for
    :func:`~seaborn.relplot`): ``{'set': {'xlabel': 'X [km]'}, 'set_titles':
    'Model {col_name}'}``.
reset_index: bool, optional (default: False)
    Put coordinate information of datasets into columns instead of (multi-)
    indices. This avoids the deletion of coordinate information if different
    groups of datasets have different dimensions but increases the memory
    footprint of this diagnostic.
write_netcdf: bool, optional (default: False)
    Output netCDF file for plotted data. What is written into the file
    is decided based on the pandas data frame and the
    seaborn kwargs: "x", "y", "hue" and "col".
    Because there is no direct way to write panda data frames into netCDF
    and to make data CF complient, the data frame is converted to an iris
    CubeList first. This is not possible for all data frames and often
    reset_index: true is required. Therefore the default is set to False.
savefig_kwargs: dict, optional
    Optional keyword arguments for :func:`matplotlib.pyplot.savefig`. By
    default, uses ``bbox_inches: tight, dpi: 300, orientation: landscape``.
seaborn_func: str
    Function used to plot the data. Must be a function of Seaborn. An overview
    of Seaborn's plotting functions is given `here
    <https://seaborn.pydata.org/tutorial/function_overview.html>`__.
seaborn_kwargs: dict, optional
    Optional keyword arguments for the plotting function given by
    `seaborn_func`. Must not include an argument called `data`. Example:
    ``{'x': 'variable_1', 'y': 'variable_2', 'hue': 'coord_1'}``. Note:
    variables (here: `variable_1` and `variable_2` are identified by their
    `variable_group` in the recipe, i.e., the keys that specify variable groups
    in `variables`.
seaborn_settings: dict, optional
    Options for :func:`seaborn.set_theme` (affects all plots).
suptitle: str or None, optional (default: None)
    Suptitle for the plot (see :func:`matplotlib.pyplot.suptitle`). If
    ``None``, do not create a suptitle. If the plot shows only a single panel,
    use `plot_object_methods` with ``{'set': {'title': 'TITLE'}}`` instead.
"""

from __future__ import annotations

import logging
from copy import deepcopy
from pathlib import Path
from pprint import pformat

import iris
import iris.pandas
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from iris.cube import CubeList
from matplotlib.colors import LogNorm, Normalize

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    io,
    run_diagnostic,
)

logger = logging.getLogger(Path(__file__).stem)

# Use the new behavior of :func:`iris.pandas.as_data_frame`
iris.FUTURE.pandas_ndim = True

# Save units of different variables
# Note: units must be unique across datasets of the same variable
UNITS: dict[str, str] = {}


def _create_plot(
    plot_func: callable,
    data_frame: pd.DataFrame,
    cfg: dict,
) -> None:
    """Create plot."""
    logger.debug(
        "Using main data frame as input for plotting:\n%s",
        data_frame,
    )

    # Plot
    plot_kwargs = cfg["seaborn_kwargs"]
    plot_func_str = cfg["seaborn_func"]
    if "data" in plot_kwargs:
        raise ValueError("'data' is an invalid argument for 'seaborn_kwargs'")
    logger.info(
        "Creating plot with\nseaborn.%s(\n    data=main_data_frame,\n%s\n)",
        plot_func_str,
        _get_str_from_kwargs(plot_kwargs),
    )

    plot_obj = plot_func(data=data_frame, **plot_kwargs)

    # Adjust plot appearance
    if cfg["plot_object_methods"]:
        for func_name, func_args in cfg["plot_object_methods"].items():
            if isinstance(func_args, dict):
                logger.debug(
                    "Running\n%s.%s(\n%s\n)",
                    type(plot_obj).__name__,
                    func_name,
                    _get_str_from_kwargs(func_args),
                )
                getattr(plot_obj, func_name)(**func_args)
            else:
                logger.debug(
                    "Running %s.%s(%r)",
                    type(plot_obj).__name__,
                    func_name,
                    func_args,
                )
                getattr(plot_obj, func_name)(func_args)
    if cfg["suptitle"] is not None:
        logger.debug("Setting `suptitle='%s'`", cfg["suptitle"])
        plt.suptitle(cfg["suptitle"], y=1.05)
    if cfg["legend_title"] is not None:
        _set_legend_title(plot_obj, cfg["legend_title"])
    if plot_func_str == "jointplot" and plot_kwargs["cbar"]:
        # Reposition colorbar so it is to the right of marginals plot
        plt.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)
        # get the current positions of the joint ax and the ax for
        # the marginal x
        pos_joint_ax = plot_obj.ax_joint.get_position()
        pos_marg_x_ax = plot_obj.ax_marg_x.get_position()
        # reposition the joint ax so it has the same width as the
        # marginal x ax
        plot_obj.ax_joint.set_position(
            [
                pos_joint_ax.x0,
                pos_joint_ax.y0,
                pos_marg_x_ax.width,
                pos_joint_ax.height,
            ],
        )
        # reposition the colorbar using new x positions and y
        # positions of the joint ax
        plot_obj.fig.axes[-1].set_position(
            [0.83, pos_joint_ax.y0, 0.07, pos_joint_ax.height],
        )

    # Save plot data
    if cfg["write_netcdf"]:
        _save_nc_data(data_frame, cfg)

    # Save plot
    plot_path = get_plot_filename(cfg["plot_filename"], cfg)
    plt.savefig(plot_path, **cfg["savefig_kwargs"])
    logger.info("Wrote %s", plot_path)
    plt.close()

    # Provenance tracking
    caption = f"Seaborn {cfg['seaborn_func']} for one or more dataset(s)"
    ancestors = [d["filename"] for d in cfg["input_data"].values()]
    provenance_record = {
        "ancestors": ancestors,
        "authors": ["schlund_manuel"],
        "caption": caption,
    }
    if cfg["write_netcdf"]:
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)
            provenance_logger.log(
                get_diagnostic_filename(cfg["plot_filename"], cfg),
                provenance_record,
            )
    else:
        with ProvenanceLogger(cfg) as provenance_logger:
            provenance_logger.log(plot_path, provenance_record)


def _get_grouped_data(cfg: dict) -> dict:
    """Get grouped input data."""
    groupby_facet = cfg["groupby_facet"]
    input_data = list(cfg["input_data"].values())

    # Check if necessary facets are present
    for dataset in input_data:
        if groupby_facet not in dataset:
            raise ValueError(
                f"groupby_facet '{groupby_facet}' is not available for "
                f"dataset {dataset['filename']}",
            )
        for facet in cfg["facets_as_columns"]:
            if facet not in dataset:
                raise ValueError(
                    f"Facet '{facet}' used for option 'facets_as_columns' is "
                    f"not available for dataset {dataset['filename']}",
                )

    # Group data accordingly
    grouped_data = group_metadata(
        input_data,
        groupby_facet,
        sort="filename",
    )

    return grouped_data


def _get_dataframe(cfg: dict) -> pd.DataFrame:
    """Get main :class:`pandas.DataFrame` used as input for plotting.

    Note
    ----
    Data is stored in long form, see also :func:`iris.pandas.as_data_frame`.
    """
    logger.info(
        "Grouping datasets by '%s' to create main data frame (data frames "
        "are merged [combined along axis 1, i.e., columns] within groups, "
        "then concatenated [combined along axis 0, i.e., rows] across groups)",
        cfg["groupby_facet"],
    )
    if cfg["add_aux_coords"]:
        logger.info("Adding aux_coords as columns")
    if cfg["add_cell_measures"]:
        logger.info("Adding cell_measures as columns")
    if cfg["add_ancillary_variables"]:
        logger.info("Adding ancillary_variables as columns")
    if cfg["facets_as_columns"]:
        logger.info("Adding facets as columns: %s", cfg["facets_as_columns"])

    grouped_data = _get_grouped_data(cfg)

    # Merge data frames within groups
    df_dict = {}
    for group, datasets in grouped_data.items():
        logger.info("Processing group '%s'", group)
        df_group = _get_df_for_group(cfg, group, datasets)
        df_dict[group] = df_group

    # Concatenate data frames across groups and use dtype 'category' for facet
    # columns to reduce memory usage and decrease computation times
    groupby_facet = cfg["groupby_facet"]
    df_main = pd.concat(df_dict.values(), ignore_index=cfg["reset_index"])
    df_main = df_main.astype(
        dict.fromkeys(cfg["facets_as_columns"] + [groupby_facet], "category"),
    )

    logger.info("Successfully retrieved main data frame from input data")
    logger.debug("Got main data frame:\n%s", df_main)
    return df_main


def _get_df_for_group(
    cfg: dict,
    group: str,
    datasets: list[dict],
) -> pd.DataFrame:
    """Extract :class:`pandas.DataFrame` for a single group of datasets.

    This merges (i.e., combines along axis 1 = columns) all data frames
    of individual datasets of a group.
    """
    df_group = pd.DataFrame()
    facets_as_columns: dict[str, str] = {}
    for dataset in datasets:
        filename = dataset["filename"]
        logger.info("Reading %s", filename)
        cube = iris.load_cube(filename)

        # Update units
        variable_group = dataset["variable_group"]
        units = dataset["units"]
        if variable_group in UNITS and UNITS[variable_group] != units:
            raise ValueError(
                f"Got duplicate units for variable '{variable_group}': "
                f"'{units}' and '{UNITS[variable_group]}'",
            )
        UNITS.setdefault(variable_group, units)

        # Get data frame for individual dataset with proper name
        df_dataset = iris.pandas.as_data_frame(
            cube,
            add_aux_coords=cfg["add_aux_coords"],
            add_cell_measures=cfg["add_cell_measures"],
            add_ancillary_variables=cfg["add_ancillary_variables"],
        )
        df_dataset = df_dataset.rename(
            {cube.name(): variable_group},
            axis="columns",
        )

        # Merge
        if df_group.empty:
            df_group = df_dataset
            facets_as_columns = {
                f: dataset[f] for f in cfg["facets_as_columns"]
            }
        else:
            # Make sure that dimensional coordinates match across cubes within
            # a group
            if not df_group.index.equals(df_dataset.index):
                raise ValueError(
                    f"Dimensions of cube {filename} differ from other cubes "
                    f"of group '{group}'. Cubes of that group:\n"
                    f"{pformat([d['filename'] for d in datasets])}",
                )

            # Make sure that facet values used as columns match across datasets
            # within a cube
            for facet, val in facets_as_columns.items():
                if dataset[facet] != val:
                    raise ValueError(
                        f"Facet value for facet '{facet}' (used by option "
                        f"'facets_as_columns') of dataset {filename} differs "
                        f"from value of other datasets of group '{group}': "
                        f"expected '{val}', got '{dataset[facet]}'. Datasets "
                        f"of that group:\n"
                        f"{pformat([d['filename'] for d in datasets])}",
                    )
            df_group = pd.merge(
                df_group,
                df_dataset,
                left_index=True,
                right_index=True,
                sort=False,
                suffixes=[None, "_DUPLICATE"],
            )

            # Assume that aux_coords, cell_measures, and ancillary_variables
            # (if requested) are equal across cubes within the group. Only add
            # them when they first appear.
            df_group = df_group.filter(regex="^(?!.*_DUPLICATE)")

    # Move dimensional coordinates from (multi-) index into columns if
    # requested
    if cfg["reset_index"]:
        df_group = df_group.reset_index()

    # Add additional information as column and save the data frame
    for facet, val in facets_as_columns.items():
        df_group[facet] = val
    if cfg["groupby_facet"] not in df_group.columns:
        df_group[cfg["groupby_facet"]] = group

    return df_group


def _get_default_cfg(cfg: dict) -> dict:
    """Get default options for configuration dictionary."""
    cfg = deepcopy(cfg)

    cfg.setdefault("add_ancillary_variables", False)
    cfg.setdefault("add_aux_coords", False)
    cfg.setdefault("add_cell_measures", False)
    cfg.setdefault("data_frame_ops", {})
    cfg.setdefault("dropna_kwargs", {})
    cfg.setdefault("facets_as_columns", [])
    cfg.setdefault("groupby_facet", "alias")
    cfg.setdefault("legend_title", None)
    cfg.setdefault("plot_object_methods", {})
    cfg.setdefault("plot_filename", f"seaborn_{cfg.get('seaborn_func', '')}")
    cfg.setdefault("reset_index", False)
    cfg.setdefault("write_netcdf", False)
    cfg.setdefault(
        "savefig_kwargs",
        {
            "bbox_inches": "tight",
            "dpi": 300,
            "orientation": "landscape",
        },
    )
    cfg.setdefault("seaborn_kwargs", {})
    cfg.setdefault("seaborn_settings", {})
    cfg.setdefault("suptitle", None)

    return cfg


def _get_str_from_kwargs(kwargs, separator="\n", prefix="    "):
    """Get overview string for kwargs."""
    return separator.join(f"{prefix}{k}={v!r}," for (k, v) in kwargs.items())


def _get_plot_func(cfg: dict) -> callable:
    """Get seaborn plot function."""
    if "seaborn_func" not in cfg:
        raise ValueError("Necessary option 'seaborn_func' missing")
    if not hasattr(sns, cfg["seaborn_func"]):
        raise AttributeError(
            f"Invalid seaborn_func '{cfg['seaborn_func']}' (must be a "
            f"function of the module seaborn; an overview of seaborn plotting "
            f"functions is given here: https://seaborn.pydata.org/tutorial/"
            f"function_overview.html)",
        )
    logger.info("Using plotting function seaborn.%s", cfg["seaborn_func"])
    return getattr(sns, cfg["seaborn_func"])


def _is_strictly_monotonic(arr):
    """Test if np.array is strictly monotonic."""
    result = np.all(np.diff(arr) > 0) | np.all(np.diff(arr) < 0)

    return result


def _modify_dataframe(data_frame: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Modify data frame according to the option ``data_frame_ops``."""
    allowed_funcs = ("query", "eval")

    # data_frame_ops
    for func, expr in cfg["data_frame_ops"].items():
        if func not in allowed_funcs:
            raise ValueError(
                f"Got invalid operation '{func}' for option 'data_frame_ops', "
                f"expected one of {allowed_funcs}",
            )
        op_str = f"'{func}' with argument '{expr}'"
        logger.info("Modifying main data frame through operation %s", op_str)
        data_frame = getattr(data_frame, func)(expr)
        logger.debug(
            "Main data frame after operation %s:\n%s",
            op_str,
            data_frame,
        )

    # dropna_kwargs
    if cfg["dropna_kwargs"]:
        logger.debug(
            "Running\ndata_frame.dropna(\n%s\n)",
            _get_str_from_kwargs(cfg["dropna_kwargs"]),
        )
        data_frame = data_frame.dropna(**cfg["dropna_kwargs"])
        logger.debug("Main data frame after dropna \n%s", data_frame)
    return data_frame


def _prepare_cube(cube, cubes_to_aux, cubes_to_coord, cfg):
    """Prepare cube to save data as netCDF."""
    cube.attributes.globals["seaborn_func"] = cfg["seaborn_func"]
    cube.attributes.globals["seaborn_kwargs"] = _get_str_from_kwargs(
        cfg["seaborn_kwargs"],
    )
    for auxcube in cubes_to_aux:
        aux_coord = iris.coords.AuxCoord(
            auxcube.data,
            long_name=auxcube.long_name,
        )
        # Add the auxiliary coordinate to the cube
        cube.add_aux_coord(aux_coord, data_dims=0)

    for dimcube in cubes_to_coord:
        dim_coord = iris.coords.DimCoord(
            dimcube.data,
            long_name=dimcube.long_name,
        )
        # Add the auxiliary coordinate to the cube
        cube.add_dim_coord(dim_coord, data_dims=0)

    return cube


def _save_nc_data(dframe: pd.DataFrame, cfg) -> None:
    """Save netCDF files for plot."""
    cubes_to_save = CubeList()
    cubes_to_aux = CubeList()
    cubes_to_coord = CubeList()

    strings_to_save = []
    for key in cfg["seaborn_kwargs"]:
        if key in ["x", "y", "hue", "col"]:
            strings_to_save.append(cfg["seaborn_kwargs"][key])

    for something in dframe:
        if something in strings_to_save:
            testcube = iris.pandas.as_cubes(dframe[something])[0]
            testcube.var_name = something
            testcube.remove_coord("unknown")

            if something in UNITS:
                testcube.units = UNITS[something]

            if something in ["shape_id", "dataset", "alias"]:
                testcube.data = testcube.data.astype(str)
                cubes_to_aux.append(testcube)
            elif something in [
                "latitude",
                "longitude",
                "height",
                "plev",
                "time",
            ] and _is_strictly_monotonic(testcube.data):
                cubes_to_coord.append(testcube)
            else:
                cubes_to_save.append(testcube)

    for cube in cubes_to_save:
        cube = _prepare_cube(
            cube,
            cubes_to_aux,
            cubes_to_coord,
            cfg,
        )

    io.iris_save(
        cubes_to_save,
        get_diagnostic_filename(cfg["plot_filename"], cfg),
    )


def _set_legend_title(plot_obj, legend_title: str) -> None:
    """Set legend title."""
    if hasattr(plot_obj, "get_legend"):  # Axes
        legend = plot_obj.get_legend()
    elif hasattr(plot_obj, "legend"):  # FacetGrid, PairGrid
        legend = plot_obj.legend
    elif isinstance(plot_obj, sns.axisgrid.JointGrid):  # JointGrid
        # Manually create a legend if needed in JointGrid
        handles, labels = plot_obj.ax_joint.get_legend_handles_labels()
        if handles and labels:
            legend = plot_obj.ax_joint.legend(
                handles=handles,
                labels=labels,
                title=legend_title,
            )
        else:
            legend = None
    else:
        raise ValueError(
            f"Cannot set legend title, `{type(plot_obj).__name__}` does not "
            f"support legends",
        )
    if legend is None:
        raise ValueError(
            "Cannot set legend title, plot does not contain legend",
        )
    logger.debug("Setting `legend_title='%s'`", legend_title)
    legend.set_title(legend_title)


def _validate_config(cfg: dict) -> dict:
    """Validate configuration dictionary."""
    cfg = deepcopy(cfg)

    # seaborn_kwargs: hue_norm
    if "hue_norm" in cfg["seaborn_kwargs"]:
        hue_norm = cfg["seaborn_kwargs"]["hue_norm"]
        if isinstance(hue_norm, str):
            vmin = cfg["seaborn_kwargs"].pop("vmin", None)
            vmax = cfg["seaborn_kwargs"].pop("vmax", None)
            if hue_norm == "linear":
                hue_norm = Normalize(vmin=vmin, vmax=vmax)
            elif hue_norm == "log":
                hue_norm = LogNorm(vmin=vmin, vmax=vmax)
            else:
                raise ValueError(
                    f"String value for `hue_norm` can only be `linear` or "
                    f"`log`, got `{hue_norm}`",
                )
            cfg["seaborn_kwargs"]["hue_norm"] = hue_norm
        if isinstance(hue_norm, list):
            cfg["seaborn_kwargs"]["hue_norm"] = tuple(hue_norm)

    return cfg


def main(cfg: dict) -> None:
    """Run diagnostic."""
    cfg = _get_default_cfg(cfg)
    cfg = _validate_config(cfg)

    sns.set_theme(**cfg["seaborn_settings"])
    plot_func = _get_plot_func(cfg)

    df_main = _get_dataframe(cfg)
    df_main = _modify_dataframe(df_main, cfg)

    _create_plot(plot_func, df_main, cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
