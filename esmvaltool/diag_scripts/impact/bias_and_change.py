"""Calculate and plot bias and change for each model.

Configuration options in recipe
-------------------------------
alias_facets: dict, optional
    Mapping from facets to CSV column names. The values of these facets are
    joined with an underscore to build a unique identifier for each model run,
    which is written to the ``dataset`` column of the CSV output file. The
    facet values are also written to the columns given in the mapping, so
    ``dataset`` cannot be used as a column name.
    Datasets that lack any of these facets (e.g. observations) are identified
    by their ``alias``. By default,
    ``{project: project, dataset: model, ensemble: member}``.
notes: list[str], optional
    Notes describing how the results were computed, e.g. the reference
    dataset and periods. These are stored in the Vega-Lite specification
    and shown as bullet points by the interactive viewer.

Output
------
Besides the figures and tables, the script writes a CSV file with the bias
and change of each model run and, for each project in the input data, a
Vega-Lite specification ``vegalite_spec_<project>.json`` with the data
embedded. These specifications are used by the interactive viewer at
https://github.com/ESMValGroup/C4I-Integration. They also record the name of
the recipe output directory, so the viewer can link to the ``index.html`` of
the recipe run that produced the data.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(Path(__file__).stem)

VARIABLE_LABELS = {
    "tas": "Temperature (K)",
    "pr": "Precipitation (kg/m2/s)",
}

METRIC_LABELS = {
    "bias": "Bias (RMSD of all gridpoints)",
    "change": "Mean change (Future - Reference)",
}

VEGALITE_TITLES = {
    "tas": "Temperature (K)",
    "pr": "Precipitation (mm/day)",
}

type DiagnosticConfig = dict[str, Any]
"""Diagnostic script configuration."""

type FacetMapping = dict[str, dict[str, Any]]
"""Mapping from dataset alias to CSV column name to facet value."""

DEFAULT_ALIAS_FACETS = {
    "project": "project",
    "dataset": "model",
    "ensemble": "member",
}


def log_provenance(
    filename: str,
    ancestors: list[str],
    caption: str,
    cfg: DiagnosticConfig,
) -> None:
    """Create a provenance record for the output file."""
    provenance = {
        "caption": caption,
        "domains": ["reg"],
        "authors": ["kalverla_peter"],
        "projects": ["isenes3"],
        "ancestors": ancestors,
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance)


def load_data(
    metadata: list[dict[str, Any]],
    alias_facets: dict[str, str],
) -> tuple[xr.DataArray, list[str], FacetMapping]:
    """Load all files from metadata into an Xarray dataset.

    ``metadata`` is a list of dictionaries with dataset descriptors.
    ``alias_facets`` is a mapping from facets that define the alias of
    each dataset to the corresponding column names in the CSV output.

    Returns the data, the ancestor files, and a mapping from alias to a
    mapping from column name to facet value.
    """
    data_arrays = []
    identifiers = []
    ancestors = []
    facets = {}

    for infodict in metadata:
        if all(infodict.get(facet) is not None for facet in alias_facets):
            alias = "_".join(str(infodict[facet]) for facet in alias_facets)
            facets[alias] = {
                column: infodict[facet]
                for facet, column in alias_facets.items()
            }
        else:
            alias = infodict["alias"]
            facets[alias] = {}
        # The project is always needed to split the output per project
        facets[alias].setdefault("project", infodict["project"])
        input_file = infodict["filename"]
        short_name = infodict["short_name"]

        xrds = xr.open_dataset(input_file)
        xrda = xrds[short_name]

        # Make sure datasets can be combined
        redundant_dims = np.setdiff1d(xrda.coords, xrda.dims)
        xrda = xrda.drop_vars(redundant_dims)

        data_arrays.append(xrda)
        identifiers.append(alias)
        ancestors.append(input_file)

    # Combine along a new dimension
    data_array = xr.concat(data_arrays, dim="dataset")
    if len(set(identifiers)) != len(identifiers):
        duplicates = sorted(
            {i for i in identifiers if identifiers.count(i) > 1},
        )
        msg = (
            f"Datasets {duplicates} are not uniquely identified by facets "
            f"{list(alias_facets)}, please add more facets to the "
            "'alias_facets' option of the diagnostic script."
        )
        raise ValueError(msg)
    data_array["dataset"] = identifiers

    return data_array, ancestors, facets


def plot_scatter(
    tidy_df: pd.DataFrame,
    ancestors: list[str],
    cfg: DiagnosticConfig,
) -> None:
    """Plot bias on one axis and change on the other."""
    grid = sns.relplot(
        data=tidy_df,
        x=METRIC_LABELS["bias"],
        y=METRIC_LABELS["change"],
        hue="dataset",
        col="variable",
        facet_kws={"sharex": False, "sharey": False},
        kind="scatter",
    )

    filename = get_plot_filename("bias_vs_change", cfg)
    grid.fig.savefig(filename, bbox_inches="tight")

    caption = "Bias and change for each variable"
    log_provenance(filename, ancestors, caption, cfg)


def plot_table(
    dataframe: pd.DataFrame,
    ancestors: list[str],
    cfg: DiagnosticConfig,
) -> None:
    """Render pandas table as a matplotlib figure."""
    dataframe = dataframe.reset_index()
    cell_text = [
        [
            f"{value:.3g}" if isinstance(value, float) else str(value)
            for value in row
        ]
        for row in dataframe.itertuples(index=False)
    ]

    # Size the figure to the table, so the text does not need to be shrunk
    nrows = len(cell_text) + 1
    fig, axes = plt.subplots(figsize=(10, 0.3 * nrows))
    axes.set_axis_off()
    table = axes.table(
        cellText=cell_text,
        colLabels=dataframe.columns,
        loc="center",
    )
    table.auto_set_font_size(value=False)
    table.set_fontsize(10)
    table.auto_set_column_width(range(len(dataframe.columns)))
    table.scale(1, 1.5)

    filename = get_plot_filename("table", cfg)
    fig.savefig(filename, bbox_inches="tight")

    caption = "Bias and change for each variable"
    log_provenance(filename, ancestors, caption, cfg)


def plot_htmltable(
    dataframe: pd.DataFrame,
    ancestors: list[str],
    cfg: DiagnosticConfig,
) -> None:
    """Render pandas table as html output.

    # https://pandas.pydata.org/pandas-docs/stable/user_guide/style.html
    """
    styles = [
        {"selector": ".index_name", "props": [("text-align", "right")]},
        {"selector": ".row_heading", "props": [("text-align", "right")]},
        {"selector": "td", "props": [("padding", "3px 25px")]},
    ]

    styled_table = (
        dataframe.pivot_table(index="dataset", columns="variable")
        .reindex(columns=VARIABLE_LABELS.values(), level="variable")
        .style.set_table_styles(styles)
        .background_gradient(cmap="RdYlGn", low=0, high=1, axis=0)
        .format("{:.2e}", na_rep="-")
        .to_html()
    )

    filename = get_diagnostic_filename("bias_vs_change", cfg, extension="html")
    Path(filename).write_text(styled_table, encoding="utf-8")

    caption = "Bias and change for each variable"
    log_provenance(filename, ancestors, caption, cfg)


def make_tidy(dataset: xr.Dataset) -> pd.DataFrame:
    """Convert xarray data to tidy dataframe."""
    return (
        dataset.rename(VARIABLE_LABELS)
        .to_dataframe()
        .reset_index()
        .melt(id_vars=["metric", "dataset"], var_name="variable")
        .pivot_table(
            index=["dataset", "variable"],
            columns="metric",
            values="value",
        )
        .reindex(VARIABLE_LABELS.values(), level="variable")
    )


def make_wide(
    dataframe: pd.DataFrame,
    facets: FacetMapping,
) -> pd.DataFrame:
    """Convert tidy dataframe to one row per dataset with facets as columns."""
    # one column per variable and metric, e.g. tas_bias
    dataframe = dataframe.pivot_table(index="dataset", columns="variable")
    short_names = {label: name for name, label in VARIABLE_LABELS.items()}
    metric_names = {label: name for name, label in METRIC_LABELS.items()}
    dataframe.columns = [
        f"{short_names[variable]}_{metric_names[metric]}"
        for metric, variable in dataframe.columns
    ]
    dataframe = dataframe[["tas_bias", "pr_bias", "tas_change", "pr_change"]]

    # metadata in separate columns
    dataframe = dataframe.join(pd.DataFrame.from_dict(facets, orient="index"))

    # kg/m2/s to mm/day
    dataframe[["pr_bias", "pr_change"]] *= 24 * 60 * 60

    return dataframe


def save_csv(
    dataframe: pd.DataFrame,
    ancestors: list[str],
    cfg: DiagnosticConfig,
) -> None:
    """Save output for use in Climate4Impact preview page."""
    filename = get_diagnostic_filename("recipe_output", cfg, extension="csv")
    caption = "Bias and change for each variable"
    dataframe.to_csv(filename)
    log_provenance(filename, ancestors, caption, cfg)


def build_vegalite_spec(
    dataframe: pd.DataFrame,
    project: str,
    alias_facets: dict[str, str],
    notes: list[str],
    recipe_output: str,
) -> dict[str, Any]:
    """Build a Vega-Lite specification for the Climate4Impact viewer.

    ``dataframe`` is the output of :func:`make_wide` for a single project.
    The data is embedded in the specification, so it can be used without
    any other files. ``recipe_output`` is the name of the recipe output
    directory.
    """
    color_field = alias_facets.get("dataset", "dataset")
    facet_columns = list(dict.fromkeys([*alias_facets.values(), "project"]))

    def scatter(variable: str) -> dict[str, Any]:
        """Plot the bias versus the change for a single variable."""
        return {
            "title": VEGALITE_TITLES[variable],
            "mark": {"type": "circle", "size": 150},
            "encoding": {
                "x": {
                    "field": f"{variable}_bias",
                    "type": "quantitative",
                    "title": METRIC_LABELS["bias"],
                    "scale": {"zero": False},
                },
                "y": {
                    "field": f"{variable}_change",
                    "type": "quantitative",
                    "title": METRIC_LABELS["change"],
                    "scale": {"zero": False},
                },
                "fill": {
                    "condition": {
                        "param": "brush",
                        "field": color_field,
                        "type": "nominal",
                    },
                    "value": "lightgray",
                },
                "stroke": {
                    "condition": {
                        "param": "query",
                        "value": "black",
                        "empty": False,
                    },
                    "value": "transparent",
                },
                "tooltip": [
                    {"field": column, "type": "nominal"}
                    for column in ["dataset", *facet_columns]
                ]
                + [
                    {
                        "field": f"{variable}_{metric}",
                        "type": "quantitative",
                        "title": metric,
                        "format": ".3g",
                    }
                    for metric in ("bias", "change")
                ],
            },
        }

    return {
        "$schema": "https://vega.github.io/schema/vega-lite/v6.json",
        "description": f"{project} data processed with ESMValTool.",
        "usermeta": {
            "project": project,
            "notes": notes,
            "recipe_output": recipe_output,
        },
        "data": {
            "values": json.loads(
                dataframe.reset_index().to_json(orient="records"),
            ),
        },
        "config": {
            "view": {"continuousWidth": 400, "continuousHeight": 300},
            "legend": {"disable": True},
        },
        "params": [
            {
                # Hold alt to select a range
                "name": "brush",
                "select": {
                    "type": "interval",
                    "on": "[pointerdown[event.altKey], pointerup] > pointermove",
                    "translate": (
                        "[pointerdown[event.altKey], pointerup] > pointermove"
                    ),
                    "zoom": "wheel![event.altKey]",
                },
            },
            {
                # Hold ctrl to pan and zoom
                "name": "panzoom",
                "select": {
                    "type": "interval",
                    "on": "[pointerdown[event.ctrlKey], pointerup] > pointermove",
                    "translate": (
                        "[pointerdown[event.ctrlKey], pointerup] > pointermove!"
                    ),
                    "zoom": "wheel![event.ctrlKey]",
                },
                "bind": "scales",
            },
            {
                # Click to select, hold shift to select multiple
                "name": "query",
                "select": {
                    "type": "point",
                    "fields": ["dataset"],
                    "clear": "dblclick",
                },
            },
        ],
        "hconcat": [scatter("tas"), scatter("pr")],
    }


def save_vegalite_specs(
    dataframe: pd.DataFrame,
    alias_facets: dict[str, str],
    notes: list[str],
    ancestors: list[str],
    cfg: DiagnosticConfig,
) -> None:
    """Save a Vega-Lite specification for each project."""
    # work_dir is <recipe output directory>/work/<diagnostic>/<script>
    recipe_output = Path(cfg["work_dir"]).parents[2].name
    for project, project_df in dataframe.groupby("project"):
        spec = build_vegalite_spec(
            project_df,
            project,
            alias_facets,
            notes,
            recipe_output,
        )
        filename = get_diagnostic_filename(
            f"vegalite_spec_{project}",
            cfg,
            extension="json",
        )
        Path(filename).write_text(json.dumps(spec, indent=2), encoding="utf-8")
        caption = f"Interactive plot of bias and change for {project}"
        log_provenance(filename, ancestors, caption, cfg)


def main(cfg: DiagnosticConfig) -> None:
    """Calculate, visualize and save the bias and change for each model."""
    metadata = cfg["input_data"].values()
    grouped_metadata = group_metadata(metadata, "variable_group")

    alias_facets = cfg.get("alias_facets", DEFAULT_ALIAS_FACETS)
    if "dataset" in alias_facets.values():
        msg = (
            "The 'dataset' column of the CSV output file is reserved for the "
            "dataset alias, please map facets to another column name in the "
            f"'alias_facets' option of the diagnostic script: {alias_facets}"
        )
        raise ValueError(msg)

    biases = {}
    changes = {}
    ancestors = []
    facets = {}
    for group, metadata in grouped_metadata.items():
        model_data, model_ancestors, model_facets = load_data(
            metadata,
            alias_facets,
        )
        ancestors.extend(model_ancestors)
        facets.update(model_facets)

        variable = model_data.name

        if group.endswith("bias"):
            # The distance_metric preprocessor prefixes the variable name
            variable = variable.removeprefix("rmse_")
            biases[variable] = model_data.rename(variable)

        elif group.endswith("change"):
            changes[variable] = model_data

        else:
            logger.warning(
                "Got input for variable group %s"
                " but I don't know what to do with it.",
                group,
            )

    # Combine all variables
    bias = xr.Dataset(biases)
    change = xr.Dataset(changes)
    combined = xr.concat([bias, change], dim="metric")
    combined["metric"] = [METRIC_LABELS["bias"], METRIC_LABELS["change"]]

    tidy_df = make_tidy(combined)
    plot_scatter(tidy_df, ancestors, cfg)
    plot_table(tidy_df, ancestors, cfg)
    plot_htmltable(tidy_df, ancestors, cfg)
    wide_df = make_wide(tidy_df, facets)
    save_csv(wide_df, ancestors, cfg)
    save_vegalite_specs(
        wide_df,
        alias_facets,
        cfg.get("notes", []),
        ancestors,
        cfg,
    )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
