"""Cross-region Quantile-MAE summary.

Reads ``quantile_mae_table.csv`` from each of its ancestor diagnostics
(one per region), tags rows with a region label, concatenates them into
a single long-form CSV ``qmae_cross_region_table.csv`` and renders a
companion heatmap showing the (method1 / method2) ratio as a
(region x derived-variable) matrix.

This merged table is produced inside the ESMValTool flow so the result is
provenance-tracked like the rest of the recipe.
"""

import logging
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared._base import ProvenanceLogger

logger = logging.getLogger(os.path.basename(__file__))

QMAE_VARIABLE_ORDER = ["wbgt", "pr", "sfcWind"]

QMAE_VARIABLE_LABELS = {
    "wbgt": "Summer WBGT",
    "pr": "Winter Precip.",
    "sfcWind": "Fall Wind Speed",
}


def _identify_region(ancestor_dir, region_labels):
    """Return the region label whose key appears in *ancestor_dir*.

    region_labels is a dict like
        {"quantileMAE_6M_CE":  "Central Europe",
         "quantileMAE_6M_IBE": "Iberic Peninsula",
         "quantileMAE_6M_SCA": "Scandinavia"}
    """
    path_str = str(ancestor_dir)
    for diag_key, label in region_labels.items():
        if diag_key in path_str:
            return label
    return None


def gather_per_region_qmae(cfg):
    """Walk each ancestor dir, read its ``quantile_mae_table.csv`` and tag
    the rows with the matching region label.
    """
    region_labels = cfg.get("region_labels", {})
    if not region_labels:
        raise ValueError(
            "qmae_cross_region requires a `region_labels` dict mapping "
            "ancestor diagnostic keys (e.g. 'quantileMAE_6M_CE') to "
            "human-readable region names."
        )

    ancestor_dirs = cfg.get("input_files", [])
    if not ancestor_dirs:
        raise RuntimeError("No ancestor input_files in cfg.")

    frames = []
    seen_csvs = set()

    for ancestor_dir in ancestor_dirs:
        region = _identify_region(ancestor_dir, region_labels)
        if region is None:
            logger.warning(
                "Could not match ancestor dir %s to any region label in %s",
                ancestor_dir, list(region_labels.keys()),
            )
            continue

        for csv_path in Path(ancestor_dir).rglob("quantile_mae_table.csv"):
            csv_str = str(csv_path.resolve())
            if csv_str in seen_csvs:
                continue
            seen_csvs.add(csv_str)
            df = pd.read_csv(csv_path)
            df.insert(0, "region", region)
            frames.append(df)
            logger.info("Loaded %s -> region=%s, %d rows",
                        csv_path, region, len(df))

    if not frames:
        raise RuntimeError(
            "No per-region quantile_mae_table.csv files found under any "
            "ancestor input_files."
        )

    combined = pd.concat(frames, ignore_index=True)
    return combined


def _coerce_numeric(df, methods):
    for m in methods:
        if m in df.columns:
            df[m] = pd.to_numeric(df[m], errors="coerce")
    return df


def save_merged_table(df, cfg):
    """Write the long-form merged table and a wide (region x var x method) view."""
    out_csv = Path(cfg["work_dir"]) / "qmae_cross_region_table.csv"
    df.to_csv(out_csv, index=False)
    logger.info("Saved merged QMAE table: %s", out_csv)

    # Convenience wide view: rows = (region, variable), columns = methods
    methods = cfg.get("methods", [])
    wide_cols = ["region", "variable"] + [m for m in methods if m in df.columns]
    if "metric" in df.columns:
        wide_cols.insert(2, "metric")
    wide = df[wide_cols].copy()
    wide_csv = Path(cfg["work_dir"]) / "qmae_cross_region_table_wide.csv"
    wide.to_csv(wide_csv, index=False)
    logger.info("Saved wide-format merged QMAE table: %s", wide_csv)
    return out_csv


def compute_ratio_matrix(df, method1, method2, regions, variables):
    """Return ratio matrix R[region, variable] = method1 / method2 .

    Lower ratio = method1 better.  NaN where either side missing.
    """
    epsilon = 1e-10
    mat = np.full((len(regions), len(variables)), np.nan)
    for i, region in enumerate(regions):
        for j, var in enumerate(variables):
            row = df[(df["region"] == region) & (df["variable"] == var)]
            if row.empty:
                continue
            row = row.iloc[0]
            v1 = row.get(method1, np.nan)
            v2 = row.get(method2, np.nan)
            try:
                v1 = float(v1)
                v2 = float(v2)
            except (ValueError, TypeError):
                continue
            if not np.isfinite(v1) or not np.isfinite(v2) or abs(v2) < epsilon:
                continue
            mat[i, j] = v1 / v2
    return mat


def _diverging_norm(vmin, vmax):
    return mcolors.TwoSlopeNorm(vmin=vmin, vcenter=1.0, vmax=vmax)


def plot_heatmap(ratio_mat, regions, variables, method1, method2, cfg,
                 abs_df=None):
    """Render the cross-region ratio heatmap.

    Annotates each cell with the ratio AND, in parentheses, the absolute
    QMAE values for method1 / method2 (helps the reader judge whether a
    ratio near 1 reflects two strong models or two weak ones).
    """
    vmin = float(cfg.get("heatmap_vmin", 0.80))
    vmax = float(cfg.get("heatmap_vmax", 1.20))
    cmap = plt.get_cmap("RdYlGn_r")
    norm = _diverging_norm(vmin, vmax)

    fig, ax = plt.subplots(
        figsize=(1.6 * len(variables) + 2.5, 0.9 * len(regions) + 2.0)
    )

    var_labels = [QMAE_VARIABLE_LABELS.get(v, v) for v in variables]

    im = ax.imshow(ratio_mat, cmap=cmap, norm=norm, aspect="auto")

    ax.set_xticks(np.arange(len(variables)))
    ax.set_xticklabels(var_labels, rotation=0, fontsize=11)
    ax.set_yticks(np.arange(len(regions)))
    ax.set_yticklabels(regions, fontsize=11)

    for i in range(len(regions)):
        for j in range(len(variables)):
            r = ratio_mat[i, j]
            if not np.isfinite(r):
                txt = "—"
                color = "0.5"
            else:
                if abs_df is not None:
                    row = abs_df[
                        (abs_df["region"] == regions[i])
                        & (abs_df["variable"] == variables[j])
                    ]
                    if not row.empty:
                        v1 = row.iloc[0].get(method1, np.nan)
                        v2 = row.iloc[0].get(method2, np.nan)
                        try:
                            v1 = float(v1)
                            v2 = float(v2)
                            txt = f"{r:.2f}\n({v1:.2g} / {v2:.2g})"
                        except (ValueError, TypeError):
                            txt = f"{r:.2f}"
                    else:
                        txt = f"{r:.2f}"
                else:
                    txt = f"{r:.2f}"
                # Pick text color for legibility on cmap background.
                rgba = cmap(norm(r))
                lum = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
                color = "white" if lum < 0.55 else "black"
            ax.text(j, i, txt, ha="center", va="center",
                    fontsize=10, color=color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.04)
    cbar.set_label(f"{method1} / {method2}  (lower = {method1} better)",
                   fontsize=10)

    ax.set_title(
        f"Quantile-MAE relative performance: {method1} vs {method2}",
        fontsize=12, pad=10,
    )

    fig.tight_layout()
    out_png = Path(cfg["plot_dir"]) / "qmae_cross_region_heatmap.png"
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved cross-region QMAE heatmap: %s", out_png)
    return out_png


def write_provenance(cfg, out_files, caption):
    ancestor_files = []
    for dataset in cfg.get("input_data", {}).values():
        ancestor_files.append(dataset["filename"])
    record = {
        "caption": caption,
        "statistics": ["other"],
        "domains": ["reg"],
        "plot_types": ["other"],
        "authors": ["debeire_kevin"],
        "ancestors": ancestor_files,
    }
    with ProvenanceLogger(cfg) as plog:
        for out in out_files:
            plog.log(str(out), record)


def main(cfg):
    df = gather_per_region_qmae(cfg)

    methods = cfg.get("methods", ["AFM-baseline", "PCAFM"])
    df = _coerce_numeric(df, methods)

    out_csv = save_merged_table(df, cfg)

    # Order regions in the same order the user listed them in the recipe
    # (preserves insertion order of region_labels).
    regions_in_order = [
        lbl for lbl in cfg.get("region_labels", {}).values()
        if lbl in set(df["region"].unique())
    ]
    available_vars = set(df["variable"].unique())
    variables = [v for v in QMAE_VARIABLE_ORDER if v in available_vars]
    if not variables:
        # Recipe might have produced different derived names; just use what we have.
        variables = sorted(available_vars)

    method1 = cfg.get("method1", "PCAFM")
    method2 = cfg.get("method2", "AFM-baseline")
    ratio_mat = compute_ratio_matrix(
        df, method1, method2, regions_in_order, variables,
    )
    out_png = plot_heatmap(
        ratio_mat, regions_in_order, variables, method1, method2, cfg,
        abs_df=df,
    )

    write_provenance(
        cfg,
        [out_csv, out_png],
        caption=(
            f"Cross-region quantile-MAE comparison ({method1} vs {method2}) "
            f"across {len(regions_in_order)} regions."
        ),
    )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
