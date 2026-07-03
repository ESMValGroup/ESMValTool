"""ML-based Downscaling Relative Performance Heatmap.

This diagnostic gathers metric values from ancestor diagnostic outputs
(comprehensive_metrics_summary.csv files) and creates a publication-quality
heatmap comparing the relative performance of two chosen methods.

Single-method mode (method1 is a string):
  The heatmap shows the ratio metric_method1 / metric_method2, where:
    - Green (ratio < 1): method1 is better
    - Red   (ratio > 1): method1 is worse
    - White (ratio = 1): methods are equivalent
  Rows = variables, columns = metrics, last row/column = averages.

Ablation mode (method1 is a list of strings):
  Activated automatically when more than one method1 is provided.
  The heatmap shows, for each (variable, method) cell, the mean ratio
  across all available metrics:
    cell = mean_m [ metric_m(method) / metric_m(method2) ]
  Rows = variables + Average row, columns = ablation methods.
  The colormap legend lists the metrics used to compute the average.

Special handling for metrics with non-zero optima:
  - SSR  (closer to 1 is better): |val - 1| used in ratio
  - acf_error (closer to 0): |val| used in ratio
"""

import logging
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared._base import ProvenanceLogger

logger = logging.getLogger(os.path.basename(__file__))

# ---------------------------------------------------------------------------
# Ordered variable / metric lists
# ---------------------------------------------------------------------------
VARIABLE_ORDER = ["pr", "huss", "tas", "ps", "uas", "vas"]

METRIC_ORDER = [
    "bias_mean",
    "crps_mean",
    "conservation_mean",
    "RALSD",
    "log_pdf_distance",
    "acf_error",
    "MCB",
]

CLOSER_TO_ONE_METRICS = {"ssr"}
CLOSER_TO_ZERO_METRICS = {"acf_error"}

METRIC_LABELS = {
    "bias_mean": "Bias",
    "crps_mean": "CRPS",
    "conservation_mean": "Conserv.\nError",
    "RALSD": "RALSD",
    "log_pdf_distance": "Log-PDF\nDist.",
    "acf_error": "ACF\nError",
    "MCB": "MCB",
}

VARIABLE_LABELS = {
    "pr": "Precip. (pr)",
    "huss": "Spec. Humid. (huss)",
    "tas": "Temp. (tas)",
    "ps": "Sfc. Press. (ps)",
    "uas": "U-Wind (uas)",
    "vas": "V-Wind (vas)",
}

QMAE_VARIABLE_ORDER = ["wbgt", "pr", "sfcWind"]

QMAE_VARIABLE_LABELS = {
    "wbgt": "Summer WBGT",
    "pr": "Winter Precip.",
    "sfcWind": "Fall Wind Speed",
}


# ---------------------------------------------------------------------------
# Shared utilities
# ---------------------------------------------------------------------------

def _get_provenance_record(cfg, plot_file, caption):
    ancestor_files = []
    for dataset in cfg.get("input_data", {}).values():
        ancestor_files.append(dataset["filename"])
    record = {
        "caption": caption,
        "statistics": ["other"],
        "domains": ["reg"],
        "plot_types": ["other"],
        "authors": ["debeire_kevin"],
        "references": [],
        "plot_file": plot_file,
        "ancestors": ancestor_files,
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(plot_file, record)


def _build_colormap():
    """Return the colorblind-friendly diverging colormap used by all heatmaps.

    Blue (improvement) → white (neutral) → orange (degradation), inspired by
    ColorBrewer "RdBu" but with the orange chosen to remain distinguishable
    under deuteranopia/protanopia. This replaces the previous green/red palette
    in response to reviewer R2-Fig3 (revision 1)."""
    colors_below = [
        (0.0, "#053061"), (0.3, "#2166ac"),
        (0.7, "#92c5de"), (1.0, "#f7f7f7"),
    ]
    colors_above = [
        (0.0, "#f7f7f7"), (0.3, "#fdb863"),
        (0.7, "#e08214"), (1.0, "#b35806"),
    ]

    n_steps = 256
    n_half = n_steps // 2

    def _interp(gradient, n):
        result = []
        for i in range(n):
            t = i / max(n - 1, 1)
            for k in range(len(gradient) - 1):
                t0, c0 = gradient[k]
                t1, c1 = gradient[k + 1]
                if t0 <= t <= t1:
                    frac = (t - t0) / (t1 - t0)
                    r0, g0, b0 = mcolors.to_rgb(c0)
                    r1, g1, b1 = mcolors.to_rgb(c1)
                    result.append((
                        r0 + frac * (r1 - r0),
                        g0 + frac * (g1 - g0),
                        b0 + frac * (b1 - b0),
                    ))
                    break
        return result

    all_colors = _interp(colors_below, n_half) + _interp(colors_above, n_steps - n_half)
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "GreenWhiteRed", all_colors, N=n_steps
    )
    cmap.set_bad(color="#d9d9d9")
    return cmap


def _draw_cells(ax, aug_data, row_labels, col_labels,
                cmap, norm, avg_rows, avg_cols, vmin, vmax,
                col_label_rotation=0):
    """Draw heatmap cells with text annotations.

    Parameters
    ----------
    avg_rows : set of int
        Row indices that should receive bold/thick styling (average rows).
    avg_cols : set of int
        Column indices that should receive bold/thick styling (average cols).
    col_label_rotation : float, optional
        Rotation angle in degrees for column labels. Use 0 (default) for
        horizontal labels in single-method mode and 40-45 for ablation mode
        where method names tend to be longer.
    """
    n_rows, n_cols = aug_data.shape

    for i in range(n_rows):
        for j in range(n_cols):
            val = aug_data[i, j]
            is_avg = (i in avg_rows) or (j in avg_cols)

            if np.isnan(val):
                fc = "#e8e8e8"
                text = "—"
                text_color = "#999999"
                fontweight = "normal"
                fontsize = 11
            else:
                rgba = cmap(norm(np.clip(val, vmin, vmax)))
                fc = rgba
                lum = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
                text_color = "#1a1a1a" if lum > 0.55 else "white"
                # R2-SI2: annotate cells whose numerical value exceeds the
                # colorbar range so the reader can see e.g. 2.01 explicitly
                # instead of just the saturated colour.
                off_scale = (val > vmax + 1e-9) or (val < vmin - 1e-9)
                if off_scale:
                    text = f"{val:.2f}*"
                else:
                    text = f"{val:.2f}"
                fontweight = "bold" if (is_avg or off_scale) else "medium"
                fontsize = 11 if (is_avg or off_scale) else 10.5

            lw = 1.8 if is_avg else 0.8
            ec_color = "#333333" if is_avg else "#b0b0b0"
            rect = plt.Rectangle(
                (j, n_rows - 1 - i), 1, 1,
                facecolor=fc, edgecolor=ec_color, linewidth=lw, zorder=2,
            )
            ax.add_patch(rect)
            ax.text(
                j + 0.5, n_rows - 1 - i + 0.5, text,
                ha="center", va="center",
                fontsize=fontsize, fontweight=fontweight,
                color=text_color, zorder=3,
            )

    # Separator lines above avg rows / left of avg cols
    for r in avg_rows:
        y = n_rows - r  # bottom of the separator
        ax.plot([0, n_cols], [y, y], color="#333333", linewidth=2.5, zorder=4)
    for c in avg_cols:
        ax.plot([c, c], [0, n_rows], color="#333333", linewidth=2.5, zorder=4)

    # Outer border
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(2.0)
        spine.set_color("#333333")

    n_rows_aug, n_cols_aug = aug_data.shape
    ax.set_xlim(0, n_cols_aug)
    ax.set_ylim(0, n_rows_aug)
    ax.set_aspect("equal")

    ax.set_xticks([j + 0.5 for j in range(n_cols_aug)])
    ha = "left" if col_label_rotation else "center"
    ax.set_xticklabels(col_labels, fontsize=12, fontweight="bold",
                       ha=ha, va="bottom", rotation=col_label_rotation,
                       rotation_mode="anchor")
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.tick_params(axis="x", which="both", length=0, pad=8)

    ax.set_yticks([n_rows_aug - 1 - i + 0.5 for i in range(n_rows_aug)])
    ax.set_yticklabels(row_labels, fontsize=12, fontweight="bold", ha="right")
    ax.tick_params(axis="y", which="both", length=0, pad=8)
    ax.minorticks_off()


def _add_colorbar(fig, cmap, norm, vmin, vmax, label):
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.15, 0.08, 0.70, 0.025])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cbar.set_ticks([vmin, 0.9, 0.95, 1.0, 1.05, 1.1, vmax])
    cbar.set_ticklabels(
        [f"{vmin:.2f}", "0.90", "0.95", "1.00", "1.05", "1.10", f"{vmax:.2f}"]
    )
    cbar.ax.tick_params(labelsize=9.5)
    # Avoid cbar.set_label here: when the label has multiple lines its
    # bounding box overlaps the asterisk footnote below. Place every caption
    # line as an explicit figtext so the spacing is fully under our control.
    lines = [ln for ln in label.split("\n") if ln.strip()]
    y = 0.045
    for ln in lines:
        fig.text(0.5, y, ln, ha="center", va="top", fontsize=10, color="black")
        y -= 0.022
    fig.text(
        0.5, y - 0.006,
        "Cells annotated with '*' have values beyond the colorbar range; "
        "the printed number gives the true value.",
        ha="center", va="top", fontsize=8, style="italic", color="#555555",
    )


# ---------------------------------------------------------------------------
# Data gathering
# ---------------------------------------------------------------------------

def gather_metrics_from_ancestors(cfg):
    """Gather all comprehensive_metrics_summary.csv files from ancestors."""
    ancestor_dirs = cfg.get("input_files", [])
    logger.info("Searching %d ancestor directories for metric CSVs",
                len(ancestor_dirs))

    all_dfs = []
    found_paths = set()

    for ancestor_dir in ancestor_dirs:
        ancestor_path = Path(ancestor_dir)
        for csv_path in ancestor_path.rglob("comprehensive_metrics_summary.csv"):
            csv_str = str(csv_path.resolve())
            if csv_str not in found_paths:
                logger.info("Found metrics CSV: %s", csv_path)
                all_dfs.append(pd.read_csv(csv_path))
                found_paths.add(csv_str)
        for search_dir in [ancestor_path.parent]:
            csv_file = search_dir / "comprehensive_metrics_summary.csv"
            csv_str = str(csv_file.resolve()) if csv_file.exists() else ""
            if csv_file.exists() and csv_str not in found_paths:
                logger.info("Found metrics CSV: %s", csv_file)
                all_dfs.append(pd.read_csv(csv_file))
                found_paths.add(csv_str)

    if not all_dfs:
        run_dir = Path(cfg["work_dir"]).parent.parent
        logger.info("Fallback: searching run dir: %s", run_dir)
        for csv_path in run_dir.rglob("comprehensive_metrics_summary.csv"):
            logger.info("Found metrics CSV: %s", csv_path)
            all_dfs.append(pd.read_csv(csv_path))

    if not all_dfs:
        raise FileNotFoundError(
            "No comprehensive_metrics_summary.csv files found. "
            "Ensure ancestor diagnostics have run successfully."
        )

    combined = pd.concat(all_dfs, ignore_index=True)
    combined = combined.drop_duplicates(subset=["variable", "metric"], keep="last")
    logger.info("Combined metrics table: %d rows", len(combined))
    return combined


def _safe_ratio(val1, val2, metric, epsilon=1e-10):
    """Return the transformed ratio for a single metric cell."""
    if metric in CLOSER_TO_ONE_METRICS:
        t1, t2 = abs(val1 - 1.0), abs(val2 - 1.0)
    elif metric in CLOSER_TO_ZERO_METRICS:
        t1, t2 = abs(val1), abs(val2)
    else:
        t1, t2 = abs(val1), abs(val2)
    return t1 / max(t2, epsilon)


# ---------------------------------------------------------------------------
# Single-method mode
# ---------------------------------------------------------------------------

def compute_ratio_matrix(df, method1, method2):
    """Compute the (variable x metric) ratio matrix for single-method mode."""
    available_metrics = set(df["metric"].unique())
    metrics_to_use = [m for m in METRIC_ORDER
                      if m in available_metrics and "quantile" not in m.lower()]
    vars_to_use = [v for v in VARIABLE_ORDER if v in set(df["variable"].unique())]

    logger.info("Variables: %s", vars_to_use)
    logger.info("Metrics:   %s", metrics_to_use)

    ratio_data = {}
    for var in vars_to_use:
        ratio_data[var] = {}
        for metric in metrics_to_use:
            row = df[(df["variable"] == var) & (df["metric"] == metric)]
            if row.empty:
                ratio_data[var][metric] = np.nan
                continue
            row = row.iloc[0]
            try:
                v1 = float(str(row.get(method1, "N/A")))
                v2 = float(str(row.get(method2, "N/A")))
            except (ValueError, TypeError):
                ratio_data[var][metric] = np.nan
                continue
            ratio_data[var][metric] = _safe_ratio(v1, v2, metric)

    ratio_df = pd.DataFrame(ratio_data).T
    ratio_df = ratio_df.reindex(index=vars_to_use, columns=metrics_to_use)
    return ratio_df, vars_to_use, metrics_to_use


def plot_heatmap(ratio_df, vars_list, metrics_list, method1, method2, cfg):
    """Create the standard (variable x metric) relative-performance heatmap."""
    vmin = cfg.get("heatmap_vmin", 0.80)
    vmax = cfg.get("heatmap_vmax", 1.20)
    cmap = _build_colormap()
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=1.0, vmax=vmax)

    data = ratio_df.values.copy()
    n_vars, n_metrics = data.shape

    var_avg = np.nanmean(data, axis=1)
    metric_avg = np.nanmean(data, axis=0)
    overall_avg = np.nanmean(data)

    aug_data = np.full((n_vars + 1, n_metrics + 1), np.nan)
    aug_data[:n_vars, :n_metrics] = data
    aug_data[:n_vars, n_metrics] = var_avg
    aug_data[n_vars, :n_metrics] = metric_avg
    aug_data[n_vars, n_metrics] = overall_avg

    row_labels = [VARIABLE_LABELS.get(v, v) for v in vars_list] + ["Average"]
    col_labels = [METRIC_LABELS.get(m, m) for m in metrics_list] + ["Average"]

    n_rows, n_cols = aug_data.shape
    fig, ax = plt.subplots(figsize=(n_cols * 1.35 + 3.2, n_rows * 0.72 + 2.8))

    _draw_cells(
        ax, aug_data, row_labels, col_labels, cmap, norm,
        avg_rows={n_rows - 1}, avg_cols={n_cols - 1}, vmin=vmin, vmax=vmax,
    )

    region_name = cfg.get("region_name", "")
    title = (
        f"{region_name}\nRelative Performance of {method1}  vs  {method2}"
        if region_name
        else f"Relative Performance of {method1}  vs  {method2}"
    )
    ax.set_title(title, fontsize=18, fontweight="bold", pad=28, loc="center")

    _add_colorbar(
        fig, cmap, norm, vmin, vmax,
        label=(
            f"Metric Ratio  {method1} / {method2}\n"
            f"(Ratio < 1: {method1} is better  |  Ratio > 1: {method2} is better)"
        ),
    )

    plt.tight_layout(rect=[0.0, 0.10, 1.0, 1.0])
    plot_file = os.path.join(
        cfg["plot_dir"],
        f"relative_performance_heatmap_{method1}_vs_{method2}.png",
    )
    plt.savefig(plot_file, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()

    caption = (
        f"Relative performance heatmap comparing {method1} vs {method2}. "
        f"Values < 1 (green) indicate {method1} outperforms {method2}."
    )
    _get_provenance_record(cfg, plot_file, caption)
    logger.info("Saved heatmap: %s", plot_file)
    return plot_file


def save_ratio_table(ratio_df, vars_list, metrics_list, method1, method2, cfg):
    """Save the ratio matrix as CSV and formatted text."""
    data = ratio_df.copy()
    data["Average"] = data.mean(axis=1, skipna=True)
    avg_row = data.mean(axis=0, skipna=True)
    avg_row.name = "Average"
    data = pd.concat([data, avg_row.to_frame().T])

    csv_file = os.path.join(
        cfg["work_dir"],
        f"relative_performance_{method1}_vs_{method2}.csv",
    )
    data.to_csv(csv_file)

    txt_file = csv_file.replace(".csv", ".txt")
    with open(txt_file, "w") as f:
        f.write("=" * 80 + "\n")
        f.write(f"RELATIVE PERFORMANCE: {method1} vs {method2}\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Ratio = metric({method1}) / metric({method2})\n")
        f.write("Ratio < 1 → method1 is better\n")
        f.write("Ratio > 1 → method2 is better\n\n")
        f.write("-" * 80 + "\n\n")
        f.write(data.to_string(float_format="%.4f"))
        f.write("\n")
    logger.info("Saved ratio table: %s", csv_file)


# ---------------------------------------------------------------------------
# Ablation mode  (method1 is a list of methods)
# ---------------------------------------------------------------------------

def compute_ablation_matrix(df, methods, method2):
    """Compute a (variable x method) matrix of metric-averaged ratios.

    For each (variable, method) cell the value is the mean of the per-metric
    ratios across all available metrics:

        cell(var, m) = mean_k [ ratio_k(m, method2) ]

    where ratio_k is computed identically to the single-method mode.

    Parameters
    ----------
    df : pd.DataFrame
        Combined metrics DataFrame.
    methods : list of str
        Ablation method names (columns of the heatmap).
    method2 : str
        Reference method (denominator).

    Returns
    -------
    avg_ratio_df : pd.DataFrame
        Shape (n_vars, n_methods); values are mean ratios (NaN when unavailable).
    vars_list : list of str
        Ordered variable names present in the data.
    metrics_used : list of str
        Metric names that were used to compute the averages.
    """
    available_metrics = set(df["metric"].unique())
    metrics_used = [m for m in METRIC_ORDER
                    if m in available_metrics and "quantile" not in m.lower()]
    vars_list = [v for v in VARIABLE_ORDER if v in set(df["variable"].unique())]

    logger.info("Ablation variables: %s", vars_list)
    logger.info("Ablation metrics:   %s", metrics_used)

    # Build a 3-D structure: per_metric_ratios[var][method] = list of ratios
    data = {}
    for var in vars_list:
        data[var] = {}
        for method in methods:
            per_metric = []
            for metric in metrics_used:
                row = df[(df["variable"] == var) & (df["metric"] == metric)]
                if row.empty:
                    continue
                row = row.iloc[0]
                try:
                    v1 = float(str(row.get(method, "N/A")))
                    v2 = float(str(row.get(method2, "N/A")))
                except (ValueError, TypeError):
                    continue
                per_metric.append(_safe_ratio(v1, v2, metric))
            data[var][method] = np.nanmean(per_metric) if per_metric else np.nan

    avg_ratio_df = pd.DataFrame(data).T          # rows=vars, cols=methods
    avg_ratio_df = avg_ratio_df.reindex(index=vars_list, columns=methods)
    return avg_ratio_df, vars_list, metrics_used


def _plot_ablation_heatmap_raw(data, row_labels, col_labels, method2,
                               metrics_label, title_suffix, filename_prefix,
                               cfg):
    """Low-level renderer for an ablation-style heatmap.

    Accepts pre-labelled numpy data directly so it can be shared between the
    main ablation heatmap (standard variables) and the quantile MAE ablation
    heatmap (derived variables) without duplicating any rendering logic.

    An Average row is appended automatically.

    Parameters
    ----------
    data : np.ndarray
        Shape (n_rows, n_cols); NaN for missing cells.
    row_labels : list of str
        Human-readable row labels (variables).
    col_labels : list of str
        Human-readable column labels (method names).
    method2 : str
        Reference method name used in axis/colorbar labels.
    metrics_label : str
        Short description of what was averaged (shown in colorbar).
    title_suffix : str
        Appended to the figure title, e.g. "Ablation Study" or "Quantile MAE".
    filename_prefix : str
        Prefix for the saved file, e.g. "ablation_heatmap" or "ablation_qmae".
    cfg : dict
        ESMValTool configuration dictionary.

    Returns
    -------
    str
        Path to the saved plot file.
    """
    vmin = cfg.get("heatmap_vmin", 0.80)
    vmax = cfg.get("heatmap_vmax", 1.20)
    cmap = _build_colormap()
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=1.0, vmax=vmax)

    n_vars, n_methods = data.shape

    # Append Average row
    method_avg = np.nanmean(data, axis=0)
    aug_data = np.full((n_vars + 1, n_methods), np.nan)
    aug_data[:n_vars, :] = data
    aug_data[n_vars, :] = method_avg
    all_row_labels = list(row_labels) + ["Average"]

    n_rows, n_cols = aug_data.shape
    cell_w = max(1.6, 10.0 / n_cols)
    cell_h = 0.72
    fig, ax = plt.subplots(
        figsize=(n_cols * cell_w + 3.2, n_rows * cell_h + 2.8)
    )

    _draw_cells(
        ax, aug_data, all_row_labels, col_labels, cmap, norm,
        avg_rows={n_rows - 1}, avg_cols=set(), vmin=vmin, vmax=vmax,
        col_label_rotation=45,
    )

    region_name = cfg.get("region_name", "")
    title_base = f"Relative Performance vs {method2}  |  {title_suffix}"
    title = f"{region_name}\n{title_base}" if region_name else title_base
    ax.set_title(title, fontsize=18, fontweight="bold", pad=60, loc="center")

    cbar_label = (
        f"{metrics_label} ratio  /  {method2}\n"
        f"(Ratio < 1: column method is better  |  Ratio > 1: {method2} is better)"
    )
    _add_colorbar(fig, cmap, norm, vmin, vmax, label=cbar_label)

    plt.tight_layout(rect=[0.0, 0.12, 1.0, 1.0])

    methods_str = "_".join(c.replace(" ", "-") for c in col_labels)
    plot_file = os.path.join(
        cfg["plot_dir"],
        f"{filename_prefix}_{methods_str}_vs_{method2}.png",
    )
    plt.savefig(plot_file, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()

    caption = (
        f"{title_suffix} ablation heatmap. Columns are ablation methods; "
        f"rows are variables. Cell values are {metrics_label} ratios "
        f"(method / {method2}). Green cells (ratio < 1) indicate "
        f"improvement over {method2}."
    )
    _get_provenance_record(cfg, plot_file, caption)
    logger.info("Saved ablation heatmap: %s", plot_file)
    return plot_file


def plot_ablation_heatmap(avg_ratio_df, vars_list, methods, metrics_used,
                          method2, cfg):
    """Create the ablation heatmap (variables x ablation methods).

    Rows = variables + Average row.
    Columns = one per method in the method1 list.
    Cell values = mean ratio across all standard metrics.

    Parameters
    ----------
    avg_ratio_df : pd.DataFrame
        Shape (n_vars, n_methods); produced by compute_ablation_matrix.
    vars_list : list of str
    methods : list of str
        Ablation method names.
    metrics_used : list of str
        Metrics averaged to produce cell values (shown in colorbar label).
    method2 : str
        Reference method name.
    cfg : dict
    """
    metric_label_list = ", ".join(
        METRIC_LABELS.get(m, m).replace("\n", " ") for m in metrics_used
    )
    row_labels = [VARIABLE_LABELS.get(v, v) for v in vars_list]

    return _plot_ablation_heatmap_raw(
        data=avg_ratio_df.values.copy(),
        row_labels=row_labels,
        col_labels=list(methods),
        method2=method2,
        metrics_label=f"Mean ratio ({metric_label_list})",
        title_suffix="Ablation Study",
        filename_prefix="ablation_heatmap",
        cfg=cfg,
    )


def save_ablation_table(avg_ratio_df, vars_list, methods, metrics_used,
                        method2, cfg):
    """Save the ablation ratio table as CSV."""
    data = avg_ratio_df.copy()
    avg_row = data.mean(axis=0, skipna=True)
    avg_row.name = "Average"
    data = pd.concat([data, avg_row.to_frame().T])

    csv_file = os.path.join(
        cfg["work_dir"],
        f"ablation_performance_vs_{method2}.csv",
    )
    data.to_csv(csv_file)

    txt_file = csv_file.replace(".csv", ".txt")
    with open(txt_file, "w") as f:
        f.write("=" * 80 + "\n")
        f.write(f"ABLATION STUDY: methods vs {method2}\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Cell = mean ratio across metrics: {', '.join(metrics_used)}\n")
        f.write("Ratio < 1 → column method is better than reference\n\n")
        f.write("-" * 80 + "\n\n")
        f.write(data.to_string(float_format="%.4f"))
        f.write("\n")
    logger.info("Saved ablation table: %s", csv_file)


# ---------------------------------------------------------------------------
# Quantile MAE heatmap (unchanged from original)
# ---------------------------------------------------------------------------

def gather_quantile_mae_data(cfg):
    ancestor_dirs = cfg.get("input_files", [])
    all_dfs = []
    found_paths = set()

    for ancestor_dir in ancestor_dirs:
        ancestor_path = Path(ancestor_dir)
        for csv_path in ancestor_path.rglob("quantile_mae_table.csv"):
            csv_str = str(csv_path.resolve())
            if csv_str not in found_paths:
                logger.info("Found quantile MAE CSV: %s", csv_path)
                all_dfs.append(pd.read_csv(csv_path))
                found_paths.add(csv_str)
        for csv_path in ancestor_path.rglob("comprehensive_metrics_summary.csv"):
            csv_str = str(csv_path.resolve())
            if csv_str not in found_paths:
                df = pd.read_csv(csv_path)
                qmae_rows = df[df["metric"] == "avg_quantile_mae"]
                if not qmae_rows.empty:
                    all_dfs.append(qmae_rows)
                    found_paths.add(csv_str)
        for search_dir in [ancestor_path.parent]:
            for fname in ["quantile_mae_table.csv", "comprehensive_metrics_summary.csv"]:
                csv_file = search_dir / fname
                csv_str = str(csv_file.resolve()) if csv_file.exists() else ""
                if csv_file.exists() and csv_str not in found_paths:
                    df = pd.read_csv(csv_file)
                    if fname == "comprehensive_metrics_summary.csv":
                        df = df[df["metric"] == "avg_quantile_mae"]
                    if not df.empty:
                        all_dfs.append(df)
                        found_paths.add(csv_str)

    if not all_dfs:
        run_dir = Path(cfg["work_dir"]).parent.parent
        for csv_path in run_dir.rglob("quantile_mae_table.csv"):
            all_dfs.append(pd.read_csv(csv_path))

    if not all_dfs:
        logger.warning("No quantile MAE data found.")
        return None

    combined = pd.concat(all_dfs, ignore_index=True)
    combined = combined.drop_duplicates(subset=["variable"], keep="last")
    if "metric" not in combined.columns:
        combined["metric"] = "avg_quantile_mae"
    return combined


def compute_qmae_ablation_matrix(df, methods, method2):
    """Compute a (qmae_variable x method) ratio matrix for ablation mode.

    Mirrors compute_qmae_ratio but accepts a list of methods, producing
    one column per method instead of a single scalar per variable.

    Parameters
    ----------
    df : pd.DataFrame
        Quantile MAE DataFrame with one row per derived variable.
    methods : list of str
        Ablation method names (columns of the heatmap).
    method2 : str
        Reference method (denominator).

    Returns
    -------
    avg_ratio_df : pd.DataFrame
        Shape (n_qmae_vars, n_methods).
    vars_list : list of str
        Ordered derived variable names present in the data.
    """
    available_vars = set(df["variable"].unique())
    vars_list = [v for v in QMAE_VARIABLE_ORDER if v in available_vars]
    epsilon = 1e-10

    data = {}
    for var in vars_list:
        row = df[df["variable"] == var]
        if row.empty:
            data[var] = {m: np.nan for m in methods}
            continue
        row = row.iloc[0]
        try:
            v2 = abs(float(str(row.get(method2, "N/A"))))
        except (ValueError, TypeError):
            data[var] = {m: np.nan for m in methods}
            continue
        data[var] = {}
        for method in methods:
            try:
                v1 = abs(float(str(row.get(method, "N/A"))))
            except (ValueError, TypeError):
                data[var][method] = np.nan
                continue
            data[var][method] = v1 / max(v2, epsilon)

    avg_ratio_df = pd.DataFrame(data).T          # rows=vars, cols=methods
    avg_ratio_df = avg_ratio_df.reindex(index=vars_list, columns=methods)
    return avg_ratio_df, vars_list


def compute_qmae_ratio(df, method1, method2):
    available_vars = set(df["variable"].unique())
    vars_to_use = [v for v in QMAE_VARIABLE_ORDER if v in available_vars]
    epsilon = 1e-10
    ratios = {}
    for var in vars_to_use:
        row = df[df["variable"] == var]
        if row.empty:
            continue
        row = row.iloc[0]
        try:
            v1 = abs(float(str(row.get(method1, "N/A"))))
            v2 = abs(float(str(row.get(method2, "N/A"))))
        except (ValueError, TypeError):
            ratios[var] = np.nan
            continue
        ratios[var] = v1 / max(v2, epsilon)
    return ratios, vars_to_use


def plot_quantile_mae_heatmap(ratios, vars_list, method1, method2, cfg):
    vmin = cfg.get("heatmap_vmin", 0.80)
    vmax = cfg.get("heatmap_vmax", 1.20)
    cmap = _build_colormap()
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=1.0, vmax=vmax)

    n_vars = len(vars_list)
    if n_vars == 0:
        logger.warning("No quantile MAE variables to plot.")
        return

    values = np.array([ratios.get(v, np.nan) for v in vars_list])
    aug_data = np.full((n_vars + 1, 1), np.nan)
    aug_data[:n_vars, 0] = values
    aug_data[n_vars, 0] = np.nanmean(values)

    row_labels = [QMAE_VARIABLE_LABELS.get(v, v) for v in vars_list] + ["Average"]
    col_labels = ["Quantile MAE\nRatio"]

    n_rows, n_cols = aug_data.shape
    fig, ax = plt.subplots(figsize=(n_cols * 2.0 + 3.6, n_rows * 0.80 + 2.6))

    _draw_cells(
        ax, aug_data, row_labels, col_labels, cmap, norm,
        avg_rows={n_rows - 1}, avg_cols=set(), vmin=vmin, vmax=vmax,
    )

    region_name = cfg.get("region_name", "")
    title = (
        f"{region_name}\nQuantile MAE relative performance:  {method1}  vs  {method2}"
        if region_name
        else f"Quantile MAE relative performance:  {method1}  vs  {method2}"
    )
    ax.set_title(title, fontsize=18, fontweight="bold", pad=28, loc="center")

    _add_colorbar(
        fig, cmap, norm, vmin, vmax,
        label=(
            f"Metric Ratio  {method1} / {method2}\n"
            f"(Ratio < 1: {method1} is better  |  Ratio > 1: {method2} is better)"
        ),
    )

    plt.tight_layout(rect=[0.0, 0.13, 1.0, 1.0])
    plot_file = os.path.join(
        cfg["plot_dir"],
        f"quantile_mae_heatmap_{method1}_vs_{method2}.png",
    )
    plt.savefig(plot_file, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()

    _get_provenance_record(
        cfg, plot_file,
        f"Quantile MAE relative performance heatmap: {method1} vs {method2}.",
    )
    logger.info("Saved quantile MAE heatmap: %s", plot_file)
    return plot_file


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(cfg):
    logger.setLevel(cfg.get("log_level", "INFO").upper())

    method1_cfg = cfg.get("method1")
    method2 = cfg.get("method2")

    if not method1_cfg or not method2:
        raise ValueError(
            "Both 'method1' and 'method2' must be specified in the recipe "
            "script configuration."
        )

    # Normalise method1 to a list regardless of how it was specified
    if isinstance(method1_cfg, str):
        methods1 = [method1_cfg]
    else:
        methods1 = list(method1_cfg)

    ablation_mode = len(methods1) > 1
    logger.info(
        "Mode: %s | method1=%s | method2=%s",
        "ablation" if ablation_mode else "single",
        methods1, method2,
    )

    # Gather metrics
    combined_df = gather_metrics_from_ancestors(cfg)
    available_cols = [
        c for c in combined_df.columns
        if c not in ("variable", "analysis_type", "metric")
    ]
    logger.info("Available methods in data: %s", available_cols)

    # Validate all requested methods are present
    for m in methods1 + [method2]:
        if m not in combined_df.columns:
            raise ValueError(
                f"Method '{m}' not found in metrics data. "
                f"Available: {available_cols}"
            )

    if ablation_mode:
        # ------------------------------------------------------------------
        # Ablation path: one column per method, rows = variables + Average
        # ------------------------------------------------------------------
        avg_ratio_df, vars_list, metrics_used = compute_ablation_matrix(
            combined_df, methods1, method2
        )
        save_ablation_table(
            avg_ratio_df, vars_list, methods1, metrics_used, method2, cfg
        )
        plot_ablation_heatmap(
            avg_ratio_df, vars_list, methods1, metrics_used, method2, cfg
        )

        # Quantile MAE ablation heatmap (same layout: qmae vars x methods)
        qmae_df = gather_quantile_mae_data(cfg)
        if qmae_df is not None and all(
            m in qmae_df.columns for m in methods1 + [method2]
        ):
            qmae_ratio_df, qmae_vars = compute_qmae_ablation_matrix(
                qmae_df, methods1, method2
            )
            if qmae_vars:
                row_labels = [QMAE_VARIABLE_LABELS.get(v, v) for v in qmae_vars]
                _plot_ablation_heatmap_raw(
                    data=qmae_ratio_df.values.copy(),
                    row_labels=row_labels,
                    col_labels=list(methods1),
                    method2=method2,
                    metrics_label="Avg. Quantile MAE",
                    title_suffix="Quantile MAE  |  Ablation Study",
                    filename_prefix="ablation_qmae",
                    cfg=cfg,
                )
            else:
                logger.warning("No quantile MAE variables found for ablation heatmap.")
        else:
            logger.warning(
                "Quantile MAE data not available or methods not found. "
                "Skipping quantile MAE ablation heatmap."
            )

    else:
        # ------------------------------------------------------------------
        # Standard single-method path (unchanged behaviour)
        # ------------------------------------------------------------------
        method1 = methods1[0]

        ratio_df, vars_list, metrics_list = compute_ratio_matrix(
            combined_df, method1, method2
        )
        logger.info(
            "Ratio matrix:\n%s", ratio_df.to_string(float_format="%.4f")
        )
        save_ratio_table(
            ratio_df, vars_list, metrics_list, method1, method2, cfg
        )
        plot_heatmap(ratio_df, vars_list, metrics_list, method1, method2, cfg)

        # Quantile MAE heatmap
        qmae_df = gather_quantile_mae_data(cfg)
        if (qmae_df is not None
                and method1 in qmae_df.columns
                and method2 in qmae_df.columns):
            qmae_ratios, qmae_vars = compute_qmae_ratio(qmae_df, method1, method2)
            if qmae_vars:
                plot_quantile_mae_heatmap(
                    qmae_ratios, qmae_vars, method1, method2, cfg
                )
            else:
                logger.warning("No quantile MAE variables found for heatmap.")
        else:
            logger.warning(
                "Quantile MAE data not available or methods not found. "
                "Skipping quantile MAE heatmap."
            )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)