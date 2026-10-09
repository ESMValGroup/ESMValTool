"""Pareto-front diagnostic for the conservation / generative-fidelity trade-off.

Reads ``comprehensive_metrics_summary.csv`` from ancestor diagnostics (one per
region) and/or from explicit paths given in the recipe, then plots each trained
variant as a point in (conservation error, generative fidelity) space.  Both
axes are "lower is better", so the attainable frontier runs toward the
lower-left and the non-dominated set is drawn as a step line.

It shows how each gradient-surgery and loss-weighting choice trades
conservation error against generative fidelity, metrics that the evaluation
diagnostics otherwise report separately.

Two ways to supply data:

1. ``ancestors`` + ``region_labels`` -- the normal ESMValTool flow, matching
   ``ml_downscaling_qmae_cross_region.py``.
2. ``external_metrics_csv`` -- a mapping of region label to an absolute path to
   a ``comprehensive_metrics_summary.csv`` written by a *previous* run.  This
   lets the figure be built from completed ablation runs without re-running the
   (expensive) evaluation diagnostics.

Both may be combined; duplicate (region, variable, metric) rows are dropped.

Colour encodes method *family* only (three slots, the maximum that validates
for all-pairs colour separation in a scatter); marker shape encodes the variant
within a family; and every point carries a direct text label, so identity is
never conveyed by colour alone.
"""

import logging
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import Bbox

from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared._base import ProvenanceLogger

logger = logging.getLogger(os.path.basename(__file__))

METRICS_CSV = "comprehensive_metrics_summary.csv"

# Categorical slots 1-3 of the validated default palette.  Only the first three
# slots clear the all-pairs colour-separation floors, which is the pairlist that
# applies to scatter plots -- hence three families, not one colour per method.
# The families split on *how the physics loss enters training*, which is the
# distinction the figure exists to make.  Lumping the naive weighted sum in with
# the unconstrained baseline would be wrong: it does carry L_phys, it just adds
# it to L_base without resolving the gradient conflict.
FAMILY_COLORS = {
    "pcafm": "#2a78d6",   # slot 1, blue   -- L_phys + ConFIG gradient surgery
    "afm": "#eb6834",     # slot 2, orange -- no L_phys at all
    "naive": "#1baf7a",   # slot 3, aqua   -- L_phys, plain weighted sum
}

FAMILY_LABELS = {
    "pcafm": "PC-AFM: physics loss + ConFIG",
    "afm": "AFM baseline: no physics loss",
    "naive": "Physics loss, naive weighted sum",
}

# Marker shapes cycle within a family; identity is carried by the direct label,
# so a repeated shape across families is harmless.
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*", "h", "<", ">", "p"]

# Text ink, never a series colour.
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID_COLOR = "#d8d8d4"

VAR_LABELS = {
    "pr": "Precipitation",
    "huss": "Specific humidity",
    "tas": "Near-surface temperature",
    "ps": "Surface pressure",
}

VAR_UNITS = {
    "pr": "mm h$^{-1}$",
    "huss": "g kg$^{-1}$",
    "tas": "$^{\\circ}$C",
    "ps": "hPa",
}

METRIC_LABELS = {
    "crps_mean": "CRPS",
    "RALSD": "RALSD (dB)",
    "mae_mean": "MAE",
    "log_pdf_distance": "Log-PDF distance",
    "avg_quantile_mae": "Avg. quantile MAE",
    "MCB": "MCB",
    "conservation_mean": "Conservation error",
}


# ---------------------------------------------------------------------------
# Data gathering
# ---------------------------------------------------------------------------

def _identify_region(ancestor_dir, region_labels):
    """Return the region label whose diagnostic key appears in *ancestor_dir*."""
    path_str = str(ancestor_dir)
    for diag_key, label in region_labels.items():
        if diag_key in path_str:
            return label
    return None


def gather_metrics(cfg):
    """Collect per-region metric tables from ancestors and/or explicit paths."""
    frames = []
    seen = set()

    region_labels = cfg.get("region_labels", {})
    for ancestor_dir in cfg.get("input_files", []):
        region = _identify_region(ancestor_dir, region_labels)
        if region is None:
            logger.debug("Ancestor %s matches no region label; skipping.",
                         ancestor_dir)
            continue
        for csv_path in Path(ancestor_dir).rglob(METRICS_CSV):
            key = str(csv_path.resolve())
            if key in seen:
                continue
            seen.add(key)
            frame = pd.read_csv(csv_path)
            frame.insert(0, "region", region)
            frames.append(frame)
            logger.info("Loaded %s -> region=%s (%d rows)",
                        csv_path, region, len(frame))

    for region, csv_path in cfg.get("external_metrics_csv", {}).items():
        csv_path = Path(csv_path)
        if not csv_path.exists():
            raise FileNotFoundError(
                f"external_metrics_csv entry for '{region}' does not exist: "
                f"{csv_path}"
            )
        key = str(csv_path.resolve())
        if key in seen:
            logger.info("Skipping duplicate external CSV %s", csv_path)
            continue
        seen.add(key)
        frame = pd.read_csv(csv_path)
        frame.insert(0, "region", region)
        frames.append(frame)
        logger.info("Loaded external %s -> region=%s (%d rows)",
                    csv_path, region, len(frame))

    if not frames:
        raise RuntimeError(
            f"No {METRICS_CSV} found. Provide `ancestors` together with "
            "`region_labels`, or `external_metrics_csv`."
        )

    combined = pd.concat(frames, ignore_index=True)
    # Dedupe *within* a region only -- the same (variable, metric) legitimately
    # recurs across regions.
    combined = combined.drop_duplicates(
        subset=["region", "variable", "metric"], keep="last"
    )
    logger.info("Combined metrics table: %d rows, regions=%s",
                len(combined), sorted(combined["region"].unique()))
    return combined


def _lookup(df, region, variable, metric, method):
    """Return a float metric value, or NaN if absent/non-numeric."""
    row = df[
        (df["region"] == region)
        & (df["variable"] == variable)
        & (df["metric"] == metric)
    ]
    if row.empty or method not in row.columns:
        return np.nan
    try:
        return float(row.iloc[0][method])
    except (ValueError, TypeError):
        return np.nan


def build_points(df, cfg):
    """Assemble a long-form table of plottable (x, y) points."""
    methods = cfg["methods"]
    x_metric = cfg.get("x_metric", "conservation_mean")
    y_metric = cfg["y_metric"]
    variables = cfg.get("variables", ["pr", "huss"])

    regions = [
        lbl for lbl in _region_order(cfg)
        if lbl in set(df["region"].unique())
    ]

    rows = []
    for region in regions:
        for variable in variables:
            for method in methods:
                x_val = _lookup(df, region, variable, x_metric, method)
                y_val = _lookup(df, region, variable, y_metric, method)
                if not (np.isfinite(x_val) and np.isfinite(y_val)):
                    logger.debug(
                        "Missing %s/%s for %s in %s (x=%s, y=%s)",
                        x_metric, y_metric, method, region, x_val, y_val,
                    )
                    continue
                rows.append({
                    "region": region,
                    "variable": variable,
                    "method": method,
                    "x_metric": x_metric,
                    "x": x_val,
                    "y_metric": y_metric,
                    "y": y_val,
                })

    if not rows:
        raise RuntimeError(
            f"No finite (x={x_metric}, y={y_metric}) pairs found for methods "
            f"{methods}. Check that the metric names and method column names "
            "match the ancestor CSVs."
        )
    return pd.DataFrame(rows), regions, variables


def _region_order(cfg):
    """Preferred region ordering: recipe order, ancestors before externals."""
    order = list(cfg.get("region_labels", {}).values())
    for region in cfg.get("external_metrics_csv", {}):
        if region not in order:
            order.append(region)
    return order


# ---------------------------------------------------------------------------
# Pareto logic
# ---------------------------------------------------------------------------

def non_dominated(points):
    """Return a boolean mask of non-dominated points (both axes minimised).

    A point is dominated when another point is at least as good on both axes
    and strictly better on at least one.
    """
    xs = points["x"].to_numpy()
    ys = points["y"].to_numpy()
    mask = np.ones(len(points), dtype=bool)
    for i in range(len(points)):
        better_or_equal = (xs <= xs[i]) & (ys <= ys[i])
        strictly_better = (xs < xs[i]) | (ys < ys[i])
        if np.any(better_or_equal & strictly_better):
            mask[i] = False
    return mask


# ---------------------------------------------------------------------------
# Styling helpers
# ---------------------------------------------------------------------------

def _family_of(method, cfg):
    """Map a method name to one of the three colour families."""
    explicit = cfg.get("method_families", {})
    if method in explicit:
        return explicit[method]
    lowered = method.lower()
    if lowered.startswith("pcafm") or lowered.startswith("pc-afm"):
        return "pcafm"
    return "afm"


def _style_map(methods, cfg):
    """Assign a (colour, marker) pair to each method, stable across panels."""
    styles = {}
    per_family_count = {}
    for method in methods:
        family = _family_of(method, cfg)
        index = per_family_count.get(family, 0)
        per_family_count[family] = index + 1
        styles[method] = {
            "color": FAMILY_COLORS.get(family, FAMILY_COLORS["afm"]),
            "marker": MARKERS[index % len(MARKERS)],
            "family": family,
        }
    return styles


# Candidate label positions, in offset points from the marker, tried in order of
# preference.  Diagonals first (they read as attached to the mark without sitting
# on the gridlines), then the cardinals, then a wider ring for crowded panels.
_LABEL_CANDIDATES = [
    (8, 8, "left", "bottom"), (-8, 8, "right", "bottom"),
    (8, -8, "left", "top"), (-8, -8, "right", "top"),
    (0, 11, "center", "bottom"), (0, -11, "center", "top"),
    (12, 0, "left", "center"), (-12, 0, "right", "center"),
    (20, 16, "left", "bottom"), (-20, 16, "right", "bottom"),
    (20, -16, "left", "top"), (-20, -16, "right", "top"),
    (0, 24, "center", "bottom"), (0, -24, "center", "top"),
]

_HA_SHIFT = {"left": 0.0, "center": 0.5, "right": 1.0}
_VA_SHIFT = {"bottom": 0.0, "center": 0.5, "top": 1.0}


def _place_labels(fig, ax, items):
    """Draw direct point labels, choosing offsets that do not collide.

    The previous version pushed each label away from the panel centre with a
    fixed 7-point offset, which collided badly whenever two variants had similar
    scores -- exactly the case this figure is built to show.  Here each label's
    rendered size is measured once (it does not depend on the offset), then every
    candidate offset is scored analytically against the already-placed labels,
    the marker keep-out boxes, and the axes rectangle.  Lowest penalty wins.

    `items` is a sequence of (x_data, y_data, text).
    """
    if not items:
        return
    renderer = fig.canvas.get_renderer()
    ax_bbox = ax.get_window_extent(renderer=renderer)
    px_per_point = fig.dpi / 72.0

    # Markers are obstacles too, so a label never lands on another point.
    obstacles = []
    for x_val, y_val, _ in items:
        px, py = ax.transData.transform((x_val, y_val))
        obstacles.append(Bbox.from_bounds(px - 8, py - 8, 16, 16))

    for x_val, y_val, text in items:
        probe = ax.annotate(
            text, xy=(x_val, y_val), xytext=(0, 0),
            textcoords="offset points", fontsize=7.5, ha="left", va="bottom",
        )
        extent = probe.get_window_extent(renderer=renderer)
        width, height = extent.width, extent.height
        probe.remove()

        px, py = ax.transData.transform((x_val, y_val))
        best, best_penalty = _LABEL_CANDIDATES[0], None

        for dx, dy, hal, val in _LABEL_CANDIDATES:
            x_0 = px + dx * px_per_point - _HA_SHIFT[hal] * width
            y_0 = py + dy * px_per_point - _VA_SHIFT[val] * height
            box = Bbox.from_bounds(x_0, y_0, width, height)

            overlap = sum(_overlap_area(box, other) for other in obstacles)
            # Spilling outside the panel is worse than a small overlap: an
            # out-of-axes label collides with the neighbouring panel instead.
            outside = (width * height) - _overlap_area(box, ax_bbox)
            penalty = overlap + 3.0 * outside

            if best_penalty is None or penalty < best_penalty:
                best, best_penalty = (dx, dy, hal, val), penalty
            if penalty <= 0.0:
                break

        dx, dy, hal, val = best
        ax.annotate(
            text, xy=(x_val, y_val), xytext=(dx, dy),
            textcoords="offset points", fontsize=7.5, color=TEXT_PRIMARY,
            ha=hal, va=val, zorder=4,
        )
        x_0 = px + dx * px_per_point - _HA_SHIFT[hal] * width
        y_0 = py + dy * px_per_point - _VA_SHIFT[val] * height
        obstacles.append(Bbox.from_bounds(x_0, y_0, width, height))


def _overlap_area(box_a, box_b):
    """Area of the intersection of two bboxes, in square pixels."""
    dx = min(box_a.x1, box_b.x1) - max(box_a.x0, box_b.x0)
    dy = min(box_a.y1, box_b.y1) - max(box_a.y0, box_b.y0)
    return dx * dy if (dx > 0 and dy > 0) else 0.0


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_pareto(points, regions, variables, cfg):
    """Render the (variable x region) grid of Pareto scatter panels."""
    methods = [m for m in cfg["methods"]
               if m in set(points["method"].unique())]
    styles = _style_map(methods, cfg)
    labels = cfg.get("method_labels", {})
    y_metric = cfg["y_metric"]
    x_metric = cfg.get("x_metric", "conservation_mean")
    annotate = cfg.get("annotate_points", True)

    n_rows, n_cols = len(variables), len(regions)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(4.1 * n_cols, 3.7 * n_rows),
        squeeze=False,
    )

    frontier_rows = []
    panel_labels = []

    for i, variable in enumerate(variables):
        for j, region in enumerate(regions):
            ax = axes[i][j]
            panel = points[
                (points["variable"] == variable) & (points["region"] == region)
            ].reset_index(drop=True)

            label_items = []

            if panel.empty:
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                        ha="center", va="center", color=TEXT_SECONDARY)
                ax.set_xticks([])
                ax.set_yticks([])
                continue

            mask = non_dominated(panel)
            panel = panel.assign(non_dominated=mask)
            frontier_rows.append(panel)

            # Frontier step line, drawn under the marks.
            front = panel[mask].sort_values("x")
            if len(front) > 1:
                ax.step(
                    front["x"], front["y"], where="post",
                    color=TEXT_SECONDARY, linewidth=1.2,
                    linestyle="--", alpha=0.7, zorder=1,
                )

            for _, row in panel.iterrows():
                style = styles[row["method"]]
                on_front = bool(row["non_dominated"])
                ax.scatter(
                    row["x"], row["y"],
                    s=95 if on_front else 70,
                    c=style["color"],
                    marker=style["marker"],
                    edgecolors="#ffffff",
                    linewidths=1.6 if on_front else 1.0,
                    zorder=3,
                    alpha=1.0 if on_front else 0.75,
                )
                if annotate:
                    label_items.append((
                        row["x"], row["y"],
                        labels.get(row["method"], row["method"]),
                    ))

            # Headroom so no marker sits on a spine and labels have somewhere to
            # go.  Must be set before the label pass, which works in pixels.
            ax.margins(0.18, 0.20)
            panel_labels.append((ax, label_items))

            ax.grid(True, color=GRID_COLOR, linewidth=0.6, alpha=0.8)
            ax.set_axisbelow(True)
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)
            for spine in ("left", "bottom"):
                ax.spines[spine].set_color(GRID_COLOR)
            ax.tick_params(colors=TEXT_SECONDARY, labelsize=8.5)
            # Conservation errors run to 4 decimals, so the default locator can
            # pack more x ticks than the panel is wide and run the labels
            # together.  Cap the count rather than shrinking the font.
            ax.xaxis.set_major_locator(MaxNLocator(nbins=5, prune=None))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=6, prune=None))

            if i == 0:
                ax.set_title(region, fontsize=11, color=TEXT_PRIMARY, pad=8)
            if i == n_rows - 1:
                unit = VAR_UNITS.get(variable, "")
                unit = f" [{unit}]" if unit else ""
                ax.set_xlabel(
                    f"{METRIC_LABELS.get(x_metric, x_metric)}{unit}"
                    "\n$\\leftarrow$ lower is better",
                    fontsize=9, color=TEXT_SECONDARY,
                )
            if j == 0:
                unit = VAR_UNITS.get(variable, "")
                suffix = f" [{unit}]" if (unit and y_metric != "RALSD") else ""
                ax.set_ylabel(
                    f"{VAR_LABELS.get(variable, variable)}\n"
                    f"{METRIC_LABELS.get(y_metric, y_metric)}{suffix}"
                    "\n$\\leftarrow$ lower is better",
                    fontsize=9, color=TEXT_SECONDARY,
                )

    # Labels are placed only once every panel's limits are final, because the
    # collision test works in pixel space.  A draw is needed first so that the
    # renderer and the data->pixel transforms are up to date.
    fig.canvas.draw()
    for ax, items in panel_labels:
        _place_labels(fig, ax, items)

    # One legend for the whole figure: families, plus the frontier line.
    handles, legend_labels = [], []
    for family in ("afm", "naive", "pcafm"):
        if any(styles[m]["family"] == family for m in methods):
            handles.append(plt.Line2D(
                [], [], linestyle="none", marker="o", markersize=8,
                markerfacecolor=FAMILY_COLORS[family],
                markeredgecolor="#ffffff",
            ))
            legend_labels.append(FAMILY_LABELS[family])
    handles.append(plt.Line2D(
        [], [], linestyle="--", color=TEXT_SECONDARY, linewidth=1.2,
    ))
    legend_labels.append("Non-dominated frontier")

    fig.legend(
        handles, legend_labels,
        loc="lower center", ncol=min(4, len(handles)),
        frameon=False, fontsize=9, labelcolor=TEXT_PRIMARY,
        bbox_to_anchor=(0.5, -0.015),
    )

    title = cfg.get(
        "figure_title",
        "Conservation / fidelity trade-off across trained variants",
    )
    fig.suptitle(title, fontsize=12.5, color=TEXT_PRIMARY, y=0.995)
    fig.tight_layout(rect=(0, 0.045, 1, 0.975))

    out_png = Path(cfg["plot_dir"]) / f"pareto_{x_metric}_vs_{y_metric}.png"
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")

    out_pdf = out_png.with_suffix(".pdf")
    fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    logger.info("Saved Pareto figure: %s", out_png)

    frontier_df = (
        pd.concat(frontier_rows, ignore_index=True)
        if frontier_rows else points.assign(non_dominated=False)
    )
    return out_png, out_pdf, frontier_df


def save_table(frontier_df, cfg):
    """Write the plotted values plus frontier membership (the table view)."""
    y_metric = cfg["y_metric"]
    x_metric = cfg.get("x_metric", "conservation_mean")
    out_csv = (
        Path(cfg["work_dir"]) / f"pareto_{x_metric}_vs_{y_metric}_table.csv"
    )
    frontier_df.to_csv(out_csv, index=False)
    logger.info("Saved Pareto table: %s", out_csv)

    dominated = frontier_df[~frontier_df["non_dominated"]]
    if not dominated.empty:
        summary = (
            dominated.groupby("method").size().sort_values(ascending=False)
        )
        logger.info(
            "Dominated in N panels (higher = worse):\n%s", summary.to_string()
        )
    return out_csv


def write_provenance(cfg, out_files, caption):
    ancestor_files = [
        dataset["filename"] for dataset in cfg.get("input_data", {}).values()
    ]
    record = {
        "caption": caption,
        "statistics": ["other"],
        "domains": ["reg"],
        "plot_types": ["scatter"],
        "authors": ["debeire_kevin"],
        "ancestors": ancestor_files,
    }
    with ProvenanceLogger(cfg) as plog:
        for out in out_files:
            plog.log(str(out), record)


def main(cfg):
    logger.setLevel(cfg.get("log_level", "INFO").upper())

    if not cfg.get("methods"):
        raise ValueError("`methods` must be listed in the script configuration.")

    df = gather_metrics(cfg)

    # One figure per y-metric, so CRPS and RALSD each get their own panel grid
    # rather than being forced onto a shared axis.
    y_metrics = cfg.get("y_metrics") or [cfg.get("y_metric", "crps_mean")]

    outputs = []
    for y_metric in y_metrics:
        run_cfg = dict(cfg)
        run_cfg["y_metric"] = y_metric
        points, regions, variables = build_points(df, run_cfg)
        out_png, out_pdf, frontier_df = plot_pareto(
            points, regions, variables, run_cfg
        )
        out_csv = save_table(frontier_df, run_cfg)
        outputs.extend([out_png, out_pdf, out_csv])

    write_provenance(
        cfg,
        outputs,
        caption=(
            "Pareto front of conservation error against generative fidelity "
            "for all trained downscaling variants; the non-dominated set is "
            "marked."
        ),
    )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
