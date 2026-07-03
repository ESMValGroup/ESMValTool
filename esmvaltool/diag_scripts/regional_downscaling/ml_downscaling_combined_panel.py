"""ML-based Downscaling Combined Panel Diagnostic.

This diagnostic gathers plots from ancestor diagnostics (spatial_metrics,
energy_spectrum, log_density, calibration) and combines them into a single
publication-quality panel figure per variable, with subplot labels (A, B, C, D).

All combined panels have the same shape: narrower spatial metrics images
(variables with fewer metric columns, e.g. tas with only bias+CRPS) are
left-aligned and padded with white on the right to match the widest variable
(typically pr, which has bias+relbias+CRPS+conservation).

Layout:
  A (top, full width)    : Spatial metrics maps (left-aligned, white-padded)
  B (bottom-left)        : Energy spectrum
  C (bottom-center)      : Log PDF
  D (bottom-right)       : Rank histograms
"""

import logging
import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec

from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared._base import ProvenanceLogger

logger = logging.getLogger(os.path.basename(__file__))

# Expected plot filename patterns per analysis type
PLOT_PATTERNS = {
    "spatial_metrics": "spatial_metrics_{var}.png",
    "energy_spectrum": "energy_spectrum_{var}.png",
    "log_density": "log_pdf_{var}.png",
    "calibration": "rank_histogram_{var}.png",
}

# Panel labels and their analysis types (in order)
PANELS = [
    ("A", "spatial_metrics", "Spatial Metrics"),
    ("B", "energy_spectrum", "Energy Spectrum"),
    ("C", "log_density", "Log PDF"),
    ("D", "calibration", "Rank Histogram"),
]

VARIABLE_LABELS = {
    "pr": "Precip. (pr)",
    "huss": "Spec. Humid. (huss)",
    "tas": "Temp. (tas)",
    "ps": "Sfc. Press. (ps)",
    "uas": "U-Wind (uas)",
    "vas": "V-Wind (vas)",
}

def _get_provenance_record(cfg, plot_file, caption):
    """Create a provenance record for the combined panel plot."""
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


def find_plot_files(cfg, var_name):
    """Search ancestor directories for plot images for a given variable.

    Parameters
    ----------
    cfg : dict
        ESMValTool configuration dictionary
    var_name : str
        Variable name (e.g., 'tas', 'pr')

    Returns
    -------
    dict
        {analysis_type: path_string} for each found plot
    """
    ancestor_dirs = cfg.get("input_files", [])
    found = {}

    # Build search targets
    targets = {}
    for analysis_type, pattern in PLOT_PATTERNS.items():
        targets[analysis_type] = pattern.format(var=var_name)

    # Search ancestor directories (plots are typically in plots/ subdirectory)
    for ancestor_dir in ancestor_dirs:
        ancestor_path = ancestor_dir.replace("work", "plots")
        ancestor_path = Path(ancestor_path)
        # ESMValTool ancestors point to work/ dirs; plots/ is a sibling
        search_dirs = [
            ancestor_path,
            ancestor_path.parent,
            ancestor_path.parent / "plots",
        ]
        # Also search recursively
        for search_dir in search_dirs:
            if not search_dir.exists():
                continue
            for analysis_type, filename in targets.items():
                if analysis_type in found:
                    continue
                for match in search_dir.rglob(filename):
                    logger.info("Found %s plot: %s", analysis_type, match)
                    found[analysis_type] = str(match)
                    break

    # Fallback: search entire run directory tree
    if len(found) < len(targets):
        run_dir = Path(cfg["work_dir"]).parent.parent
        logger.info("Searching run directory for missing plots: %s", run_dir)
        for analysis_type, filename in targets.items():
            if analysis_type in found:
                continue
            for match in run_dir.rglob(filename):
                logger.info("Found %s plot (fallback): %s", analysis_type, match)
                found[analysis_type] = str(match)
                break

    return found


def pad_image_to_aspect(img, target_ar):
    """Pad an image with white on the right to match a target aspect ratio.

    The image is left-aligned; white padding is added on the right side only.
    The height stays the same.

    Parameters
    ----------
    img : numpy.ndarray
        Image array (H, W, C) with values in [0, 1] (float) or [0, 255] (uint8)
    target_ar : float
        Target aspect ratio (width / height). Must be >= current AR.

    Returns
    -------
    numpy.ndarray
        Padded image with the target aspect ratio
    """
    h, w = img.shape[0], img.shape[1]
    current_ar = w / h

    if current_ar >= target_ar - 1e-3:
        # Already wide enough (or wider), no padding needed
        return img

    # Compute new width to match target AR
    new_w = int(round(h * target_ar))

    # Create white canvas
    if img.dtype == np.uint8:
        white_val = 255
    else:
        white_val = 1.0

    if img.ndim == 3:
        n_channels = img.shape[2]
        padded = np.full((h, new_w, n_channels), white_val, dtype=img.dtype)
        padded[:, :w, :] = img
    else:
        padded = np.full((h, new_w), white_val, dtype=img.dtype)
        padded[:, :w] = img

    logger.info("Padded image from %dx%d (AR=%.3f) to %dx%d (AR=%.3f)",
                w, h, current_ar, new_w, h, target_ar)

    return padded


def create_combined_panel(plot_files, var_name, cfg, reference_ar=None):
    """Create a combined panel figure from individual plot images.

    Layout:
        Row 0 (tall): A — spatial metrics (left-aligned, padded to reference_ar)
        Row 1:        B — energy spectrum | C — log PDF | D — rank histogram

    Parameters
    ----------
    plot_files : dict
        {analysis_type: filepath} for each available plot
    var_name : str
        Variable name
    cfg : dict
        Configuration dictionary
    reference_ar : float, optional
        Reference aspect ratio for the spatial metrics panel. If provided,
        narrower images are padded with white on the right to match this AR.
        This ensures all combined panels have the same shape.
    """
    # Load available images
    images = {}
    for label, analysis_type, title in PANELS:
        if analysis_type in plot_files:
            img = mpimg.imread(plot_files[analysis_type])
            images[analysis_type] = img
            logger.info("Loaded %s: %s  (shape %s)",
                        label, analysis_type, img.shape)
        else:
            logger.warning("Missing plot for panel %s (%s) — %s",
                           label, analysis_type, var_name)

    if not images:
        logger.error("No plot images found for variable %s, skipping.", var_name)
        return

    # =========================================================================
    # Pad spatial metrics to reference aspect ratio if needed
    # =========================================================================
    has_A = "spatial_metrics" in images
    if has_A and reference_ar is not None:
        images["spatial_metrics"] = pad_image_to_aspect(
            images["spatial_metrics"], reference_ar
        )

    # =========================================================================
    # Determine layout proportions
    # =========================================================================
    bottom_panels = [
        (label, atype, title)
        for label, atype, title in PANELS[1:]
        if atype in images
    ]
    n_bottom = len(bottom_panels)

    if n_bottom == 0 and not has_A:
        logger.error("No panels available for %s", var_name)
        return

    # =========================================================================
    # Build figure with GridSpec
    # =========================================================================
    fig_width = cfg.get("panel_fig_width", 18)

    # Panel A aspect ratio (after padding)
    if has_A:
        img_a = images["spatial_metrics"]
        ar_a = img_a.shape[1] / img_a.shape[0]
        h_a = fig_width / ar_a
    else:
        h_a = 0

    # Bottom row: width_ratios from image aspect ratios
    if n_bottom > 0:
        bottom_aspects = []
        for _, atype, _ in bottom_panels:
            img = images[atype]
            bottom_aspects.append(img.shape[1] / img.shape[0])

        # Optionally widen the rank-histogram panel (D) so the per-model
        # histograms are larger and more legible (R2-Fig4-7-10). The aspect
        # used to size the row height (h_b) is kept at the true image aspect,
        # so panel D itself is not distorted; only its share of the row width
        # grows. Default 1.0 leaves all other recipes unchanged.
        rank_width_scale = cfg.get("rank_width_scale", 1.0)
        width_ratios = list(bottom_aspects)
        for i, (_, atype, _) in enumerate(bottom_panels):
            if atype == "calibration":
                width_ratios[i] = bottom_aspects[i] * rank_width_scale

        total_ratio = sum(width_ratios)
        h_b = max(
            (fig_width * (wr / total_ratio)) / ar
            for wr, ar in zip(width_ratios, bottom_aspects)
        )
    else:
        h_b = 0
        width_ratios = [1]

    # Total figure height
    gap_inches = 0.6
    fig_height = h_a + h_b + gap_inches + 1.5

    # Height ratios
    if has_A and n_bottom > 0:
        height_ratios = [h_a, h_b]
        n_rows = 2
    elif has_A:
        height_ratios = [1]
        n_rows = 1
    else:
        height_ratios = [1]
        n_rows = 1

    n_cols = max(n_bottom, 1)

    fig = plt.figure(figsize=(fig_width, fig_height))

    gs = GridSpec(
        n_rows, n_cols, figure=fig,
        height_ratios=height_ratios if n_rows > 1 else [1],
        width_ratios=width_ratios[:n_cols],
        hspace=0.08,
        wspace=0.08,
        left=0.02, right=0.98,
        top=0.96, bottom=0.02,
    )

    # =========================================================================
    # Panel A — Spatial metrics (full width, top row)
    # =========================================================================
    if has_A:
        if n_cols > 1:
            ax_a = fig.add_subplot(gs[0, :])
        else:
            ax_a = fig.add_subplot(gs[0, 0])
        ax_a.imshow(images["spatial_metrics"], aspect="auto")
        ax_a.axis("off")

        # Add label "A"
        ax_a.text(
            0.01, 0.98, "A",
            transform=ax_a.transAxes,
            fontsize=22, fontweight="bold",
            va="top", ha="left",
            bbox=dict(
                facecolor="white", edgecolor="black",
                boxstyle="round,pad=0.3", alpha=0.9, linewidth=1.5,
            ),
            zorder=10,
        )

    # =========================================================================
    # Panels B, C, D — Bottom row
    # =========================================================================
    if n_bottom > 0:
        row_idx = 1 if has_A else 0
        label_start = ord("B") if has_A else ord("A")
        for col, (_, atype, title) in enumerate(bottom_panels):
            ax = fig.add_subplot(gs[row_idx, col])
            ax.imshow(images[atype], aspect="auto")
            ax.axis("off")

            label = chr(label_start + col)
            ax.text(
                0.01, 0.98, label,
                transform=ax.transAxes,
                fontsize=22, fontweight="bold",
                va="top", ha="left",
                bbox=dict(
                    facecolor="white", edgecolor="black",
                    boxstyle="round,pad=0.3", alpha=0.9, linewidth=1.5,
                ),
                zorder=10,
            )
    # =========================================================================
    # Supertitle
    # =========================================================================
    region_name = cfg.get("region_name", "")
    variable_name = VARIABLE_LABELS.get(var_name, var_name)
    suptitle = f"{variable_name}, {region_name}" if region_name else variable_name
    fig.suptitle(suptitle, fontsize=20, fontweight="bold", y=0.99)
    # =========================================================================
    # Save
    # =========================================================================
    plot_file = os.path.join(
        cfg["plot_dir"],
        f"combined_panel_{var_name}.png",
    )
    plt.savefig(plot_file, dpi=300, bbox_inches="tight", facecolor="white")

    plot_file_pdf = plot_file.replace(".png", ".pdf")
    plt.savefig(plot_file_pdf, bbox_inches="tight", facecolor="white")
    plt.close()

    caption = (
        f"Combined evaluation panel for {var_name}: "
        f"(A) spatial metrics, (B) energy spectrum, "
        f"(C) log PDF, (D) rank histograms."
    )
    _get_provenance_record(cfg, plot_file, caption)

    logger.info("Saved combined panel: %s", plot_file)
    logger.info("Saved combined panel PDF: %s", plot_file_pdf)

    return plot_file


def main(cfg):
    """Run the combined panel diagnostic."""
    logger.setLevel(cfg.get("log_level", "INFO").upper())

    variables = cfg.get("variables", ["pr", "huss", "tas", "ps", "uas", "vas"])
    logger.info("Creating combined panels for variables: %s", variables)

    # =========================================================================
    # First pass: find all spatial metrics images and determine the widest one
    # (highest aspect ratio). This becomes the reference — all narrower images
    # will be padded with white on the right to match.
    # =========================================================================
    all_plot_files = {}
    reference_ar = 0.0

    for var_name in variables:
        plot_files = find_plot_files(cfg, var_name)
        all_plot_files[var_name] = plot_files

        if "spatial_metrics" in plot_files:
            img = mpimg.imread(plot_files["spatial_metrics"])
            ar = img.shape[1] / img.shape[0]  # width / height
            logger.info("Spatial metrics %s: %dx%d pixels, AR=%.3f",
                        var_name, img.shape[1], img.shape[0], ar)
            if ar > reference_ar:
                reference_ar = ar

    if reference_ar > 0:
        logger.info("Reference spatial metrics AR: %.3f (widest variable)",
                     reference_ar)
    else:
        reference_ar = None
        logger.info("No spatial metrics images found, skipping AR normalization")

    # =========================================================================
    # Second pass: create combined panels with consistent shape
    # =========================================================================
    for var_name in variables:
        plot_files = all_plot_files[var_name]

        if not plot_files:
            logger.warning("No plots found for %s, skipping.", var_name)
            continue

        logger.info("Found %d/%d plots for %s: %s",
                     len(plot_files), len(PLOT_PATTERNS), var_name,
                     list(plot_files.keys()))

        create_combined_panel(plot_files, var_name, cfg,
                              reference_ar=reference_ar)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)