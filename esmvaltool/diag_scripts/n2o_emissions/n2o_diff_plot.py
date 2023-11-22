"""N2O emissions plot diagnostic."""
import iris
import logging
import matplotlib
import numpy as np
from pathlib import Path

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
)
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(Path(__file__).stem)


COLOR_MAPS = {
    "wbor": matplotlib.colors.LinearSegmentedColormap(
        "wbor",
        {
            "red": [
                (0.0, 0.968, 0.968),
                (0.375, 0.28, 0.28),
                (0.625, 0.925, 0.925),
                (0.875, 1.0, 1.0),
                (1.0, 0.5, 0.5),
            ],
            "green": [
                (0.0, 0.957, 0.957),
                (0.125, 0.7978, 0.79787),
                (0.375, 0.6383, 0.6383),
                (0.625, 0.7127, 0.7127),
                (0.875, 0.0, 0.0),
                (1.0, 0.0, 0.0),
            ],
            "blue": [
                (0, 0.97, 0.97),
                (0.125, 0.86, 0.86),
                (0.375, 0.77, 0.77),
                (0.625, 0.0, 0.0),
                (1.0, 0.0, 0.0),
            ],
        },
    ),
    "GMT_polar": matplotlib.colors.LinearSegmentedColormap(
        "GMT_polar",
        {
            "red": [(0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)],
            "green": [(0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0)],
            "blue": [(0.0, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0)],
        },
        matplotlib.rcParams["image.lut"],
    ),
}


def get_provenance_record(caption, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        "caption": caption,
        "ancestors": ancestor_files,
    }
    return record


def prep_config(cfg):
    """Prepare plot configuration."""
    for plot_type in ["data_plot", "diff_plot"]:
        plot_cfg = cfg.get(plot_type)

        if plot_cfg is not None:
            # Check if color map is set to one of the custom ones
            # defined by this diagnostic.
            cmap = plot_cfg.get("cmap")
            if cmap in COLOR_MAPS:
                plot_cfg["cmap"] = COLOR_MAPS[cmap]

            # Get plot limits and levels if they are defined
            # in the recipe.
            vmin = plot_cfg.get("vmin")
            vmax = plot_cfg.get("vmax")
            levs = plot_cfg.get("levs")
            if None not in [vmin, vmax, levs]:
                plot_cfg["levels"] = np.linspace(vmin, vmax, levs)
                plot_cfg["extend"] = "both"


def save_data_plot(cube, basename, provenance_record, plot_type, cfg):
    """Save diagnostic data and plot figure."""
    save_data(basename, provenance_record, cfg, cube)
    if cfg.get(plot_type):
        quickplot(cube, **cfg[plot_type])
        save_figure(basename, provenance_record, cfg)


def plot_diff(cube, input_file, obs_cube, obs_input_file, cfg):
    input_file_stem = Path(input_file).stem
    obs_input_file_stem = Path(obs_input_file).stem

    # Produce caption and provenance.
    caption = f"Difference between {input_file_stem} and {obs_input_file_stem}"
    provenance_record = get_provenance_record(
        caption, ancestor_files=[input_file, obs_input_file]
    )

    # Compute difference and plot.
    cube_diff = cube - obs_cube
    output_basename = f"diff_{input_file_stem}_{obs_input_file_stem}"
    save_data_plot(cube_diff, output_basename, provenance_record, "diff_plot", cfg)


def plot_dataset(cube, input_file, cfg):
    input_file_stem = Path(input_file).stem

    # Produce caption and provenance.
    caption = input_file_stem
    provenance_record = get_provenance_record(
        caption,
        ancestor_files=[input_file],
    )

    output_basename = input_file_stem
    save_data_plot(cube, output_basename, provenance_record, "data_plot", cfg)


def main(cfg):
    """Compute the time average for each input dataset."""

    # Modify configuration if necessary.
    prep_config(cfg)

    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg["input_data"].values()
    datasets = group_metadata(input_data, "dataset")

    # Extract reference dataset
    ref_dataset = datasets.pop(cfg["reference_dataset"])
    ref_filename = ref_dataset[0]["filename"]
    ref_cube = iris.load_cube(ref_filename)
    plot_dataset(ref_cube, ref_filename, cfg)

    # Loop over datasets.
    for model_dataset, group in datasets.items():
        logger.info(f"Processing dataset {model_dataset}")
        logger.info(group)

        for attributes in group:
            filename = attributes["filename"]
            cube = iris.load_cube(filename)
            plot_dataset(cube, filename, cfg)
            plot_diff(cube, filename, ref_cube, ref_filename, cfg)


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
