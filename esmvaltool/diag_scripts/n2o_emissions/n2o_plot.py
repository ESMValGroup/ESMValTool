"""N2O emissions diagnostic."""
import logging
import matplotlib
from pathlib import Path

import iris

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
)
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(Path(__file__).stem)


COLOR_MAPS = {
    "wbor_r": matplotlib.colors.LinearSegmentedColormap(
        "wbor_r",
        {
            "red": (
                (0.0, 0.5, 0.5),
                (0.125, 1.0, 1.0),
                (0.375, 0.925, 0.925),
                (0.625, 0.28, 0.28),
                (1.0, 0.968, 0.968),
            ),
            "green": (
                (0.0, 0.0, 0.0),
                (0.125, 0.0, 0.0),
                (0.375, 0.7127, 0.7127),
                (0.625, 0.6383, 0.6383),
                (0.875, 0.7978, 0.79787),
                (1.0, 0.957, 0.957),
            ),
            "blue": (
                (0.0, 0.0, 0.0),
                (0.375, 0.0, 0.0),
                (0.625, 0.77, 0.77),
                (0.875, 0.86, 0.86),
                (1.0, 0.97, 0.97),
            ),
        }
    )
}


def get_provenance_record(attributes, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    caption = attributes['caption'].format(**attributes)

    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_types': ['zonal'],
        'authors': [
            'andela_bouwe',
            'righi_mattia',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': ancestor_files,
    }
    return record


def prep_cube(filename):
    """Load and squeeze cube."""
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)
    cube = iris.util.squeeze(cube)
    return cube


def prep_config(cfg):
    if cfg.get('quickplot'):
        quickplot_cfg = cfg.get('quickplot')

        cmap = quickplot_cfg.get("cmap")
        if cmap in COLOR_MAPS:
            quickplot_cfg["cmap"] = COLOR_MAPS[cmap]

        vmin = quickplot_cfg.pop("vmin", None)
        vmax = quickplot_cfg.pop("vmax", None)
        if vmin is not None or vmax is not None:
            quickplot_cfg["norm"] = matplotlib.colors.Normalize(
                vmin, vmax, clip=True
            )


def plot_diagnostic(cube, basename, provenance_record, cfg):
    """Create diagnostic data and plot it."""

    # Save the data used for the plot
    save_data(basename, provenance_record, cfg, cube)

    if cfg.get('quickplot'):
        # Create the plot
        quickplot(cube, **cfg['quickplot'])
        # And save the plot
        save_figure(basename, provenance_record, cfg)


def main(cfg):
    """Compute the time average for each input dataset."""

    prep_config(cfg)

    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    # Loop over variables/datasets in alphabetical order
    groups = group_metadata(input_data, 'variable_group', sort='activity')
    for group_name in groups:
        logger.info("Processing variable %s", group_name)
        for attributes in groups[group_name]:
            logger.info("Processing dataset %s", attributes['dataset'])
            input_file = attributes['filename']
            cube = prep_cube(input_file)

            output_basename = Path(input_file).stem
            if group_name != attributes['short_name']:
                output_basename = group_name + '_' + output_basename

            if "caption" not in attributes:
                attributes['caption'] = input_file

            provenance_record = get_provenance_record(
                attributes, ancestor_files=[input_file]
            )

            plot_diagnostic(cube, output_basename, provenance_record, cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
