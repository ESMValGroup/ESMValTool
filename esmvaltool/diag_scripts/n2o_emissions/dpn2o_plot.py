"""N2O emissions diagnostic."""
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
    select_metadata,
)
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(Path(__file__).stem)


COLOR_MAPS = {
    'wbor': matplotlib.colors.LinearSegmentedColormap(
        'wbor',
        {
            'red': [
                (0.0, 0.968, 0.968),
                (0.375, 0.28, 0.28),
                (0.625, 0.925, 0.925),
                (0.875, 1.0, 1.0),
                (1.0, 0.5, 0.5),
            ],
            'green': [
                (0.0, 0.957, 0.957),
                (0.125, 0.7978, 0.79787),
                (0.375, 0.6383, 0.6383),
                (0.625, 0.7127, 0.7127),
                (0.875, 0.0, 0.0),
                (1.0, 0.0, 0.0),
            ],
            'blue': [
                (0, 0.97, 0.97),
                (0.125, 0.86, 0.86),
                (0.375, 0.77, 0.77),
                (0.625, 0.0, 0.0),
                (1.0, 0.0, 0.0),
            ],
        }
    ),
    'GMT_polar': matplotlib.colors.LinearSegmentedColormap(
        'GMT_polar',
        {
            'red': [(0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 1.0, 1.0)],
            'green': [(0.0, 0.0, 0.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0)],
            'blue': [(0.0, 1.0, 1.0), (0.5, 1.0, 1.0), (1.0, 0.0, 0.0)],
        },
        matplotlib.rcParams['image.lut'],
    )
}


def get_provenance_record(caption, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_types': ['zonal'],
        'authors': ['righi_mattia'],
        'references': ['acknow_project'],
        'ancestors': ancestor_files,
    }
    return record


def prep_config(cfg):
    """Additional plot configuration."""
    for plot_type in ['data_plot', 'diff_plot']:
        plot_cfg = cfg.get(plot_type)

        if plot_cfg is not None:
            plot_cfg = cfg.get(plot_type)

            cmap = plot_cfg.get('cmap')
            if cmap in COLOR_MAPS:
                plot_cfg['cmap'] = COLOR_MAPS[cmap]

            vmin = plot_cfg.get('vmin')
            vmax = plot_cfg.get('vmax')
            levs = plot_cfg.get('levs')
            if vmin is not None and vmax is not None and levs is not None:
                plot_cfg['levels'] = np.linspace(vmin, vmax, levs)
                plot_cfg['extend'] = 'both'


def prep_cube(filename):
    """Load and squeeze cube."""
    logger.debug('Loading %s', filename)
    cube = iris.load_cube(filename)
    cube = iris.util.squeeze(cube)
    return cube


def _plot(cube, basename, provenance_record, plot_type, cfg):
    """Save diagnostic data and plot it."""
    save_data(basename, provenance_record, cfg, cube)
    if cfg.get(plot_type):
        quickplot(cube, **cfg[plot_type])
        save_figure(basename, provenance_record, cfg)


def plot_diff(cube, attributes, obs_cube, obs_attributes, cfg):
    input_file = attributes['filename']
    obs_input_file = obs_attributes['filename']
    input_file_stem = Path(input_file).stem
    obs_input_file_stem = Path(obs_input_file).stem

    # Produce caption and provenance.
    caption = f'Difference between {input_file_stem} and {obs_input_file_stem}'
    provenance_record = get_provenance_record(
        caption, ancestor_files=[input_file, obs_input_file]
    )

    # Compute difference and plot.
    cube_diff = cube - obs_cube
    output_basename = f'diff_{input_file_stem}_{obs_input_file_stem}'
    _plot(cube_diff, output_basename, provenance_record, 'diff_plot', cfg)


def plot_dataset(cube, attributes, cfg):
    input_file = attributes['filename']
    input_file_stem = Path(input_file).stem

    # Produce caption and provenance.
    caption = input_file_stem
    provenance_record = get_provenance_record(
        caption, ancestor_files=[input_file]
    )

    output_basename = input_file_stem
    _plot(cube, output_basename, provenance_record, 'data_plot', cfg)


def main(cfg):
    """Compute the time average for each input dataset."""

    prep_config(cfg)

    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    # Get Yang2020 climatology from among the datasets.
    obs_clim_list = select_metadata(input_data, dataset='Yang2020')
    attributes_yang2020 = obs_clim_list[0] if len(obs_clim_list) == 1 else None

    cube_yang2020 = prep_cube(attributes_yang2020['filename'])
    plot_dataset(cube_yang2020, attributes_yang2020, cfg)

    # Loop over datasets.
    grouped_datasets = group_metadata(input_data, 'dataset', sort='activity')
    for (dataset, attributes_by_activity) in grouped_datasets.items():
        if dataset != attributes_yang2020['dataset']:
            for attributes in attributes_by_activity:

                logger.info(
                    f'Processing dataset {dataset} ({attributes["activity"]})'
                )
                cube = prep_cube(attributes['filename'])
                plot_dataset(cube, attributes, cfg)
                plot_diff(
                    cube, attributes, cube_yang2020, attributes_yang2020, cfg
                )


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
