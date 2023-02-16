"""Plot time series."""
import logging
from copy import deepcopy
from pathlib import Path

import iris
import matplotlib.pyplot as plt

from esmvaltool.diag_scripts.mlr.plot import unify_time_coord
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(Path(__file__).stem)


def get_plot_kwargs(cfg, dataset, **kwargs):
    """Get plot kwargs for different datasets."""
    plot_kwargs = dict(kwargs)
    plot_kwargs_for_keys = {
        'ERA5': dict(color='black'),
        'ICON': dict(color='C1'),
        # TODO REMOVE
        # 'IRIS REMAPNN': dict(color='black'),
        # 'CDO REMAPNN': dict(color='C0'),
        # 'CDO REMAPCON': dict(color='C1'),
        # 'CDO REMAPDIS': dict(color='C2'),
    }
    key = dataset.get(cfg['title_key'])
    if key is None:
        return plot_kwargs
    plot_kwargs.update(plot_kwargs_for_keys.get(key, {}))
    return plot_kwargs


def load_and_preprocess(dataset):
    """Load and preprocess data."""
    filename = dataset['filename']
    logger.info("Loading %s", filename)
    cube = iris.load_cube(filename)
    unify_time_coord(cube)

    # Global mean
    if all([cube.coords('latitude', dim_coords=True),
            cube.coords('longitude', dim_coords=True)]):
        area_weights = iris.analysis.cartography.area_weights(cube)
        cube = cube.collapsed(['latitude', 'longitude'], iris.analysis.MEAN,
                              weights=area_weights)

    # Convert units of some variables
    if cube.var_name == 'tas':
        cube.convert_units('celsius')
    if cube.var_name == 'pr':
        cube.units = 'mm s-1'
        cube.convert_units('mm day-1')
    if cube.var_name == 'hus':
        cube.convert_units('g kg-1')
    return cube


def main(cfg):
    """Run diagnostic."""
    cfg = deepcopy(cfg)
    cfg.setdefault('title_key', 'dataset')
    logger.info("Using key '%s' to create titles for datasets",
                cfg['title_key'])

    # This diagnostic supports only a single variable, raise error otherwise
    input_data = list(cfg['input_data'].values())
    all_vars = list(group_metadata(input_data, 'short_name'))
    if len(all_vars) != 1:
        raise ValueError(
            f"Expected exactly 1 variable, got {len(all_vars):d}: {all_vars}")

    # Set fixed size for figure
    plt.figure(figsize=(8, 4), dpi=150)

    # Plot all datasets in single figure
    ancestors = []
    for dataset in input_data:
        ancestors.append(dataset['filename'])
        cube = load_and_preprocess(dataset)

        # Plot monthly means
        plot_kwargs = get_plot_kwargs(cfg, dataset)
        iris.plot.plot(cube, label=dataset[cfg['title_key']], **plot_kwargs)

        # Plot annual means
        iris.coord_categorisation.add_year(cube, 'time')
        cube = cube.aggregated_by('year', iris.analysis.MEAN)
        plot_kwargs['linestyle'] = ':'
        iris.plot.plot(cube, **plot_kwargs)

    # Plot appearance
    long_name = input_data[0]['long_name']
    short_name = input_data[0]['short_name']
    units = cube.units
    plt.title(f"Global Mean {long_name}")
    plt.xlabel("Year")
    plt.ylabel(f"{short_name} [{units}]")
    plt.legend()

    # Save plot
    plot_path = get_plot_filename(short_name, cfg)
    plt.savefig(plot_path, bbox_inches='tight', orientation='landscape')
    logger.info("Wrote %s", plot_path)
    plt.close()

    # Provenance tracking
    caption = (f"Monthly mean (solid lines) and annual mean (dashed lines) "
               f"time series of {input_data[0]['long_name']} for various "
               f"datasets.")
    provenance_record = {
        'ancestors': ancestors,
        'authors': ['schlund_manuel'],
        'caption': caption,
        'plot_types': ['line'],
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(plot_path, provenance_record)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
