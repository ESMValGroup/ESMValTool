"""Python example diagnostic."""
import logging
import os
from copy import deepcopy
from pathlib import Path
from pprint import pformat

import cf_units
import iris
import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_diagnostic_filename,
    group_metadata,
    run_diagnostic,
    save_data,
    save_figure,
    select_metadata,
    sorted_metadata,
    io,
)
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(Path(__file__).stem)

LINE_LEGEND = {
    'ECS_high_hist': 'ECS_high',
    'ECS_med_hist': 'ECS_med',
    'ECS_low_hist': 'ECS_low',
}

LINE_COLOR = {
    'ECS_high_hist': 'royalblue',
    'ECS_high_scen': 'royalblue',
    'ECS_med_hist': 'green',
    'ECS_med_scen': 'green',
    'ECS_low_hist': 'orange',
    'ECS_low_scen': 'orange',
    'CMIP6': 'firebrick',
    'CMIP5': 'royalblue',
    'CMIP3': 'darkcyan',
    'OBS': 'black'
}

LINE_DASH = {
    'ECS_high_hist': 'solid',
    'ECS_high_scen': 'dashed',
    'ECS_med_hist': 'solid',
    'ECS_med_scen': 'dashed',
    'ECS_low_hist': 'solid',
    'ECS_low_scen': 'dashed',
    'CMIP6': 'solid',
    'CMIP5': 'solid',
    'CMIP3': 'solid',
    'OBS': 'solid'
}

def _get_provenance_record(caption):
    """Create a provenance record."""

    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_types': ['zonal'],
        'authors': [
            'bock_lisa',
        ],
        'references': [
            'acknow_project',
        ],
    }
    return record


def plot_icefrac(ice_frac, legend, cfg):
    """Create diagnostic data and plot it."""

    cube_label = legend
    line_color = LINE_COLOR.get(legend, legend)
    line_dash = LINE_DASH.get(legend, legend)

    plt.plot(ice_frac[0], ice_frac[1], label=cube_label, color=line_color,
              linestyle=line_dash)
    plt.xlabel('Temperature (K)')
    plt.ylabel('Ice fraction')

    logger.info("Plotting %s", legend)


def plot_errorband(x_axis, y_axis_p5, y_axis_p95, legend, cfg):
    """Create diagnostic data and plot it."""

    line_color = LINE_COLOR.get(legend, legend)
    line_dash = LINE_DASH.get(legend, legend)

    plt.fill_between(x_axis, y_axis_p5, y_axis_p95, color=line_color,
                     linestyle=line_dash, alpha=.1)

    logger.info("Plotting %s", legend)


def set_default_cfg(cfg):
    """Set default values for cfg."""
    cfg = deepcopy(cfg)
    cfg.setdefault('title_key', 'dataset')
    return cfg


def main(cfg):
    """Run diagnostic."""
    cfg = set_default_cfg(cfg)

    input_files = io.get_all_ancestor_files(cfg, pattern='ice_fraction.nc')

    for group_name in cfg['groups']:

        logger.info("Processing group %s", group_name)

        input_file = [filename for filename in input_files if group_name in filename]

        print(input_file)

        cube_ice_frac = iris.load_cube(input_file)

        print(cube_ice_frac)

        coords = cube_ice_frac.coord('dataset').points

        print(coords)

        ice_frac = []

        ice_frac.append(cube_ice_frac.data[0, :])
        ice_frac.append(cube_ice_frac.data[1, :])
        ice_frac.append(cube_ice_frac.data[2, :])
        ice_frac.append(cube_ice_frac.data[3, :])

        plot_icefrac(ice_frac, group_name, cfg)

        plot_errorband(ice_frac[0], ice_frac[2], ice_frac[3], group_name, cfg)

    caption = ("Ice fraction")

    path = get_diagnostic_filename('ice_fraction_all', cfg)

    # Provenance
    provenance_record = _get_provenance_record(caption)
    provenance_record['ancestors'] = cfg['input_files']

    basename = 'ice_fraction_all'

    # Save the data used for the plot
    save_data(basename, provenance_record, cfg, cube_ice_frac)

    title = 'Ice fraction'

    plt.title(title)
    plt.legend(ncol=1)
    plt.axhline(y=0.5, xmin=230., xmax=300., color='black', linewidth=3)
    plt.grid(True)

    # And save the plot
    save_figure(basename, provenance_record, cfg)



if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
