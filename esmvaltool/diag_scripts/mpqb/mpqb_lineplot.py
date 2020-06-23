#!/usr/bin/env python
"""Python example diagnostic."""
import logging
import os
from pprint import pformat

import iris
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvaltool.diag_scripts.shared._base import get_plot_filename


DATASET_PLOTNAMES = {
    'ERA-Interim-Land': 'ERA-Interim/Land',
    'CDS-SATELLITE-SOIL-MOISTURE': 'ESA-CCI',
    'cds-era5-land-monthly': 'ERA5-Land',
    'cds-era5-monthly': 'ERA5',
    'MERRA2': 'MERRA-2',
    'cds-satellite-lai-fapar': 'SPOT-VGT',
    'CDS-SATELLITE-ALBEDO': 'SPOT-VGT',
}

YLIMS = {'sm': (0.22, 0.28), 'sm1m': (0.22, 0.28)}

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    """Create lineplot."""
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    grouped_input_data = group_metadata(input_data, 'dataset', sort='dataset')

    logger.info(
        "Example of how to group and sort input data by standard_name:"
        "\n%s", pformat(grouped_input_data))

    # In order to get the right line colors for MPQB soil moisture
    # here we put ERA-Interim-Land at the end of the dictionary if
    # it is included.
    if 'ERA-Interim-Land' in grouped_input_data.keys():
        grouped_input_data.move_to_end('ERA-Interim-Land')

    plt.clf()
    fig = plt.figure(figsize=(10, 4))
    ax1 = fig.add_subplot()
    for dataset in grouped_input_data:
        dataset_cfg = grouped_input_data[dataset][0]
        logger.info("Opening dataset: %s", dataset)
        cube = iris.load_cube(dataset_cfg['filename'])
        iris.quickplot.plot(cube, label=DATASET_PLOTNAMES[dataset])
    plt.legend()
    plt.xticks(rotation=90)
    # Add the zero line when plotting anomalies
    if 'ano' in dataset_cfg['preprocessor']:
        plt.axhline(y=0, linestyle=':', color='k')
    plt.tight_layout()
    # Time axis formatting
    years = mdates.YearLocator()  # every year
    years_fmt = mdates.DateFormatter('%Y')
    ax1 = plt.gca()
    ax1.xaxis.set_major_locator(years)
    ax1.xaxis.set_major_formatter(years_fmt)
    ax1.grid(True, which='major', axis='x')

    # set ylims if specified
    if dataset_cfg['variable_group'] in YLIMS:
        ax1.set_ylim(YLIMS[dataset_cfg['variable_group']])

    baseplotname = f"lineplot_{dataset_cfg['variable_group']}\
                    _{dataset_cfg['start_year']}\
                    -{dataset_cfg['end_year']}"

    filename = get_plot_filename(baseplotname, cfg)
    logger.info("Saving as %s", filename)
    fig.savefig(filename)
    plt.close(fig)
    logger.info("Finished!")


if __name__ == '__main__':
    with run_diagnostic() as config:
        if config['write_plots']:
            main(config)
        else:
            logger.warning("This diagnostic wants to plot,\
                            but isn't allowed to")
