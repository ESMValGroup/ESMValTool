#!/usr/bin/env python
"""Python example diagnostic."""
import logging
import os
from pprint import pformat

import iris
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

from esmvalcore.preprocessor import multi_model_statistics
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvaltool.diag_scripts.shared._base import get_plot_filename
from mpqb_plots import read_mpqb_cfg

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    """Create lineplot."""
    ylims = [cfg.pop('y0', None), cfg.pop('y1', None)]

    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    mpqb_cfg = read_mpqb_cfg()
    datasetnames = mpqb_cfg['datasetnames']

    grouped_input_data = group_metadata(input_data, 'alias', sort='alias')

    logger.info(
        "Example of how to group and sort input data by standard_name:"
        "\n%s", pformat(grouped_input_data))

    # In order to get the right line colors for MPQB soil moisture
    # here we put ERA-Interim-Land at the end of the dictionary if
    # it is included.
    if 'ERA-Interim-Land' in grouped_input_data.keys():
        grouped_input_data.move_to_end('ERA-Interim-Land')

    # Get the multimodel statistics from the preprocessor
    statistics_data = multi_model_statistics([cube1,...,cubeN], 'overlap', ['mean'])
    #statistics_data = multi_model_statistics
    mean_cube = statistics_data['mean']
        
    plt.clf()
    fig = plt.figure(figsize=(10, 4))
    ax1 = fig.add_subplot()

    for dataset in grouped_input_data:
        dataset_cfg = grouped_input_data[dataset][0]
        alias = dataset_cfg['alias']

        logger.info("Opening dataset: %s", dataset)
        cube = iris.load_cube(dataset_cfg['filename'])

        anom_cube = mean_cube - cube
        
        iris.quickplot.plot(anom_cube, label=datasetnames[alias],
                            color=mpqb_cfg['datasetcolors'][alias])
    plt.legend()
    plt.xticks(rotation=90)
    # Add the zero line when plotting anomalies
    #if 'ano' in dataset_cfg['preprocessor']:
    plt.axhline(y=0, linestyle=':', color='k')
    plt.tight_layout()
    # Time axis formatting
    years = mdates.YearLocator()  # every year
    years_fmt = mdates.DateFormatter('%Y')
    ax1 = plt.gca()
    ax1.xaxis.set_major_locator(years)
    ax1.xaxis.set_major_formatter(years_fmt)
    ax1.grid(True, which='major', axis='x')

    ax1.set_ylim(ylims)

    baseplotname = f"lineplot_{dataset_cfg['variable_group']}_{dataset_cfg['start_year']}-{dataset_cfg['end_year']}"

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
