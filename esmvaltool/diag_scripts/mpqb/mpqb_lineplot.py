#!/usr/bin/env python
"""Python example diagnostic."""
import logging
import os
from pprint import pformat
from esmvaltool.diag_scripts.shared.trend_mpqb_common.sharedutils import parallel_apply_along_axis
from esmvaltool.diag_scripts.shared.trend_mpqb_common.diag1d import *
import numpy as np

import matplotlib.pyplot as plt
import iris
import itertools as it
import warnings

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata, sorted_metadata)
from esmvaltool.diag_scripts.shared._base import (
    ProvenanceLogger, get_diagnostic_filename, get_plot_filename)
from esmvaltool.diag_scripts.shared.plot import quickplot


logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    grouped_input_data = group_metadata(
        input_data, 'dataset', sort='dataset')

    logger.info(
        "Example of how to group and sort input data by standard_name:"
        "\n%s", pformat(grouped_input_data))

    # In order to get the right line colors for MPQB soil moisture
    # here we put ERA-Interim-Land at the end of the dictionary if
    # it is included. 
    if 'ERA-Interim-Land' in grouped_input_data.keys():
        grouped_input_data.move_to_end('ERA-Interim-Land')

    if cfg['write_plots']:
        plt.clf()
        fig = plt.figure(figsize=(10,4))
        ax = fig.add_subplot()
        for dataset in grouped_input_data:
            logger.info("Opening dataset: {0}".format(dataset))
            cube = iris.load_cube(grouped_input_data[dataset][0]['filename'])
            iris.quickplot.plot(cube, label=grouped_input_data[dataset][0].pop(
                               'alias',
                               grouped_input_data[dataset][0].pop('dataset')))
        plt.legend()
        plt.xticks(rotation=90)
        # Add the zero line when plotting anomalies
        if 'ano' in grouped_input_data[dataset][0]['preprocessor']:
            plt.axhline(y=0, linestyle=':', color='k')
        plt.tight_layout()
        # Time axis formatting
        import matplotlib.dates as mdates
        years = mdates.YearLocator()   # every year
        years_fmt = mdates.DateFormatter('%Y')
        ax = plt.gca()
        ax.xaxis.set_major_locator(years)
        ax.xaxis.set_major_formatter(years_fmt)
        ax.grid(True, which='major', axis='x')
        filename = get_plot_filename('lineplot', cfg)
        logger.info("Saving as %s", filename)
        fig.savefig(filename)
        plt.close(fig)
    else:
        logger.warning("This diagnostic wants to plot, but isn't allowed to")
    logger.info("Finished!")


if __name__ == '__main__':
    with run_diagnostic() as cfg:
        main(cfg)
