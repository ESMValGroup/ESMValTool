#!/usr/bin/env python
"""Python example diagnostic."""
import logging
import os

import cf_units
import iris
import numpy as np

from esmvaltool.diag_scripts.shared.trend_mpqb_common.diag1d import mannkendall1d, theilslopes1d
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvaltool.diag_scripts.shared._base import get_plot_filename
from mpqb_plots import get_ecv_plot_config, mpqb_mapplot
from esmvaltool.diag_scripts.shared.trend_mpqb_common.sharedutils import parallel_apply_along_axis

logger = logging.getLogger(os.path.basename(__file__))


def theilsenmk(cube):
    """Estimate theil-sen slope and mask based on mann-kendall test."""
    template = cube.collapsed('time', iris.analysis.MEAN)
    indata = cube.data.filled(np.nan)
    theilsendata = parallel_apply_along_axis(theilslopes1d, 0, indata)
    outcube = template.copy()
    outcube.data = np.ma.fix_invalid(theilsendata)
    # Now mask based on MK test
    mkdata = parallel_apply_along_axis(mannkendall1d, 0, indata)
    mkmask = (mkdata == 0)  # create mask where MK test is zero
    outcube.data.mask |= mkmask
    # Set name
    outcube.rename('theil-sen trend of ' + outcube.name())
    outcube.units = cf_units.Unit(str(outcube.units) + ' year-1')
    return outcube


def timemean(cube):
    """Calculate mean over time."""
    fullmean = cube.collapsed('time', iris.analysis.MEAN)
    return fullmean


def _array2cube(array_in, cube_template):
    newcube = cube_template.copy()
    newcube.data = np.ma.fix_invalid(array_in)
    return newcube


def compute_diagnostic(filename):
    """Compute an example diagnostic."""
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)

    logger.debug("Running example computation")
    return cube.collapsed('time', iris.analysis.MEAN)


def main(cfg):
    """Compute the time average for each input dataset."""
    metrics_to_calculate = ['theilsenmk', 'timemean']

    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    grouped_input_data = group_metadata(
        input_data, 'dataset', sort='dataset')

    # Loop through all datasets
    for dataset in grouped_input_data.keys():
        dataset_cfg = grouped_input_data[dataset]

        logger.info("Opening dataset: %s", dataset)
        # Opening the pair
        cube = iris.load_cube(dataset_cfg[0]['filename'])

        for metricname in metrics_to_calculate:
            try:
                # Functions are defined in this script
                resultcube = globals()[metricname](cube)
            except KeyError:
                logger.error("Metric %s is not defined. ", metricname)
                continue
            # Plot the results (if configured to plot)
            plot_file = get_plot_filename(metricname + '_' + dataset, cfg)
            if cfg['write_plots']:
                metrics_plot_dictionary = get_ecv_plot_config(
                    dataset_cfg[0]['short_name'])
                mpqb_mapplot(resultcube, plot_file, **
                             metrics_plot_dictionary[metricname])
            logger.info("Finished aux plots for dataset: %s", dataset)
    logger.info("Finished!")


if __name__ == '__main__':
    with run_diagnostic() as global_cfg:
        main(global_cfg)
