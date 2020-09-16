#!/usr/bin/env python
"""Python example diagnostic."""
import logging
import os

import cf_units
import iris
import numpy as np
import yaml

from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvaltool.diag_scripts.shared._base import get_plot_filename
from esmvaltool.diag_scripts.shared.trend_mpqb_common.diag1d import (
    mannkendall1d, theilslopes1d)
from esmvaltool.diag_scripts.shared.trend_mpqb_common.sharedutils import \
    parallel_apply_along_axis
from mpqb_plots import mpqb_mapplot
from mpqb_utils import get_mpqb_cfg

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
    mean_days = np.mean(np.diff(cube.coord('time').points))
    days_per_year = 365.25
    outcube.data = outcube.data * (days_per_year / mean_days)
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


def main(cfg):
    """Compute the time average for each input dataset."""
    metrics_to_calculate = ['theilsenmk', 'timemean']

    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    grouped_input_data = group_metadata(input_data, 'alias', sort='alias')

    # Loop through all datasets
    for alias in grouped_input_data.keys():
        dataset_cfg = grouped_input_data[alias][0]
        dataset = grouped_input_data[alias][0]['dataset']

        logger.info("Opening dataset: %s", dataset)
        # Opening the pair
        cube = iris.load_cube(dataset_cfg['filename'])

        for metricname in metrics_to_calculate:
            try:
                # Functions are defined in this script
                resultcube = globals()[metricname](cube)
            except KeyError:
                logger.error("Metric %s is not defined. ", metricname)
                continue
            # Plot the results (if configured to plot)
            if cfg['write_plots']:
                baseplotname = f"{alias}_{metricname}" \
                               f"_{dataset_cfg['variable_group']}" \
                               f"_{dataset_cfg['start_year']}-" \
                               f"{dataset_cfg['end_year']}"
                plot_file = get_plot_filename(baseplotname, cfg)
                metrics_plot_dictionary = get_mpqb_cfg('colormap',
                    dataset_cfg['short_name'])
                plot_kwargs = metrics_plot_dictionary[metricname]
                # Overwrite plot title to be dataset name
                plot_kwargs['title'] = alias
                # Specify to add small text with field mean for timemean
                if metricname == 'timemean':
                    plot_kwargs['addglobmeanvalue'] = True
                mpqb_mapplot(resultcube, cfg, plot_file, **plot_kwargs)
            logger.info("Finished aux plots for dataset: %s", dataset)
    logger.info("Finished!")




if __name__ == '__main__':
    with run_diagnostic() as global_cfg:
        main(global_cfg)
