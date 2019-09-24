"""Python example diagnostic."""
import logging
import os
from pprint import pformat
from sharedutils import parallel_apply_along_axis
from diag1d import theilslopes1d, mannkendall1d
from mpqb_plots import mpqb_mapplot, get_ecv_plot_config
import numpy as np
from esmvalcore.preprocessor._time import annual_statistics
from esmvaltool.diag_scripts.shared._base import get_diagnostic_filename


import cf_units
import iris
import itertools as it
import warnings
from basaux import freqchecker

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata, sorted_metadata)
from esmvaltool.diag_scripts.shared._base import (
    ProvenanceLogger, get_diagnostic_filename, get_plot_filename)
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))


def theilsenmk(cube):
    template = cube.collapsed('time', iris.analysis.MEAN)
    indata = cube.data.filled(np.nan)
    theilsendata = parallel_apply_along_axis(theilslopes1d, 0, indata)
    outcube = template.copy()
    outcube.data = np.ma.fix_invalid(theilsendata)
    # Now mask based on MK test
    mkdata = parallel_apply_along_axis(mannkendall1d, 0, indata)
    mkmask = (mkdata == 0) # create mask where MK test is zero
    outcube.data.mask |= mkmask
    # Set name
    outcube.rename('theil-sen trend of '+outcube.name())
    outcube.units = cf_units.Unit(str(outcube.units)+' year-1')
    return outcube


def timemean(cube):
    fullmean = cube.collapsed('time', iris.analysis.MEAN)
    return fullmean

def array2cube(array_in,cube_template):
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

    #TODO move these parameters to config file
    reference_dataset = 'CDS-SATELLITE-SOIL-MOISTURE'
    metrics_to_calculate = ['theilsenmk', 'timemean']

    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    grouped_input_data = group_metadata(
        input_data, 'dataset', sort='dataset')

    # Loop through all datasets
    for dataset in grouped_input_data.keys():
        dataset_cfg = grouped_input_data[dataset]
        
        logger.info("Opening dataset: {0}".format(dataset))
        # Opening the pair 
        cube = iris.load_cube(dataset_cfg[0]['filename'])
        
        for metricname in metrics_to_calculate:
            try:
                resultcube = globals()[metricname](cube) # Functions are defined in this script
            except KeyError:
                logger.error("Metric %s is not defined. ", metricname)
                continue
            # Plot the results (if configured to plot)
            plot_file = get_plot_filename(metricname+'_'+dataset, cfg)
            if cfg['write_plots']:
                metrics_plot_dictionary = get_ecv_plot_config(dataset_cfg[0]['short_name'])
                mpqb_mapplot(resultcube, plot_file, **metrics_plot_dictionary[metricname])
            logger.info("Finished aux plots for dataset: {0}".format(dataset))
    logger.info("Finished!")


if __name__ == '__main__':
    with run_diagnostic() as cfg:
        main(cfg)
