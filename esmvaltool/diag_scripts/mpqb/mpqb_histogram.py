#!/usr/bin/env python
"""Python example diagnostic."""
import logging
import os
from pprint import pformat
from esmvaltool.diag_scripts.shared.trend_mpqb_common.sharedutils import parallel_apply_along_axis
from esmvaltool.diag_scripts.shared.trend_mpqb_common.diag1d import *
import numpy as np
import dask.array as da

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

def calculate_histograms(cube, numbars, lower_upper):
    # returns a dictionary for histogram plotting
    hist, bins = da.histogram(cube.core_data(), bins=numbars, range=lower_upper)
    hist = hist.compute()
    hist = hist/np.sum(hist)
    return {"hist": hist, "bins": bins}
    
def calculate_min_max(cube):
    # returns a list of the minimum and maximum value of a cube
    vmin = da.min(cube.core_data().flatten()).compute()
    vmax = da.max(cube.core_data().flatten()).compute()
    return [vmin, vmax]

def main(cfg):
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()
    number_of_bars = cfg.pop('number_of_bars', 10) # default number of bins set to 10
    lower_upper = [cfg.pop('vmin', None), cfg.pop('vmax', None)]
    logarithmic = cfg.pop('logarithmic', False)
    grouped_input_data = group_metadata(
        input_data, 'dataset', sort='dataset')
    logger.info(
        "Example of how to group and sort input data by standard_name:"
        "\n%s", pformat(grouped_input_data))
    
    # calculate common range if vmin or vmax not given
    if lower_upper[0] is None or lower_upper[1] is None:
        logger.info("Common range for histogram is calculated " + 
                    "due to vmin and/or vmax specification missing in recipe.")
        lower_upper = [np.inf, -np.inf]
        for dataset in grouped_input_data:
            logger.info("Opening dataset: {0}".format(dataset))
            cube = iris.load_cube(grouped_input_data[dataset][0]['filename'])
            l_u = calculate_min_max(cube)
            if l_u[0] < lower_upper[0]:
                lower_upper[0] = l_u[0]
            if l_u[1] > lower_upper[1]:
                lower_upper[1] = l_u[1]
                
    hists = dict()
    logger.info("Calculating the histograms.")
    for dataset in grouped_input_data:
        logger.info("Opening dataset: {0}".format(dataset))
        cube = iris.load_cube(grouped_input_data[dataset][0]['filename'])    
        hists.update({dataset: calculate_histograms(cube,
                                                    number_of_bars, 
                                                    lower_upper)})

    if cfg['write_plots']:
        plt.clf()
        fig, ax = plt.subplots(1,1)
        for dataset, hist in hists.items():
            
            width = np.diff(hist["bins"]) / len(hists)
            x = hist["bins"][:-1] + width * (list(
                    hists.keys()).index(dataset) + 0.5)
                
            ax.bar(x,
                   hist["hist"],
                   width,
                   label = dataset,
                   )
            
        plt.vlines(hist["bins"], 0, 1, linestyles='dashed', alpha = 0.3)
        plt.legend()
        plt.xlabel(grouped_input_data[dataset][0]['long_name'] +
                   " [" + grouped_input_data[dataset][0]['units'] + "]")
        plt.ylabel("frequency")
        if logarithmic:
            ax.set_yscale('log', nonposy='clip')
#        if logarithmic:
#            plt.ylim(0, 5*np.max(
#                     [np.max(content["hist"]) for _, content in hists.items()]))
#        else:
        plt.ylim(0, 1.1*np.max(
                 [np.max(content["hist"]) for _, content in hists.items()]))
        plt.tight_layout()
#        filename = get_plot_filename('histogram', cfg)
        filename = get_plot_filename('histogram_' + cfg["script"], cfg)
        logger.info("Saving as %s", filename)
        fig.savefig(filename)
        plt.close(fig)
    else:
        logger.warning("This diagnostic wants to plot, but isn't allowed to")
    logger.info("Finished!")


if __name__ == '__main__':
    with run_diagnostic() as cfg:
        main(cfg)