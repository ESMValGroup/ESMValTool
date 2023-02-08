#!/usr/bin/env python
"""Python example diagnostic."""
import logging
import os
from pprint import pformat

import dask.array as da
import iris
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvaltool.diag_scripts.shared._base import get_plot_filename
from mpqb_utils import get_mpqb_cfg

logger = logging.getLogger(os.path.basename(__file__))


def _calculate_histograms(cube, numbars, lower_upper):
    # returns a dictionary for histogram plotting
    weights = iris.analysis.cartography.area_weights(cube, normalize=True)
    # Now convert to dask array, needed for da.histogram
    weights = da.from_array(weights, chunks=cube.core_data().chunks)
    hist, bins = da.histogram(cube.core_data(),
                              bins=numbars,
                              range=lower_upper, weights=weights)
    hist = hist.compute()
    hist = hist / np.sum(hist)
    return {"hist": hist, "bins": bins}


def _get_lower_upper(input_data):
    lower_upper = [np.inf, -np.inf]
    for dataset in input_data:
        logger.info("Opening dataset: %s", dataset)
        cube = iris.load_cube(input_data[dataset][0]['filename'])
        l_u = [
            da.min(cube.core_data().flatten()).compute(),
            da.max(cube.core_data().flatten()).compute()
        ]
        if l_u[0] < lower_upper[0]:
            lower_upper[0] = l_u[0]
        if l_u[1] > lower_upper[1]:
            lower_upper[1] = l_u[1]
    return lower_upper


def _plot_histograms(hists, cfg, grouped_input_data):
    histtype = cfg.pop('histtype', 'bar')  # default type is 'bar'
    dataset = None  # needed to avoid pylint error

    plt.clf()
    fig, ax1 = plt.subplots(1, 1)

    for alias, hist in hists.items():
        dataset_cfg = grouped_input_data[alias][0]
        dataset = dataset_cfg['dataset']


        width = np.diff(hist["bins"]) / len(hists)
        xvals = hist["bins"][:-1] + width * (
            list(hists.keys()).index(alias) + 0.5)

        label = get_mpqb_cfg('datasetname', alias)
        color = get_mpqb_cfg('datasetcolor', alias)
        if histtype == 'bar':
            ax1.bar(
                xvals,
                hist["hist"],
                width,
                label=label,
                color=color,
            )
            plt.vlines(hist["bins"], 0, 1, linestyles='dashed', alpha=0.3)
            plt.legend()
        elif histtype == 'step':
            ax1.hist(hist["bins"][:-1],
                     hist["bins"],
                     weights=hist["hist"],
                     label=label,
                     histtype='step',
                     color=color,
                     linewidth=2)
            handles, labels = ax1.get_legend_handles_labels()
            handles = [Line2D([], [], c=h.get_edgecolor(), lw=2.0) for h in handles]
            plt.legend(handles=handles, labels=labels)
        else:
            logger.warning("Unsupported argument for histtype: %s", histtype)
    ax1.tick_params(axis='both', which='major', labelsize='large')
    plt.xlabel(grouped_input_data[alias][0]['long_name'] + " [" +
               grouped_input_data[alias][0]['units'] + "]", fontsize='large')
    plt.ylabel("frequency",  fontsize='x-large')
    if cfg.pop('logarithmic', False):
        ax1.set_yscale('log', nonposy='clip')
    ax1.tick_params(axis='both', which='major', labelsize=12)

    plt.ylim(
        0,
        1.1 * np.max([np.max(content["hist"])
                      for _, content in hists.items()]))
    plt.ylim(cfg.pop('y0', None), cfg.pop('y1', None))
    plt.tight_layout()
    filename = get_plot_filename('histogram', cfg)
    logger.info("Saving as %s", filename)
    fig.savefig(filename)
    plt.close(fig)


def main(cfg):
    """Plot a histogram."""
    # Get a description of the preprocessed data that we will use as input.
    lower_upper = [cfg.pop('vmin', None), cfg.pop('vmax', None)]
    nbars = cfg.pop('number_of_bars', 10)

    grouped_input_data = group_metadata(cfg['input_data'].values(),
                                        'alias',
                                        sort='alias')

    # calculate common range if vmin or vmax not given
    if lower_upper[0] is None or lower_upper[1] is None:
        logger.info("Common range for histogram is calculated")
        lower_upper = _get_lower_upper(grouped_input_data)

    hists = dict()
    logger.info("Calculating the histograms.")
    for alias in grouped_input_data:
        dataset_cfg = grouped_input_data[alias][0]

        logger.info("Opening dataset: %s", alias)
        cube = iris.load_cube(grouped_input_data[alias][0]['filename'])
        hists.update(
            {alias: _calculate_histograms(cube, nbars, lower_upper)})
    logger.info("Plotting the histograms.")
    _plot_histograms(hists, cfg, grouped_input_data)
    logger.info("Finished!")


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
