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
from mpqb_plots import read_mpqb_cfg

logger = logging.getLogger(os.path.basename(__file__))


def _calculate_histograms(cube, numbars, lower_upper):
    # returns a dictionary for histogram plotting
<<<<<<< HEAD
    hist, bins = da.histogram(cube.core_data(), bins=numbars, range=lower_upper)
    print(hist)
    hist = hist.compute()
    print(hist)
    hist = hist/np.sum(hist)
    print(hist)
    print(bins)
=======
    hist, bins = da.histogram(cube.core_data(),
                              bins=numbars,
                              range=lower_upper)
    hist = hist.compute()
    hist = hist / np.sum(hist)
>>>>>>> d8e576043103bff703315f0a7266e5d9a3421f21
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


    mpqb_cfg = read_mpqb_cfg()
    datasetnames = mpqb_cfg['datasetnames']

    plt.clf()
    fig, ax1 = plt.subplots(1, 1)

    for alias, hist in hists.items():
        dataset_cfg = grouped_input_data[alias][0]
        dataset = dataset_cfg['dataset']


        width = np.diff(hist["bins"]) / len(hists)
        xvals = hist["bins"][:-1] + width * (
            list(hists.keys()).index(alias) + 0.5)

        if histtype == 'bar':
            ax1.bar(
                xvals,
                hist["hist"],
                width,
                label=datasetnames[alias],
                color=mpqb_cfg['datasetcolors'][alias],
            )
            plt.vlines(hist["bins"], 0, 1, linestyles='dashed', alpha=0.3)
            plt.legend()
        elif histtype == 'step':
            ax1.hist(hist["bins"][:-1],
                     hist["bins"],
                     weights=hist["hist"],
                     label=datasetnames[alias],
                     histtype='step',
                     color=mpqb_cfg['datasetcolors'][alias],
                     linewidth=2)
            handles, labels = ax1.get_legend_handles_labels()
            handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
            plt.legend(handles=handles, labels=labels)
        else:
            logger.warning("Unsupported argument for histtype: %s", histtype)

    plt.xlabel(grouped_input_data[alias][0]['long_name'] + " [" +
               grouped_input_data[alias][0]['units'] + "]")
    plt.ylabel("frequency")
    if cfg.pop('logarithmic', False):
        ax1.set_yscale('log', nonposy='clip')

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
<<<<<<< HEAD
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
                
            if histtype=='bar':
                ax.bar(x,
                       hist["hist"],
                       width,
                       label = grouped_input_data[dataset][0]['alias'],
                       )
                plt.vlines(hist["bins"], 0, 1, linestyles='dashed', alpha = 0.3)
                plt.legend()
            elif histtype=='step':
                ax.hist(hist["bins"][:-1],
                        hist["bins"],
                        weights = hist["hist"],
                        label = grouped_input_data[dataset][0]['alias'],
                        histtype = 'step',
                        linewidth = 2
                        )
                handles, labels = ax.get_legend_handles_labels()
                new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
                plt.legend(handles=new_handles, labels=labels)
            else:
                logger.warning(f"Unsupported argument for histtype: {histtype}")

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
        plt.ylim(y0,y1)
        plt.tight_layout()
#        filename = get_plot_filename('histogram', cfg)
        filename = get_plot_filename('histogram_' + cfg["script"], cfg)
        logger.info("Saving as %s", filename)
        fig.savefig(filename)
        plt.close(fig)
    else:
        logger.warning("This diagnostic wants to plot, but isn't allowed to")
=======
    for alias in grouped_input_data:
        dataset_cfg = grouped_input_data[alias][0]

        logger.info("Opening dataset: %s", alias)
        cube = iris.load_cube(grouped_input_data[alias][0]['filename'])
        hists.update(
            {alias: _calculate_histograms(cube, nbars, lower_upper)})
    logger.info("Plotting the histograms.")
    _plot_histograms(hists, cfg, grouped_input_data)
>>>>>>> d8e576043103bff703315f0a7266e5d9a3421f21
    logger.info("Finished!")


if __name__ == '__main__':
<<<<<<< HEAD
    with run_diagnostic() as cfg:
        main(cfg)
=======
    with run_diagnostic() as config:
        if config['write_plots']:
            main(config)
        else:
            logger.warning("This diagnostic wants to plot,\
                            but isn't allowed to")
>>>>>>> d8e576043103bff703315f0a7266e5d9a3421f21
