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

logger = logging.getLogger(os.path.basename(__file__))

DATASET_PLOTNAMES = {
    'ERA-Interim-Land': 'ERA-Interim/Land',
    'CDS-SATELLITE-SOIL-MOISTURE': 'ESA-CCI',
    'cds-era5-land-monthly': 'ERA5-Land',
    'cds-era5-monthly': 'ERA5',
    'MERRA2': 'MERRA-2',
    'cds-satellite-lai-fapar': 'SPOT-VGT',
}


def _calculate_histograms(cube, numbars, lower_upper):
    # returns a dictionary for histogram plotting
    hist, bins = da.histogram(
        cube.core_data(), bins=numbars, range=lower_upper)
    hist = hist.compute()
    hist = hist / np.sum(hist)
    return {"hist": hist, "bins": bins}


def _get_lower_upper(input_data):
    lower_upper = [np.inf, -np.inf]
    for dataset in input_data:
        logger.info("Opening dataset: %s", dataset)
        cube = iris.load_cube(input_data[dataset][0]['filename'])
        l_u = [da.min(cube.core_data().flatten()).compute(),
               da.max(cube.core_data().flatten()).compute()]
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

    for dataset, hist in hists.items():

        width = np.diff(hist["bins"]) / len(hists)
        xvals = hist["bins"][:-1] + width * (list(
            hists.keys()).index(dataset) + 0.5)

        if histtype == 'bar':
            ax1.bar(xvals,
                    hist["hist"],
                    width,
                    label=DATASET_PLOTNAMES[dataset],
                    )
            plt.vlines(hist["bins"], 0, 1, linestyles='dashed', alpha=0.3)
            plt.legend()
        elif histtype == 'step':
            ax1.hist(hist["bins"][:-1],
                     hist["bins"],
                     weights=hist["hist"],
                     label=DATASET_PLOTNAMES[dataset],
                     histtype='step',
                     linewidth=2
                     )
            handles, labels = ax1.get_legend_handles_labels()
            handles = [Line2D([], [], c=h.get_edgecolor())
                       for h in handles]
            plt.legend(handles=handles, labels=labels)
        else:
            logger.warning(
                "Unsupported argument for histtype: %s", histtype)

    plt.xlabel(grouped_input_data[dataset][0]['long_name'] +
               " [" + grouped_input_data[dataset][0]['units'] + "]")
    plt.ylabel("frequency")
    if cfg.pop('logarithmic', False):
        ax1.set_yscale('log', nonposy='clip')

    plt.ylim(0, 1.1 * np.max(
        [np.max(content["hist"]) for _, content in hists.items()]))
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

    grouped_input_data = group_metadata(
        cfg['input_data'].values(), 'dataset', sort='dataset')
    logger.info(
        "Example of how to group and sort input data by standard_name:"
        "\n%s", pformat(grouped_input_data))

    # calculate common range if vmin or vmax not given
    if lower_upper[0] is None or lower_upper[1] is None:
        logger.info("Common range for histogram is calculated")
        lower_upper = _get_lower_upper(grouped_input_data)

    hists = dict()
    logger.info("Calculating the histograms.")
    for dataset in grouped_input_data:
        logger.info("Opening dataset: %s", dataset)
        cube = iris.load_cube(grouped_input_data[dataset][0]['filename'])
        hists.update({dataset: _calculate_histograms(cube,
                                                     nbars,
                                                     lower_upper)})
    logger.info("Plotting the histograms.")
    _plot_histograms(hists, cfg, grouped_input_data)
    logger.info("Finished!")


if __name__ == '__main__':
    with run_diagnostic() as config:
        if config['write_plots']:
            main(config)
        else:
            logger.warning("This diagnostic wants to plot,\
                            but isn't allowed to")
