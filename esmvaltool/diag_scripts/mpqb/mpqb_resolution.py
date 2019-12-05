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
import matplotlib.cm as mpl_cm
import iris
import itertools as it
import warnings
from pprint import pformat

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
    
    # calculate resolution
    logger.info("Calculating the resolution.")
    datasets = []
    names = []
    for dataset in grouped_input_data:
        logger.info("Opening dataset: {0}".format(dataset))
        datasets.append(
                iris.load_cube(grouped_input_data[dataset][0]['filename']))
        names.append(dataset)
        
    res_data = __get_tim_res__(datasets, names)
    
    logger.info("The extracted resolution information is:")
    logger.info(pformat(res_data))

    if cfg['write_plots']:
        plt.clf()
        fig, ax = plt.subplots(1,1)
        fig.set_figwidth(1.3*fig.get_figwidth())
        ax.plot()
        plt.quiver([0] * len(res_data["units"]),
                   range(len(res_data["units"])),
                   [1.15] * len(res_data["units"]),
                   [0] * len(res_data["units"]),
                   angles='xy', scale_units='xy', scale=0.95)
        plt.quiver([1.] * len(res_data["units"]),
                   range(len(res_data["units"])),
                   [-1.15] * len(res_data["units"]),
                   [0] * len(res_data["units"]),
                   angles='xy', scale_units='xy', scale=0.95)
        plt.xlim(-0.6, 1.25)
        plt.ylim(-1, len(res_data["units"]))

        coords = res_data.copy()
        coords.pop("names")
        coords.pop("units")
        arrow = 0.
        cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
        for c in coords.keys():
            plt.text(-0.6, arrow, c + " [" + res_data["units"][c] + "]")
            plot_info = shinescale_0_1(res_data[c])
            plt.scatter(plot_info[1],
                        [arrow] * len(plot_info[1]),
                        color="black",
                        marker='|',
                        s=1000,
                        alpha=1.)
            
            for l, v in enumerate(plot_info[1]):
                plt.text(v, arrow + 0.3, str(plot_info[2][l]), ha="center")
                
            plt.scatter(list(reversed(plot_info[0])),
                        [arrow - 0.15] * len(plot_info[0]),
                        color=list(reversed(cycle[0:len(plot_info[0])])),
                        marker='^',
                        s=list(reversed(60+50*np.arange(len(plot_info[0])))),
                        )
            arrow += 1

        for ind, lab in enumerate(res_data["names"]):
            plt.plot(-100, -100, 
                     marker="^",
                     markersize=7+1.8*ind,
                     linewidth=0,
                     label=lab,
                     )

        ax.axis("off")

        box = ax.get_position()
        ax.set_position([box.x0, box.y0,
                         box.width, box.y1 - box.height * 0.1])

        # Put a legend below current axis
        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=3, mode="expand", borderaxespad=0.)

        plt.tight_layout()
        filename = get_plot_filename('resolution', cfg)
        logger.info("Saving as %s", filename)
        fig.savefig(filename)
        plt.close(fig)
    else:
        logger.warning("This diagnostic wants to plot, but isn't allowed to")
    logger.info("Finished!")

def shinescale_0_1(l):
    l0 = np.array(l)

    vmin = np.min(l0)
    vmax = np.max(l0)
    diff = vmax - vmin
    if diff == 0:
        diff = 1.
    rounder = int(np.ceil(-np.log10(diff)))
    vmin = (np.floor(vmin * 10**rounder) - 1) / 10**rounder
    vmax = (np.ceil(vmax * 10**rounder) + 1) / 10**rounder
    levels = np.round(np.linspace(vmin, vmax, num=4), rounder)

    l0 = (l0 - vmin) / (vmax - vmin)

    labels = list(levels)
    levels = (levels - np.min(levels)) / \
        (np.max(levels) - np.min(levels))

    return((list(l0), levels, labels))
    
def __get_tim_res__(datasets, names):
        
    tim_res_dict = {}
    # We assume respective coords of all datasets in the same unit!
    unitdict={}
    
    loc_mpdata = datasets.copy()
    [__remove_all_aux_coords__(mpdats) for mpdats in loc_mpdata]
    
    coords = list(set("/".join([c.name() for c in mpd.coords()]) 
        for mpd in datasets))[0].split("/")
    
    for c in coords:

        tim_res_dict.update({c:[np.mean(np.diff(mpd.coord(c).points)) 
            for mpd in datasets]})
        unitdict.update({c:[str(datasets[0].coord(c).units) 
            if not datasets[0].coord(c).units.is_time_reference() 
            else str(datasets[0].coord(c).units).split(" ")[0]][0]})
        
    tim_res_dict.update({"names":names})
    
    tim_res_dict.update({"units":unitdict})
    
    return tim_res_dict

def __remove_all_aux_coords__(cube):
    """
    remove all auxiliary coordinates from cube
    """
    for dim in cube.coords():
        if isinstance(dim, iris.coords.AuxCoord):
             cube.remove_coord(dim)
            
    return

if __name__ == '__main__':
    with run_diagnostic() as cfg:
        main(cfg)