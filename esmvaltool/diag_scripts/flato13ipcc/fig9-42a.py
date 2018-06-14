#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###############################################################################
flato13ipcc/fig9-42a.py
Author: Manuel Schlund (DLR, Germany)
CRESCENDO project
###############################################################################

Description
    Calculate and plot the equilibrium climate sensitivity (ECS) vs. the global
    mean surface temperature (GMSAT) for several CMIP5 models (see IPCC AR5 WG1
    ch. 9, fig. 9.42a).

Configuration options
    plot_ecs_regression : Switch to plot the linear regressions needed for the
                          ECS calculations

Optional diag_script_info attributes (diagnostic specific)
    None

Caveats

Modification history
    20180522-A_schl_ma: ported to v2.0
    20171109-A_schl_ma: written

###############################################################################
"""


from esmvaltool.diag_scripts.shared import *

import iris

from collections import OrderedDict
from scipy import stats
import logging
import numpy as np
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    """
    Arguments
        cfg : Dictionary containing project information

    Description
        This is the main routine of the diagnostic.
    """

    ###########################################################################
    # Setup diagnostic
    ###########################################################################

    logging.info(cfg)

    # Model data containers
    MODELS = Models(cfg)
    logging.info("Found models:\n{}".format(MODELS))

    # Variables
    VARS = Variables(cfg)
    logging.info("Found variables:\n{}".format(VARS))

    # Experiments
    PICONTROL = 'piControl'
    HISTORICAL = 'historical'
    ABRUPT4XCO2 = 'abrupt4xCO2'

    # Matplotlib instance
    fig, axes = plt.subplots()

    ###########################################################################
    # Read data
    ###########################################################################

    # Create iris cube for each model
    for model_path in MODELS:
        cube = iris.load(model_path, VARS.standard_names())[0]

        # Global annual mean
        for coord in [cube.coord(LAT), cube.coord(LON)]:
            if (not coord.has_bounds()):
                coord.guess_bounds()
        area_weights = iris.analysis.cartography.area_weights(cube)
        cube = cube.collapsed([LAT, LON], iris.analysis.MEAN,
                              weights=area_weights)
        cube = cube.aggregated_by(YEAR, iris.analysis.MEAN)
        MODELS.set_data(cube.data, path_to_model=model_path)

    print(MODELS.get_data())

    ###########################################################################
    # Process data
    ###########################################################################

    ###########################################################################
    # Plot data
    ###########################################################################


if __name__ == '__main__':

    with run_diagnostic() as cfg:
        main(cfg)
