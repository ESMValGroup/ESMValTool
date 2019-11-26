#!/usr/bin/env python
# -*- coding: utf-8 -*-


""".

Calculates radiative constraint on hydrologic cycle intensification
following DeAngelis et al. (2015).


###############################################################################
testkw/deangelis1b.py
Author: Katja Weigel (IUP, Uni Bremen, Germany)
EVal4CMIP project
###############################################################################

Description
-----------

    Calculates radiative constraint on hydrologic cycle intensification
    following DeAngelis et al. (2015).
    Based on diag_scripts/climate_metrics/ecs.py by Manuel Schlund

Configuration options
---------------------
    output_name     : Name of the output files.

###############################################################################

"""


import logging
import os
from collections import OrderedDict
import matplotlib.pyplot as plt
import iris
import numpy as np
import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n

logger = logging.getLogger(os.path.basename(__file__))


def plot_bar_kw(cfg, dd1, dd2, dd3):
    """Plot linear regression used to calculate ECS."""
    if not cfg[n.WRITE_PLOTS]:
        return

    file_name = os.path.join(cfg[n.PLOT_DIR],
                             'bar_all.' + cfg[n.OUTPUT_FILE_TYPE])

    # Plot data
    n_groups = 4

    fig, axx = plt.subplots()

    index = np.arange(n_groups)
    bar_width = 0.25

    axx.bar(index, dd1, bar_width, color='tab:gray',
            label='Pre-industrial')
    axx.bar(index + bar_width, dd2, bar_width, color='g',
            label='RCP8.5')
    axx.bar(index + (bar_width * 2), dd3, bar_width, color='r',
            label=r'4xCO$_2$')

    axx.set_xlabel(' ')
    axx.set_ylabel(r'Model mean (W m$^{-2}$)')
    axx.set_title(' ')
    axx.set_xticks(index + bar_width)
    axx.set_xticklabels((r'L$_{\rm v}$P', 'LWC', 'SWA', 'SH'))
    axx.legend(loc=1)

    fig.tight_layout()
    fig.savefig(file_name)
    plt.close()


###############################################################################
# Setup diagnostic
###############################################################################

# Variables


def main(cfg):
    """Run the diagnostic.

    Parameters :

    ----------
    cfg : dict
        Configuration dictionary of the recipe.

    """
    ###########################################################################
    # Read recipe data
    ###########################################################################

    # Dataset data containers
    data = e.Datasets(cfg)
    logging.debug("Found datasets in recipe:\n%s", data)

    # Variables
    var = e.Variables(cfg)
    logging.debug("Found variables in recipe:\n%s", var)

    # Check for tas and rlnst
    if not var.vars_available('lvp', 'rlnst', 'rsnst', 'hfss'):
        raise ValueError("This diagnostic needs 'lvp', 'rlnst', 'rsnst' and" +
                         " 'hfss' variables")

    ###########################################################################
    # Read data
    ###########################################################################

    # Create iris cube for each dataset and save annual means
    for dataset_path in data:
        cube = iris.load(dataset_path)[0]
        # cube = iris.load(dataset_path, var.standard_names())[0]
        cube = cube.collapsed('time', iris.analysis.MEAN)

        data.set_data(cube.data, dataset_path)

    ###########################################################################
    # Process data
    ###########################################################################

    data_var_pi = OrderedDict()
    data_var_4co2 = OrderedDict()
    data_var_rcp85 = OrderedDict()
    varvar = ['lvp', 'rlnst', 'rsnst', 'hfss']   # var.short_names()

    for jvar in varvar:
        data_var_pi[jvar] = 0.0
        data_var_4co2[jvar] = 0.0
        data_var_rcp85[jvar] = 0.0

    pathlist = data.get_path_list(short_name='lvp', exp='piControl')

    for dataset_path in pathlist:

        # Substract piControl experiment from abrupt4xCO2 experiment
        dataset = data.get_info(n.DATASET, dataset_path)

        for jvar in varvar:
            data_var_pi[jvar] = data_var_pi[jvar] + \
                data.get_data(short_name=jvar,
                              exp='piControl', dataset=dataset)
            data_var_4co2[jvar] = data_var_4co2[jvar] + \
                data.get_data(short_name=jvar, exp='abrupt4xCO2',
                              dataset=dataset)
            data_var_rcp85[jvar] = data_var_rcp85[jvar] + \
                data.get_data(short_name=jvar, exp='rcp85',
                              dataset=dataset)

    pi_all = np.fromiter(data_var_pi.values(), dtype=float) / \
        float(len(pathlist))
    rcp85_all = np.fromiter(data_var_rcp85.values(), dtype=float) / \
        float(len(pathlist))
    co2_all = np.fromiter(data_var_4co2.values(), dtype=float) / \
        float(len(pathlist))

    # Plot ECS regression if desired
    plot_bar_kw(cfg, pi_all, rcp85_all, co2_all)


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
