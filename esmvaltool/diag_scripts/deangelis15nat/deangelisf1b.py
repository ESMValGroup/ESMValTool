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
from pprint import pformat
import matplotlib.pyplot as plt
import iris
import numpy as np
import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared._base import (
    ProvenanceLogger, get_diagnostic_filename, get_plot_filename,
    select_metadata)

logger = logging.getLogger(os.path.basename(__file__))


def _get_sel_files_var(cfg, varnames):
    """Get filenames from cfg for all model mean and differen variables."""
    selection = []

    for var in varnames:
        for hlp in select_metadata(cfg['input_data'].values(), short_name=var):
            selection.append(hlp['filename'])

    return selection


def cube_to_save_vars(list_dict):
    """Create cubes to prepare bar plot data for saving to netCDF."""
    # cubes = iris.cube.CubeList()
    for iii, var in enumerate(list_dict["data"]):
        if iii == 0:
            cubes = iris.cube.CubeList([
                iris.cube.Cube(var,
                               var_name=list_dict["name"][iii]['var_name'],
                               long_name=list_dict["name"][iii]['long_name'],
                               units=list_dict["name"][iii]['units'])])
        else:
            cubes.append(
                iris.cube.Cube(var,
                               var_name=list_dict["name"][iii]['var_name'],
                               long_name=list_dict["name"][iii]['long_name'],
                               units=list_dict["name"][iii]['units']))

    return cubes


def get_provenance_record(ancestor_files, caption, statistics,
                          domains, plot_type='bar'):
    """Get Provenance record."""
    record = {
        'caption': caption,
        'statistics': statistics,
        'domains': domains,
        'plot_type': plot_type,
        'themes': ['phys'],
        'authors': [
            'weigel_katja',
        ],
        'references': [
            'deangelis15nat',
        ],
        'ancestors': ancestor_files,
    }
    return record


def plot_bar_deangelis(cfg, pi_all, rcp85_all, co2_all):
    """Plot linear regression used to calculate ECS."""
    if not cfg[n.WRITE_PLOTS]:
        return

    # Plot data
    n_groups = 4

    fig, axx = plt.subplots()

    index = np.arange(n_groups)
    bar_width = 0.25

    axx.bar(index, pi_all, bar_width, color='cornflowerblue',
            label='Pre-industrial')
    axx.bar(index + bar_width, rcp85_all, bar_width, color='orange',
            label='RCP8.5')
    axx.bar(index + (bar_width * 2), co2_all, bar_width, color='silver',
            label=r'4xCO$_2$')

    axx.set_xlabel(' ')
    axx.set_ylabel(r'Model mean (W m$^{-2}$)')
    axx.set_title(' ')
    axx.set_xticks(index + bar_width)
    axx.set_xticklabels(('lvp', 'rlnst', 'rsnst', 'hfss'))
    axx.legend(loc=1)

    fig.tight_layout()
    fig.savefig(get_plot_filename('bar_all', cfg), dpi=300)
    plt.close()

    caption = '...'

    provenance_record = get_provenance_record(_get_sel_files_var(cfg,
                                                                 ['rlnst',
                                                                  'rsnst',
                                                                  'lvp',
                                                                  'hfss']),
                                              caption, ['mean'], ['global'])

    diagnostic_file = get_diagnostic_filename('bar_all', cfg)

    logger.info("Saving analysis results to %s", diagnostic_file)

    list_dict = {}
    list_dict["data"] = [pi_all, rcp85_all, co2_all]
    list_dict["name"] = [{'var_name': 'pi_all',
                          'long_name': 'Fluxes for piControl scenario',
                          'units': 'W m-2'},
                          {'var_name': 'rcp85_all',
                          'long_name': 'Fluxes for future scenario',
                          'units': 'W m-2'},
                          {'var_name': 'co2_all',
                          'long_name': 'Fluxes for 4 times CO2 scenario',
                          'units': 'W m-2'}]

    iris.save(cube_to_save_vars(list_dict), target=diagnostic_file)

    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                pformat(provenance_record))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)


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
    plot_bar_deangelis(cfg, pi_all, rcp85_all, co2_all)


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
