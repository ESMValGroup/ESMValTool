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

import iris
import matplotlib.pyplot as plt
import numpy as np

import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            get_plot_filename, group_metadata,
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


def plot_bar_deangelis(cfg, data_var_sum, available_exp, available_vars):
    """Plot linear regression used to calculate ECS."""
    if not cfg[n.WRITE_PLOTS]:
        return

    # Plot data
    fig, axx = plt.subplots()

    set_colors = ['cornflowerblue', 'orange', 'silver', 'limegreen',
                  'rosybrown', 'orchid']
    bar_width = 1.0 / float(len(available_vars))

    for iii, iexp in enumerate(available_exp):
        axx.bar(np.arange(len(available_vars)) + bar_width * float(iii),
                data_var_sum[iexp],
                bar_width, color=set_colors[iii], label=iexp)

    axx.set_xlabel(' ')
    axx.set_ylabel(r'Model mean (W m$^{-2}$)')
    axx.set_title(' ')
    axx.set_xticks(np.arange(len(available_vars)) + bar_width)
    axx.set_xticklabels(available_vars)
    axx.legend(loc=1)

    fig.tight_layout()
    fig.savefig(get_plot_filename('bar_all', cfg), dpi=300)
    plt.close()

    caption = 'Global average multi-model mean comparing different ' + \
              'model experiments and flux variables.'

    provenance_record = get_provenance_record(
        _get_sel_files_var(cfg, available_vars), caption, ['mean'], ['global'])

    diagnostic_file = get_diagnostic_filename('bar_all', cfg)

    logger.info("Saving analysis results to %s", diagnostic_file)

    list_dict = {}
    list_dict["data"] = []
    list_dict["name"] = []
    for iexp in available_exp:
        list_dict["data"].append(data_var_sum[iexp])
        list_dict["name"].append({'var_name': iexp + '_all',
                                  'long_name': 'Fluxes for ' + iexp +
                                               ' experiment',
                                  'units': 'W m-2'})

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
    # var = e.Variables(cfg)
    available_vars = list(group_metadata(cfg['input_data'].values(),
                                         'short_name'))
    logging.debug("Found variables in recipe:\n%s", available_vars)

    available_exp = list(group_metadata(cfg['input_data'].values(), 'exp'))

    if len(available_exp) > 6:
        raise ValueError("The diagnostic can only plot up to 6 different " +
                         "model experiments.")

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

    data_var = OrderedDict()
    for iexp in available_exp:
        data_var[iexp] = OrderedDict()
        for jvar in available_vars:
            # data_var[iexp] = OrderedDict()
            data_var[iexp][jvar] = 0.0

    pathlist = data.get_path_list(short_name=available_vars[0],
                                  exp=available_exp[0])

    for dataset_path in pathlist:

        # Substract piControl experiment from abrupt4xCO2 experiment
        dataset = data.get_info(n.DATASET, dataset_path)

        for jvar in available_vars:
            for iexp in available_exp:
                print(data_var[iexp])
                print((data_var[iexp].values()))
                (data_var[iexp])[jvar] = (data_var[iexp])[jvar] + \
                    data.get_data(short_name=jvar, exp=iexp,
                                  dataset=dataset)

    data_var_sum = {}
    for iexp in available_exp:
        data_var_sum[iexp] = np.fromiter(data_var[iexp].values(),
                                         dtype=float) / float(len(pathlist))

    # Plot ECS regression if desired
    plot_bar_deangelis(cfg, data_var_sum, available_exp, available_vars)


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
