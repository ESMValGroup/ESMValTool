#!/usr/bin/env python
# -*- coding: utf-8 -*-


""".

Calculates radiative constraint on hydrologic cycle intensification
following DeAngelis et al. (2015).

###############################################################################
testkw/deangelis2.py
Author: Katja Weigel (IUP, Uni Bremen, Germany)
EVal4CMIP project
###############################################################################

Description
-----------
    Calculates radiative constraint on hydrologic cycle intensification
    following DeAngelis et al. (2015).
    Creates figure 2 and extended data figure 1 and figure 2
    Based on diag_scripts/climate_metrics/ecs.py by Manuel Schlund

Configuration options
---------------------

###############################################################################

"""


import logging
import os
from collections import OrderedDict
from pprint import pformat
import iris
import iris.coord_categorisation as cat
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared._base import (
    ProvenanceLogger, get_diagnostic_filename, get_plot_filename,
    select_metadata, group_metadata)

logger = logging.getLogger(os.path.basename(__file__))


def _get_sel_files_var(cfg, varnames):
    """Get filenames from cfg for all model mean and differen variables."""
    selection = []

    for var in varnames:
        for hlp in select_metadata(cfg['input_data'].values(), short_name=var):
            selection.append(hlp['filename'])

    return selection


def cube_to_save_matrix(var1, name):
    """Create cubes to prepare scatter plot data for saving to netCDF."""
    cubes = iris.cube.CubeList([iris.cube.Cube(var1,
                                               var_name=name['var_name'],
                                               long_name=name['long_name'],
                                               units=name['units'])])

    return cubes


def get_provenance_record(ancestor_files, caption, statistics,
                          plot_type='scatter'):
    """Get Provenance record."""
    record = {
        'caption': caption,
        'statistics': statistics,
        'domains': ['global'],
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


def plot_slope_regression(cfg, data_dict):
    """Scatter plot of linear regression slope, some variables (fig2a)."""
    if not cfg[n.WRITE_PLOTS]:
        return

    data_model = data_dict['regressions']
    sa_lvpdt = data_model[:, 3]
    sa_rsnstdt = data_model[:, 1]
    sa_rsnstcsdt = data_model[:, 5]
    datasets = data_dict['datasets']

    m_all = np.array([np.mean(sa_lvpdt), np.mean(sa_rsnstdt),
                      np.mean(sa_rsnstcsdt)])

    # Regression between LvdP/dtas and the clr-dSWA/dtas and all-dSWA/dtas
    reg_sa = stats.linregress(sa_rsnstcsdt, sa_lvpdt)
    reg_sa_all = stats.linregress(sa_rsnstdt, sa_lvpdt)
    y_reg_sa = reg_sa.slope * np.linspace(0.2, 1.4, 2) + reg_sa.intercept
    text_sa = '{:.2f}'.format(reg_sa.rvalue)
    text_sa_all = '{:.2f}'.format(reg_sa_all.rvalue)

    # Regression between clr-dSWA/dtas and all-dSWA/dtas
    reg_rsnst = stats.linregress(sa_rsnstcsdt, sa_rsnstdt)
    y_reg_rsnst = reg_rsnst.slope * np.linspace(0.2, 1.4, 2) + \
        reg_rsnst.intercept

    # plt.style.use('/work/bd0854/b380216/esmvaltool/v2private/' +
    #               'ESMValTool-private/esmvaltool/diag_scripts/testkw/' +
    #               'style_kw_deangelis2.mplstyle')

    fig, axx = plt.subplots(figsize=(7, 7))
    axx.plot(np.arange(len(m_all)) + 1, m_all,
             linestyle='none', marker='x',
             markersize=25, markeredgewidth=4.0, markerfacecolor='r',
             markeredgecolor='r', label='Model mean')

    if not cfg[n.OUTPUT_FILE_TYPE] == 'eps':
        axx.plot(np.ones((len(sa_lvpdt))), sa_lvpdt,
                 linestyle='none', marker='o',
                 markersize=15, markeredgewidth=1.0, fillstyle='none',
                 markeredgecolor='b', label='Individual models')
        axx.plot(np.ones((len(sa_lvpdt))) + 1, sa_rsnstdt,
                 linestyle='none', marker='o',
                 markersize=15, markeredgewidth=1.0, fillstyle='none',
                 markeredgecolor='b')
        axx.plot(np.ones((len(sa_lvpdt))) + 2, sa_rsnstcsdt,
                 linestyle='none',
                 marker='o', markersize=15, markeredgewidth=1.0,
                 fillstyle='none', markeredgecolor='b')
    else:
        axx.plot(np.ones((len(sa_lvpdt))), sa_lvpdt,
                 linestyle='none', marker='o',
                 markersize=15, markeredgewidth=1.0,
                 markerfacecolor='w', markeredgecolor='b',
                 label='Individual models')
        axx.plot(np.ones((len(sa_lvpdt))) + 1, sa_rsnstdt,
                 linestyle='none', marker='o',
                 markersize=15, markeredgewidth=1.0,
                 markerfacecolor='w', markeredgecolor='b')
        axx.plot(np.ones((len(sa_lvpdt))) + 2, sa_rsnstcsdt,
                 linestyle='none',
                 marker='o', markersize=15, markeredgewidth=1.0,
                 markerfacecolor='w', markeredgecolor='b')

    axx.plot(np.arange(len(m_all)) + 1, m_all, linestyle='none', marker='x',
             markersize=25, markeredgewidth=4.0, markeredgecolor='r')

    axx.set_xlabel(' ')
    axx.set_title(' ')
    axx.set_ylabel(r'Temperature-mediated response (W m$^{-2}$ K$^{-1}$)')
    axx.set_xlim([0.5, 3.5])
    axx.set_xticks(np.linspace(1.0, 3.0, 3))
    axx.set_xticklabels(("dlvp/dtas", "drsnst/dtas", "rsnstcs/dtas"),
                        rotation=45, ha='right', rotation_mode='anchor')
    axx.set_ylim([0, 3.0])
    axx.set_yticks(np.linspace(0.5, 2.5, 5))
    axx.text(1.9, 0.2, text_sa)
    axx.text(2.9, 0.2, text_sa_all)
    axx.legend(loc=2)

    fig.tight_layout()
    fig.savefig(get_plot_filename('fig2a', cfg), dpi=300)
    plt.close()
    
    caption = 'The temperature-mediated response of each atmospheric ' + \
              'energy budget term for each model as blue circles and ' + \
              'the model mean as a red cross. The numbers above the ' + \
              'abscissa are the cross-model correlations between ' + \
              'dlvp/dtas and each other temperature-mediated response.'

    provenance_record = get_provenance_record(
        _get_sel_files_var(cfg, ['lvp', 'rsnst', 'rsnstcs', 'tas']),
                           caption, ['mean'])

    diagnostic_file = get_diagnostic_filename('fig2a', cfg)

    logger.info("Saving analysis results to %s", diagnostic_file)

    iris.save(cube_to_save_matrix(data_model, {'var_name': 'all',
                                               'long_name': 'dlvp/dtas, ' + \
                                                            'drsnst/dtas, ' + \
                                                            'drsnstcs/dtas',
                                               'units': 'W m-2 K-1'}),
              target=diagnostic_file)

    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                pformat(provenance_record))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)

    fig, axx = plt.subplots(figsize=(7, 7))

    axx.plot(np.linspace(0.2, 1.4, 2), y_reg_sa, color='r')

    for iii, model in enumerate(datasets):
        style = e.plot.get_dataset_style(model)
        axx.plot(
            sa_rsnstcsdt[iii],
            sa_lvpdt[iii],
            marker=style['mark'],
            color=style['color'],
            markerfacecolor=style['facecolor'],
            linestyle='none',
            markersize=10,
            markeredgewidth=2.0,
            label=model)

    axx.set_xlabel(r'drsnstcs/dtas (W m$^{-2}$ K$^{-1}$)')
    axx.set_title(' ')
    axx.set_ylabel(r'dlvp/dtas (W m$^{-2}$ K$^{-1}$)')
    axx.set_xlim([0.3, 1.35])
    axx.set_xticks(np.linspace(0.4, 1.2, 5))
    axx.set_ylim([1.75, 2.8])
    axx.set_yticks(np.linspace(1.8, 2.8, 6))
    axx.text(0.9, 2.75, 'Fit (r={:.2f}, '.format(reg_sa.rvalue) +
             ' slope = {:.2f}, '.format(reg_sa.slope) +
             ')')
    axx.legend(loc=3)

    fig.tight_layout()
    fig.savefig(get_plot_filename('fig2b', cfg), dpi=300)
    plt.close()

    fig, axx = plt.subplots(figsize=(7, 7))

    axx.plot(np.linspace(0.2, 1.4, 2), y_reg_rsnst, color='r')

    for iii, model in enumerate(datasets):
        style = e.plot.get_dataset_style(model)
        axx.plot(
            sa_rsnstcsdt[iii],
            sa_rsnstdt[iii],
            marker=style['mark'],
            color=style['color'],
            markerfacecolor=style['facecolor'],
            linestyle='none',
            markersize=10,
            markeredgewidth=2.0,
            label=model)

    axx.set_xlabel(r'drsnstcs/dtas (W m$^{-2}$ K$^{-1}$)')
    axx.set_title(' ')
    axx.set_ylabel(r'drsnst/dtas (W m$^{-2}$ K$^{-1}$)')
    axx.set_xlim([0.45, 1.15])
    axx.set_xticks(np.linspace(0.5, 1.1, 7))
    axx.set_ylim([0.45, 1.15])
    axx.set_yticks(np.linspace(0.5, 1.1, 7))
    axx.text(0.85, 1.1, 'Fit (r={:.2f}, '.format(reg_rsnst.rvalue) +
             ' slope = {:.2f}, '.format(reg_rsnst.slope) +
             ')')
    axx.legend(loc=2)

    fig.tight_layout()
    fig.savefig(get_plot_filename('exfig2b', cfg), dpi=300)
    plt.close()


def plot_slope_regression_all(cfg, data_dict, available_vars, available_exp):
    """Scatter plot of linear regression slope, all variables (exfig2a)."""
    if not cfg[n.WRITE_PLOTS]:
        return

    data_model = data_dict['regressions']
    m_all = np.array([np.mean(data_model[:, 3]), np.mean(data_model[:, 0]),
                      np.mean(data_model[:, 1]), np.mean(data_model[:, 2]),
                      np.mean(data_model[:, 4]), np.mean(data_model[:, 5])])

    reg_rsnstcsdt = stats.linregress(data_model[:, 5], data_model[:, 3])
    reg_rsnstdt = stats.linregress(data_model[:, 1], data_model[:, 3])
    reg_rlnstcsdt = stats.linregress(data_model[:, 4], data_model[:, 3])
    reg_rlnstdt = stats.linregress(data_model[:, 0], data_model[:, 3])
    reg_hfssdt = stats.linregress(data_model[:, 2], data_model[:, 3])

    text_rsnstcsdt = '{:.2f}'.format(reg_rsnstcsdt.rvalue)
    text_rsnstdt = '{:.2f}'.format(reg_rsnstdt.rvalue)
    text_rlnstcsdt = '{:.2f}'.format(reg_rlnstcsdt.rvalue)
    text_rlnstdt = '{:.2f}'.format(reg_rlnstdt.rvalue)
    text_hfssdt = '{:.2f}'.format(reg_hfssdt.rvalue)

    fig, axx = plt.subplots(figsize=(7, 7))
    axx.plot(np.arange(len(m_all)) + 1, m_all,
             linestyle='none', marker='x',
             markersize=25, markeredgewidth=4.0, markerfacecolor='r',
             markeredgecolor='r', label='Model mean')
    if not cfg[n.OUTPUT_FILE_TYPE] == 'eps':
        axx.plot(np.ones((len(data_model[:, 2]))), data_model[:, 3],
                 linestyle='none', marker='o',
                 markersize=15, markeredgewidth=1.0, fillstyle='none',
                 markeredgecolor='b', label='Individual models')
        axx.plot(np.ones((len(data_model[:, 2]))) + 1, data_model[:, 0],
                 linestyle='none', marker='o',
                 markersize=15, markeredgewidth=1.0, fillstyle='none',
                 markeredgecolor='b')
        axx.plot(np.ones((len(data_model[:, 2]))) + 2, data_model[:, 1],
                 linestyle='none', marker='o',
                 markersize=15, markeredgewidth=1.0, fillstyle='none',
                 markeredgecolor='b')
        axx.plot(np.ones((len(data_model[:, 2]))) + 3, data_model[:, 2],
                 linestyle='none', marker='o',
                 markersize=15, markeredgewidth=1.0, fillstyle='none',
                 markeredgecolor='b')
        axx.plot(np.ones((len(data_model[:, 2]))) + 4, data_model[:, 4],
                 linestyle='none', marker='o',
                 markersize=15, markeredgewidth=1.0, fillstyle='none',
                 markeredgecolor='b')
        axx.plot(np.ones((len(data_model[:, 2]))) + 5, data_model[:, 5],
                 linestyle='none', marker='o',
                 markersize=15, markeredgewidth=1.0, fillstyle='none',
                 markeredgecolor='b')
    else:
        axx.plot(np.ones((len(data_model[:, 2]))), data_model[:, 2],
                 linestyle='none', marker='o',
                 markersize=15, markeredgewidth=1.0, markerfacecolor='w',
                 markeredgecolor='b', label='Individual models')
        axx.plot(np.ones((len(data_model[:, 2]))) + 1, data_model[:, 0],
                 linestyle='none', marker='o',
                 markersize=15, markeredgewidth=1.0, markerfacecolor='w',
                 markeredgecolor='b')
        axx.plot(np.ones((len(data_model[:, 2]))) + 2, data_model[:, 1],
                 linestyle='none', marker='o',
                 markersize=15, markeredgewidth=1.0, markerfacecolor='w',
                 markeredgecolor='b')
        axx.plot(np.ones((len(data_model[:, 2]))) + 3, data_model[:, 2],
                 linestyle='none', marker='o',
                 markersize=15, markeredgewidth=1.0, markerfacecolor='w',
                 markeredgecolor='b')
        axx.plot(np.ones((len(data_model[:, 2]))) + 4, data_model[:, 4],
                 linestyle='none', marker='o',
                 markersize=15, markeredgewidth=1.0, markerfacecolor='w',
                 markeredgecolor='b')
        axx.plot(np.ones((len(data_model[:, 2]))) + 5, data_model[:, 5],
                 acecolor='w',
                 markeredgecolor='b')

    axx.plot(np.arange(len(m_all)) + 1, m_all, linestyle='none', marker='x',
             markersize=25, markeredgewidth=4.0, markerfacecolor='r',
             markeredgecolor='r')

    axx.set_xlabel(' ')
    axx.set_title(' ')
    axx.set_ylabel(r'Temperature-mediated response (W m$^{-2}$ K$^{-1}$)')
    axx.set_xlim([0.5, 6.5])
    axx.set_xticks(np.linspace(1.0, 6.0, 6))
    axx.set_xticklabels(("dlvp/dtas", "drlnst/dtas", "drsnst/dtas",
                         "dhfss/dtas", "drlnstcs/dtas", "drsnstcs/dtas"),
                        rotation=45, ha='right', rotation_mode='anchor')
    axx.set_ylim([-1.5, 4.5])
    axx.set_yticks(np.linspace(-1.0, 4.0, 11))
    axx.vlines([1.5, 4.5], -2, 5, colors='k', linestyle='solid')
    axx.hlines(0, 0, 7, colors='k', linestyle='dashed')
    axx.text(1.9, 0.2, text_rlnstdt)
    axx.text(2.9, 0.2, text_rsnstdt)
    axx.text(3.9, 0.2, text_hfssdt)
    axx.text(4.9, 0.2, text_rlnstcsdt)
    axx.text(5.9, 0.2, text_rsnstcsdt)
    axx.legend(loc=2)

    fig.tight_layout()
    fig.savefig(get_plot_filename('exfig2a', cfg), dpi=300)
    plt.close()

    caption = 'The temperature-mediated response of each atmospheric ' + \
              'energy budget term for each model as blue circles and ' + \
              'the model mean as a red cross. The numbers above the ' + \
              'abscissa are the cross-model correlations between ' + \
              'dlvp/dtas and each other temperature-mediated response.'

    provenance_record = get_provenance_record(
        _get_sel_files_var(cfg, available_vars), caption, ['mean'])

    diagnostic_file = get_diagnostic_filename('exfig2a', cfg)

    logger.info("Saving analysis results to %s", diagnostic_file)

    iris.save(cube_to_save_matrix(data_model, {'var_name': 'all',
                                               'long_name': 'dlvp/dtas, ' + \
                                                            'drlnst/dtas, ' + \
                                                            'drsnst/dtas, ' + \
                                                            'dhfss/dtas, ' + \
                                                            'drlnstcs/' + \
                                                            'dtas, ' + \
                                                            'drsnstcs/dtas',
                                               'units': 'W m-2 K-1'}),
              target=diagnostic_file)

    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                pformat(provenance_record))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)


def plot_rlnst_regression(cfg, dataset_name, data, variables, regs):
    """Plot linear regression used to calculate ECS."""
    if not cfg[n.WRITE_PLOTS]:
        return

    filepath = get_plot_filename(dataset_name, cfg)

    # Regression line
    # x_reg = np.linspace(0.0, 7.0, 2)
    y_reg = regs["rlnst"].slope * np.linspace(0.0, 7.0, 2) + \
        regs["rlnst"].intercept
    y_reg2 = regs["rsnst"].slope * np.linspace(0.0, 7.0, 2) + \
        regs["rsnst"].intercept
    y_reg3 = regs["hfss"].slope * np.linspace(0.0, 7.0, 2) + \
        regs["hfss"].intercept
    y_reg4 = regs["lvp"].slope * np.linspace(0.0, 7.0, 2) + \
        regs["lvp"].intercept

    lab4 = 'dlvp/dtas = {:.2f}, '.format(regs["lvp"].slope) + \
        'y−int  = {:.2f}, '.format(regs["lvp"].intercept) + \
        'r = {:.2f}'.format(regs["lvp"].rvalue)
    lab = 'drlnst/dtas = {:.2f}, '.format(regs["rlnst"].slope) + \
        'y−int  = {:.2f}, '.format(regs["rlnst"].intercept) + \
        'r = {:.2f}'.format(regs["rlnst"].rvalue)
    lab2 = 'drsnst/dtas = {:.2f}, '.format(regs["rsnst"].slope) + \
        'y−int  = {:.2f}, '.format(regs["rsnst"].intercept) + \
        'r = {:.2f}'.format(regs["rsnst"].rvalue)
    lab3 = 'dhfss/dtas = {:.2f}, '.format(regs["hfss"].slope) + \
        'y−int  = {:.2f}, '.format(regs["hfss"].intercept) + \
        'r = {:.2f}'.format(regs["hfss"].rvalue)

    axhline_dict = {'y': 0, 'linestyle': 'dashed', 'color': 'k',
                    'linewidth': 2.0}

    e.plot.scatterplot(
        [data[0], np.linspace(0.0, 7.0, 2), data[0], np.linspace(0.0, 7.0, 2),
         data[0], np.linspace(0.0, 7.0, 2), data[0], np.linspace(0.0, 7.0, 2)],
        [data[4], y_reg4, data[1], y_reg, data[2], y_reg2, data[3], y_reg3],
        filepath,
        plot_kwargs=[{'linestyle': 'none',
                      'marker': 'o',
                      'markerfacecolor': 'g',
                      'markeredgecolor': 'g',
                      'label': lab4},
                     {'color': 'g',
                      'linestyle': '-'},
                     {'linestyle': 'none',
                      'marker': '^',
                      'markerfacecolor': 'b',
                      'markeredgecolor': 'b',
                      'label': lab},
                     {'color': 'b',
                      'linestyle': '-'},
                     {'linestyle': 'none',
                      'marker': 's',
                      'markerfacecolor': 'r',
                      'markeredgecolor': 'r',
                      'label': lab2},
                     {'color': 'r',
                      'linestyle': '-'},
                     {'linestyle': 'none',
                      'marker': '*',
                      'markerfacecolor': 'tab:gray',
                      'markeredgecolor': 'tab:gray',
                      'label': lab3},
                     {'color': 'tab:gray',
                      'linestyle': '-'}],
        save_kwargs={'bbox_inches': 'tight',
                     'orientation': 'landscape'},
        axes_functions={'set_title': dataset_name,
                        'set_xlabel': '2−m temperature (T)' +
                                      'global−mean annual anomaly (' +
                                      variables.units('tas') + ')',
                        'set_ylabel': r'Energy budget term global - ' +
                                      'mean annual anomalies (W m$^{-2}$)',
                        'set_xlim': [0, 7.0],
                        'set_ylim': [-5.0, 17.0],
                        'set_yticks': np.linspace(-4, 16, 11),
                        'axhline': axhline_dict,
                        'legend': {'loc': 2}})


def substract_and_reg_deangelis2(cfg, data, var):
    """Substract piControl from abrupt4xCO2 for all models and variables."""
    pathlist = data.get_path_list(short_name='tas', exp=PICONTROL)
    regressions = np.zeros((len(pathlist), 6))
    datasets = []

    data_var = OrderedDict()
    reg_var = OrderedDict()
    varvar = var.short_names()

    for iii, dataset_path in enumerate(pathlist):

        # Substract piControl experiment from abrupt4xCO2 experiment
        dataset = data.get_info(n.DATASET, dataset_path)
        datasets.append(dataset)

        for jvar in varvar:
            data_var[jvar] = data.get_data(short_name=jvar, exp=ABRUPT4XCO2,
                                           dataset=dataset) - \
                data.get_data(short_name=jvar, exp=PICONTROL,
                              dataset=dataset)

        # Perform linear regression
        for jvar in varvar:
            if jvar != 'tas':
                reg_var[jvar] = stats.linregress(data_var["tas"],
                                                 data_var[jvar])

        # Plot ECS regression if desired
        plot_rlnst_regression(cfg, dataset, [data_var["tas"],
                                             data_var["rlnst"],
                                             data_var["rsnst"],
                                             data_var["hfss"],
                                             data_var["lvp"]],
                              var, reg_var)

        # Save data
        regressions[iii] = [reg_var["rlnst"].slope, reg_var["rsnst"].slope,
                            reg_var["hfss"].slope, reg_var["lvp"].slope,
                            reg_var["rlnstcs"].slope,
                            reg_var["rsnstcs"].slope]

    return dict([('regressions', regressions), ('datasets', datasets)])


###############################################################################
# Setup diagnostic
###############################################################################

# Variables
# ECS = e.Variable('ecs',

# Experiments
PICONTROL = 'piControl'
ABRUPT4XCO2 = 'abrupt4xCO2'


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
    # logging.debug("Found variables in recipe:\n%s", var)

    available_vars = list(group_metadata(cfg['input_data'].values(),
                                         'short_name'))
    logging.debug("Found variables in recipe:\n%s", available_vars)

    available_exp = list(group_metadata(cfg['input_data'].values(), 'exp'))

    # Check for available variables
    required_vars = ('tas', 'lvp', 'rlnst', 'rsnst', 'rlnstcs',
                     'rsnstcs', 'hfss')
    if not e.variables_available(cfg, required_vars):
        raise ValueError("This diagnostic needs {required_vars} variables")

    # Check for experiments
    if 'abrupt-4xCO2' not in available_exp:
        if 'abrupt4xCO2' not in available_exp:
            raise ValueError("The diagnostic needs an experiment with " +
                             "4 times CO2.")

    if 'piControl' not in available_exp:
        raise ValueError("The diagnostic needs a pre industrial control " +
                         "experiment.")

    ###########################################################################
    # Read data
    ###########################################################################

    # Create iris cube for each dataset and save annual means
    for dataset_path in data:
        cube = iris.load(dataset_path)[0]
        cat.add_year(cube, 'time', name='year')
        cube = cube.aggregated_by('year', iris.analysis.MEAN)
        experiment = data.get_info(n.EXP, dataset_path)
        if experiment == PICONTROL:
            # DeAngelis use a 21 year running mean on piControl but the
            # full extend of 150 years abrupt4xCO2. I could not find out,
            # how they tread the edges, currently I just skip the mean for
            # the edges. This is not exacly the same as done in the paper,
            # small differences remain in extended data Fig 1,
            # but closer than other methods I
            # tried, e.g. skipping the edges.
            # For most data sets it would be also possible to
            # extend the piControl for 20 years, but then it would
            # not be centered means of piControl for each year of
            # abrupt4xCO2 any more.
            cube_new = cube.rolling_window('time', iris.analysis.MEAN, 21)
            endm10 = len(cube.coord('time').points) - 10
            cube.data[10:endm10] = cube_new.data

        data.set_data(cube.data, dataset_path)

    ###########################################################################
    # Process data
    ###########################################################################

    data_dict = substract_and_reg_deangelis2(cfg, data, var)

    plot_slope_regression(cfg, data_dict)
    plot_slope_regression_all(cfg, data_dict, available_vars, available_exp)


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
