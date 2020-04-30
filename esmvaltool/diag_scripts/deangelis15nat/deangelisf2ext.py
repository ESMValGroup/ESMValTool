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
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            get_plot_filename, group_metadata,
                                            select_metadata)

logger = logging.getLogger(os.path.basename(__file__))


def _set_list_dict1(sa_dict):
    list_dict = {}
    list_dict["data"] = [sa_dict["rsnstcsdt"], sa_dict["rsnstdt"]]
    list_dict["name"] = [{'var_name': 'drsnstcs_divby_dtas',
                          'long_name': 'Temperature mediated ' +
                                       'shortwave absorption for clear skye',
                          'units': 'W m-2 K-1'},
                         {'var_name': 'drsnst_divby_dtas',
                          'long_name': 'Temperature mediated ' +
                                       'shortwave absorption for all skye',
                          'units': 'W m-2 K-1'}]
    return list_dict


def _set_list_dict2(sa_dict):
    list_dict = {}
    list_dict["data"] = [sa_dict["lvpdt"], sa_dict["rsnstcsdt"]]
    list_dict["name"] = [{'var_name': 'dlvp_divby_dtas',
                          'long_name': 'Temperature mediated latent heat ' +
                                       'release from precipitation',
                          'units': 'W m-2 K-1'},
                         {'var_name': 'drsnstcs_divby_dtas',
                          'long_name': 'Temperature mediated ' +
                                       'shortwave absorption for clear skye',
                          'units': 'W m-2 K-1'}]

    return list_dict


def _calculate_regression_sa(sa_dict):
    """Regression between dlvp/dtas, drsnstcs/dtas, drsnst/dtas."""
    # Regression between LvdP/dtas and the clr-dSWA/dtas and all-dSWA/dtas
    reg_dict = {}
    reg_dict["sa"] = stats.linregress(sa_dict["rsnstcsdt"], sa_dict["lvpdt"])
    reg_dict["sa_all"] = stats.linregress(sa_dict["rsnstdt"], sa_dict["lvpdt"])
    reg_dict["y_sa"] = reg_dict["sa"].slope * np.linspace(0.2, 1.4, 2) + \
        reg_dict["sa"].intercept

    # Regression between clr-dSWA/dtas and all-dSWA/dtas
    reg_dict["rsnst"] = stats.linregress(sa_dict["rsnstcsdt"],
                                         sa_dict["rsnstdt"])
    reg_dict["y_rsnst"] = reg_dict["rsnst"].slope * \
        np.linspace(0.2, 1.4, 2) + reg_dict["rsnst"].intercept

    return reg_dict


def _set_axx_fig2a(cfg, axx, m_all, reg_dict, sa_dict):
    """Text for fig2a."""
    text_sa = '{:.2f}'.format(reg_dict["sa"].rvalue)
    text_sa_all = '{:.2f}'.format(reg_dict["sa_all"].rvalue)

    axx.plot(np.arange(len(m_all)) + 1, m_all,
             linestyle='none', marker='x',
             markersize=25, markeredgewidth=4.0, markerfacecolor='r',
             markeredgecolor='r', label='Model mean')
    axx.plot(np.arange(len(m_all)) + 1, m_all, linestyle='none', marker='x',
             markersize=25, markeredgewidth=4.0, markeredgecolor='r')

    if not cfg[n.OUTPUT_FILE_TYPE] == 'eps':
        axx.plot(np.ones((len(sa_dict["lvpdt"]))), sa_dict["lvpdt"],
                 linestyle='none', marker='o',
                 markersize=15, markeredgewidth=1.0, fillstyle='none',
                 markeredgecolor='b', label='Individual models')
        axx.plot(np.ones((len(sa_dict["lvpdt"]))) + 1, sa_dict["rsnstdt"],
                 linestyle='none', marker='o',
                 markersize=15, markeredgewidth=1.0, fillstyle='none',
                 markeredgecolor='b')
        axx.plot(np.ones((len(sa_dict["lvpdt"]))) + 2, sa_dict["rsnstcsdt"],
                 linestyle='none',
                 marker='o', markersize=15, markeredgewidth=1.0,
                 fillstyle='none', markeredgecolor='b')
    else:
        axx.plot(np.ones((len(sa_dict["lvpdt"]))), sa_dict["lvpdt"],
                 linestyle='none', marker='o',
                 markersize=15, markeredgewidth=1.0,
                 markerfacecolor='w', markeredgecolor='b',
                 label='Individual models')
        axx.plot(np.ones((len(sa_dict["lvpdt"]))) + 1, sa_dict["rsnstdt"],
                 linestyle='none', marker='o',
                 markersize=15, markeredgewidth=1.0,
                 markerfacecolor='w', markeredgecolor='b')
        axx.plot(np.ones((len(sa_dict["lvpdt"]))) + 2, sa_dict["rsnstcsdt"],
                 linestyle='none',
                 marker='o', markersize=15, markeredgewidth=1.0,
                 markerfacecolor='w', markeredgecolor='b')

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

    return axx


def _set_axx_fig2b(axx, cfg, reg_dict, datasets, sa_dict):
    """Text for fig2b."""
    axx.plot(np.linspace(0.2, 1.4, 2), reg_dict["y_sa"], color='r')

    for iii, model in enumerate(datasets):
        proj = (select_metadata(cfg['input_data'].values(),
                                dataset=model))[0]['project']
        style = e.plot.get_dataset_style(model, style_file=proj.lower())
        axx.plot(
            sa_dict["rsnstcsdt"][iii],
            sa_dict["lvpdt"][iii],
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
    axx.text(0.9, 2.75, 'Fit (r={:.2f}, '.format(reg_dict["sa"].rvalue) +
             ' slope = {:.2f}, '.format(reg_dict["sa"].slope) +
             ')')
    axx.legend(loc=3)
    return axx


def _set_text_exfig2a(axx, text_dict):
    """Text for exfig2a."""
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
    axx.text(1.9, 0.2, text_dict["rlnstdt"])
    axx.text(2.9, 0.2, text_dict["rsnstdt"])
    axx.text(3.9, 0.2, text_dict["hfssdt"])
    axx.text(4.9, 0.2, text_dict["rlnstcsdt"])
    axx.text(5.9, 0.2, text_dict["rsnstcsdt"])
    axx.legend(loc=2)

    return axx


def _set_axx_exfig2b(axx, cfg, datasets, reg_dict, sa_dict):
    """Text for exfig2b."""
    axx.plot(np.linspace(0.2, 1.4, 2), reg_dict["y_rsnst"], color='r')

    for iii, model in enumerate(datasets):
        proj = (select_metadata(cfg['input_data'].values(),
                                dataset=model))[0]['project']
        style = e.plot.get_dataset_style(model, style_file=proj.lower())
        axx.plot(
            sa_dict["rsnstcsdt"][iii],
            sa_dict["rsnstdt"][iii],
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
    axx.text(0.85, 1.1, 'Fit (r={:.2f}, '.format(reg_dict["rsnst"].rvalue) +
             ' slope = {:.2f}, '.format(reg_dict["rsnst"].slope) +
             ')')
    axx.legend(loc=2)

    return axx


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

    sa_dict = {}
    sa_dict["lvpdt"] = data_dict['regressions'][:, 3]
    sa_dict["rsnstdt"] = data_dict['regressions'][:, 1]
    sa_dict["rsnstcsdt"] = data_dict['regressions'][:, 5]
    datasets = data_dict['datasets']

    m_all = np.array([np.mean(sa_dict["lvpdt"]), np.mean(sa_dict["rsnstdt"]),
                      np.mean(sa_dict["rsnstcsdt"])])

    reg_dict = _calculate_regression_sa(sa_dict)

    # plt.style.use('/work/bd0854/b380216/esmvaltool/v2private/' +
    #               'ESMValTool-private/esmvaltool/diag_scripts/testkw/' +
    #               'style_kw_deangelis2.mplstyle')

    fig, axx = plt.subplots(figsize=(7, 7))

    axx = _set_axx_fig2a(cfg, axx, m_all, reg_dict, sa_dict)

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
        caption, ['corr', 'mean'])

    diagnostic_file = get_diagnostic_filename('fig2a', cfg)

    logger.info("Saving analysis results to %s", diagnostic_file)

    list_dict = {}
    list_dict["data"] = [sa_dict["lvpdt"], sa_dict["rsnstdt"],
                         sa_dict["rsnstcsdt"]]
    list_dict["name"] = [{'var_name': 'dlvp_divby_dtas',
                          'long_name': 'Temperature mediated latent heat ' +
                                       'release from precipitation',
                          'units': 'W m-2 K-1'},
                         {'var_name': 'drsnst_divby_dtas',
                          'long_name': 'Temperature mediated ' +
                                       'shortwave absorption',
                          'units': 'W m-2 K-1'},
                         {'var_name': 'drsnstcs_divby_dtas',
                          'long_name': 'Temperature mediated ' +
                                       'shortwave absorption for clear skye',
                          'units': 'W m-2 K-1'}]

    iris.save(cube_to_save_vars(list_dict), target=diagnostic_file)

    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                pformat(provenance_record))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)

    fig, axx = plt.subplots(figsize=(7, 7))

    axx = _set_axx_fig2b(axx, cfg, reg_dict, datasets, sa_dict)

    fig.tight_layout()
    fig.savefig(get_plot_filename('fig2b', cfg), dpi=300)
    plt.close()

    caption = 'Scatterplot of dlvp/dtas versus drsnstcs/dtas with ' + \
              'corresponding least-squares linear fit (red line).'

    provenance_record = get_provenance_record(
        _get_sel_files_var(cfg, ['lvp', 'rsnstcs', 'tas']),
        caption, ['corr'])

    diagnostic_file = get_diagnostic_filename('fig2b', cfg)

    logger.info("Saving analysis results to %s", diagnostic_file)

    list_dict = _set_list_dict2(sa_dict)

    iris.save(cube_to_save_vars(list_dict), target=diagnostic_file)

    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                pformat(provenance_record))

    fig, axx = plt.subplots(figsize=(7, 7))

    axx = _set_axx_exfig2b(axx, cfg, datasets, reg_dict, sa_dict)

    fig.tight_layout()
    fig.savefig(get_plot_filename('exfig2b', cfg), dpi=300)
    plt.close()

    caption = 'Scatterplot of drsnstcs/dtas versus drsnst/dtas with ' + \
              'corresponding least-squares linear fit (red line).'

    provenance_record = get_provenance_record(
        _get_sel_files_var(cfg, ['rsnstcs', 'rsnst', 'tas']),
        caption, ['corr'])

    diagnostic_file = get_diagnostic_filename('exfig2b', cfg)

    logger.info("Saving analysis results to %s", diagnostic_file)

    list_dict = _set_list_dict1(sa_dict)

    iris.save(cube_to_save_vars(list_dict), target=diagnostic_file)

    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                pformat(provenance_record))


def plot_slope_regression_all(cfg, data_dict, available_vars):
    """Scatter plot of linear regression slope, all variables (exfig2a)."""
    if not cfg[n.WRITE_PLOTS]:
        return

    data_model = data_dict['regressions']
    m_all = np.array([np.mean(data_model[:, 3]), np.mean(data_model[:, 0]),
                      np.mean(data_model[:, 1]), np.mean(data_model[:, 2]),
                      np.mean(data_model[:, 4]), np.mean(data_model[:, 5])])

    reg_dict = {}
    reg_dict["rsnstcsdt"] = stats.linregress(data_model[:, 5],
                                             data_model[:, 3])
    reg_dict["rsnstdt"] = stats.linregress(data_model[:, 1], data_model[:, 3])
    reg_dict["rlnstcsdt"] = stats.linregress(data_model[:, 4],
                                             data_model[:, 3])
    reg_dict["rlnstdt"] = stats.linregress(data_model[:, 0], data_model[:, 3])
    reg_dict["hfssdt"] = stats.linregress(data_model[:, 2], data_model[:, 3])

    text_dict = {}
    text_dict["rsnstcsdt"] = '{:.2f}'.format(reg_dict["rsnstcsdt"].rvalue)
    text_dict["rsnstdt"] = '{:.2f}'.format(reg_dict["rsnstdt"].rvalue)
    text_dict["rlnstcsdt"] = '{:.2f}'.format(reg_dict["rlnstcsdt"].rvalue)
    text_dict["rlnstdt"] = '{:.2f}'.format(reg_dict["rlnstdt"].rvalue)
    text_dict["hfssdt"] = '{:.2f}'.format(reg_dict["hfssdt"].rvalue)

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

    axx = _set_text_exfig2a(axx, text_dict)

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
                                               'long_name': 'dlvp, ' +
                                                            'drlnst, ' +
                                                            'drsnst, ' +
                                                            'dhfss, ' +
                                                            'drlnstcs, and,' +
                                                            'drsnstcs ' +
                                                            'divided by ' +
                                                            'dtas',
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
    yreg_dict = {}
    yreg_dict["rlnst"] = regs["rlnst"].slope * np.linspace(0.0, 7.0, 2) + \
        regs["rlnst"].intercept
    yreg_dict["rsnst"] = regs["rsnst"].slope * np.linspace(0.0, 7.0, 2) + \
        regs["rsnst"].intercept
    yreg_dict["hfss"] = regs["hfss"].slope * np.linspace(0.0, 7.0, 2) + \
        regs["hfss"].intercept
    yreg_dict["lvp"] = regs["lvp"].slope * np.linspace(0.0, 7.0, 2) + \
        regs["lvp"].intercept

    yreg_dict["lab_lvp"] = 'dlvp/dtas = {:.2f}, '.format(regs["lvp"].slope) + \
        'y−int  = {:.2f}, '.format(regs["lvp"].intercept) + \
        'r = {:.2f}'.format(regs["lvp"].rvalue)
    yreg_dict["lab_rlnst"] = 'drlnst' +  \
        '/dtas = {:.2f}, '.format(regs["rlnst"].slope) + \
        'y−int  = {:.2f}, '.format(regs["rlnst"].intercept) + \
        'r = {:.2f}'.format(regs["rlnst"].rvalue)
    yreg_dict["lab_rsnst"] = 'drsnst' + \
        '/dtas = {:.2f}, '.format(regs["rsnst"].slope) + \
        'y−int  = {:.2f}, '.format(regs["rsnst"].intercept) + \
        'r = {:.2f}'.format(regs["rsnst"].rvalue)
    yreg_dict["lab_hfss"] = 'dhfss' + \
        '/dtas = {:.2f}, '.format(regs["hfss"].slope) + \
        'y−int  = {:.2f}, '.format(regs["hfss"].intercept) + \
        'r = {:.2f}'.format(regs["hfss"].rvalue)

    axhline_dict = {'y': 0, 'linestyle': 'dashed', 'color': 'k',
                    'linewidth': 2.0}

    e.plot.scatterplot(
        [data["tas"], np.linspace(0.0, 7.0, 2), data["tas"],
         np.linspace(0.0, 7.0, 2),
         data["tas"], np.linspace(0.0, 7.0, 2), data["tas"],
         np.linspace(0.0, 7.0, 2)],
        [data["lvp"], yreg_dict["lvp"], data["rlnst"], yreg_dict["rlnst"],
         data["rsnst"], yreg_dict["rsnst"], data["hfss"], yreg_dict["hfss"]],
        filepath,
        plot_kwargs=[{'linestyle': 'none',
                      'marker': 'o',
                      'markerfacecolor': 'g',
                      'markeredgecolor': 'g',
                      'label': yreg_dict["lab_lvp"]},
                     {'color': 'g',
                      'linestyle': '-'},
                     {'linestyle': 'none',
                      'marker': '^',
                      'markerfacecolor': 'b',
                      'markeredgecolor': 'b',
                      'label': yreg_dict["lab_rsnst"]},
                     {'color': 'b',
                      'linestyle': '-'},
                     {'linestyle': 'none',
                      'marker': 's',
                      'markerfacecolor': 'r',
                      'markeredgecolor': 'r',
                      'label': yreg_dict["lab_rlnst"]},
                     {'color': 'r',
                      'linestyle': '-'},
                     {'linestyle': 'none',
                      'marker': '*',
                      'markerfacecolor': 'tab:gray',
                      'markeredgecolor': 'tab:gray',
                      'label': yreg_dict["lab_hfss"]},
                     {'color': 'tab:gray',
                      'linestyle': '-'}],
        save_kwargs={'bbox_inches': 'tight',
                     'orientation': 'landscape'},
        axes_functions={'set_title': dataset_name,
                        'set_xlabel': '2−m temperature (tas)' +
                                      'global−mean annual anomaly (' +
                                      variables.units('tas') + ')',
                        'set_ylabel': r'Energy budget term global - ' +
                                      'mean annual anomalies (W m$^{-2}$)',
                        'set_xlim': [0, 7.0],
                        'set_ylim': [-5.0, 17.0],
                        'set_yticks': np.linspace(-4, 16, 11),
                        'axhline': axhline_dict,
                        'legend': {'loc': 2}})

    caption = ' Demonstration of the Gregory method for ' + dataset_name + \
              '. Global-mean annual anomalies (' + ABRUPT4XCO2 + '—' + \
              PICONTROL + ' in atmospheric energy budget terms ' + \
              '(latent heat release from precipitation (lvp), ' + \
              'net longwave cooling (rlnst), shortwave absorption ' + \
              '(rsnst), and sensible heating (hfss)) are regressed ' + \
              'against those in 2-m air temperature. For lvp, ' + \
              'precipitation anomalies are multiplied by the latent ' + \
              'heat of vaporization. Radiative terms are computed with ' + \
              'all-sky fluxes. The statistics of the linear regression ' + \
              '(slope, y-intercept, and correlation coefficient, r) ' + \
              'are displayed in the key.'

    provenance_record = get_provenance_record(
        _get_sel_files_var(cfg, ['lvp', 'rlnst', 'rsnst', 'hfss', 'tas']),
        caption, ['corr'])

    diagnostic_file = get_diagnostic_filename(dataset_name, cfg)

    logger.info("Saving analysis results to %s", diagnostic_file)

    list_dict = {}
    list_dict["data"] = [data["tas"], data["lvp"], data["rlnst"],
                         data["rsnst"], data["hfss"]]
    list_dict["name"] = [{'var_name': 'tas',
                          'long_name': '2-m air temperature',
                          'units': 'K'},
                         {'var_name': 'lvp',
                          'long_name': 'Latent heat release ' +
                                       'from precipitation',
                          'units': 'W m-2'},
                         {'var_name': 'rlnst',
                          'long_name': 'Net longwave cooling',
                          'units': 'W m-2'},
                         {'var_name': 'rsnst',
                          'long_name': 'Shortwave absorption',
                          'units': 'W m-2'},
                         {'var_name': 'hfss',
                          'long_name': 'Sensible heating',
                          'units': 'W m-2'}]

    iris.save(cube_to_save_vars(list_dict), target=diagnostic_file)

    logger.info("Recording provenance of %s:\n%s", diagnostic_file,
                pformat(provenance_record))
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(diagnostic_file, provenance_record)


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
        plot_rlnst_regression(cfg, dataset, data_var, var, reg_var)

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
    plot_slope_regression_all(cfg, data_dict, available_vars)


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
