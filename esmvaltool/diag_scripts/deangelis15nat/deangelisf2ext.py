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
import iris
import iris.coord_categorisation as cat
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import esmvaltool.diag_scripts.shared as e
import esmvaltool.diag_scripts.shared.names as n

logger = logging.getLogger(os.path.basename(__file__))


def plot_slope_regression(cfg, data_model):
    """Scatter plot of linear regression slope, some variables (fig2a)."""
    # def plot_slope_regression(cfg, dlvpdt, drsnstdt, drsnstcsdt):

    if not (cfg[n.WRITE_PLOTS] and cfg.get('plot_ecs_regression')):
        return

    sa_lvpdt = data_model[:, 3]
    sa_rsnstdt = data_model[:, 1]
    sa_rsnstcsdt = data_model[:, 5]

    m_all = np.array([np.mean(sa_lvpdt), np.mean(sa_rsnstdt),
                      np.mean(sa_rsnstcsdt)])

    # Regression between LvdP/dT and the clr-dSWA/dT and all-dSWA/dT
    reg_sa = stats.linregress(sa_rsnstcsdt, sa_lvpdt)
    reg_sa_all = stats.linregress(sa_rsnstdt, sa_lvpdt)
    y_reg_sa = reg_sa.slope * np.linspace(0.2, 1.4, 2) + reg_sa.intercept
    text_sa = '{:.2f}'.format(reg_sa.rvalue)
    text_sa_all = '{:.2f}'.format(reg_sa_all.rvalue)

    # Regression between clr-dSWA/dT and all-dSWA/dT
    reg_rsnst = stats.linregress(sa_rsnstcsdt, sa_rsnstdt)
    y_reg_rsnst = reg_rsnst.slope * np.linspace(0.2, 1.4, 2) + \
        reg_rsnst.intercept

    plt.style.use('/work/bd0854/b380216/esmvaltool/v2private/' +
                  'ESMValTool-private/esmvaltool/diag_scripts/testkw/' +
                  'style_kw_deangelis2.mplstyle')

    fig, axx = plt.subplots()
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
    axx.set_xticklabels((r"L$_v$dP/dT", "all-dSWA/dT", "clr-dSWA/dT"),
                        rotation=45, ha='right', rotation_mode='anchor')
    axx.set_ylim([0, 3.0])
    axx.set_yticks(np.linspace(0.5, 2.5, 5))
    axx.text(1.9, 0.2, text_sa)
    axx.text(2.9, 0.2, text_sa_all)
    axx.legend(loc=2)

    fig.tight_layout()
    fig.savefig(os.path.join(cfg[n.PLOT_DIR], 'fig2.' +
                             cfg[n.OUTPUT_FILE_TYPE]))
    plt.close()

    # filepath2 = os.path.join(cfg[n.PLOT_DIR], 'fig2b.' +
    #                          cfg[n.OUTPUT_FILE_TYPE])

    e.plot.scatterplot(
        [sa_rsnstcsdt, np.linspace(0.2, 1.4, 2)],
        [sa_lvpdt, y_reg_sa],
        os.path.join(cfg[n.PLOT_DIR], 'fig2b.' +
                     cfg[n.OUTPUT_FILE_TYPE]),
        mpl_style_file='default.mplstyle',
        plot_kwargs=[{'linestyle': 'none',
                      'marker': '*',
                      'markersize': 15,
                      'markeredgewidth': 2.0,
                      'markerfacecolor': 'k',
                      'markeredgecolor': 'k'},
                     {'color': 'r',
                      'linestyle': '-',
                      'label': 'Fit (r={:.2f}, '.format(reg_sa.rvalue) +
                               ' slope = {:.2f}, '.format(reg_sa.slope) +
                               ')'}],
        save_kwargs={
            'bbox_inches': 'tight',
            'orientation': 'portrait'},
        axes_functions={
            'set_title': '',
            'set_xlabel': r'clr-dSWA/dT (W m$^{-2}$ K$^{-1}$)',
            'set_ylabel': r'dL$_{\rm v}$P/dT (W m$^{-2}$ K$^{-1}$)',
            'set_xlim': [0.3, 1.35],
            'set_xticks': np.linspace(0.4, 1.2, 5),
            'set_ylim': [1.75, 2.8],
            'set_yticks': np.linspace(1.8, 2.8, 6),
            'legend': {'loc': 3}})

    # filepath2 = os.path.join(cfg[n.PLOT_DIR], 'exfig2b.' +
    #                          cfg[n.OUTPUT_FILE_TYPE])

    e.plot.scatterplot(
        [sa_rsnstcsdt, np.linspace(0.2, 1.4, 2), np.linspace(0.2, 1.4, 2)],
        [sa_rsnstdt, y_reg_rsnst, np.linspace(0.2, 1.4, 2)],
        os.path.join(cfg[n.PLOT_DIR], 'exfig2b.' +
                     cfg[n.OUTPUT_FILE_TYPE]),
        mpl_style_file='default.mplstyle',
        plot_kwargs=[{'linestyle': 'none', 'marker': '*',
                      'markersize': 15, 'markeredgewidth': 2.0,
                      'markerfacecolor': 'k', 'markeredgecolor': 'k'},
                     {'color': 'r', 'linestyle': '-',
                      'label': 'Fit (r={:.2f}, '.format(reg_rsnst.rvalue) +
                               ' slope = {:.2f}, '.format(reg_rsnst.slope) +
                               ')'},
                     {'color': 'tab:gray', 'linestyle': '-',
                      'label': '1:1 line'}],
        save_kwargs={
            'bbox_inches': 'tight',
            'orientation': 'portrait'},
        axes_functions={
            'set_title': '',
            'set_xlabel': r'clr-dSWA/dT (W m$^{-2}$ K$^{-1}$)',
            'set_ylabel': r'all-dSWA/dT (W m$^{-2}$ K$^{-1}$)',
            'set_xlim': [0.45, 1.15],
            'set_xticks': np.linspace(0.5, 1.1, 7),
            'set_ylim': [0.45, 1.15],
            'set_yticks': np.linspace(0.5, 1.1, 7),
            'legend': {'loc': 2}})


def plot_slope_regression_all(cfg, data_model):
    """Scatter plot of linear regression slope, all variables (exfig2a)."""
    if not (cfg[n.WRITE_PLOTS] and cfg.get('plot_ecs_regression')):
        return

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

    plt.style.use('/work/bd0854/b380216/esmvaltool/v2private/' +
                  'ESMValTool-private/esmvaltool/diag_scripts/' +
                  'testkw/style_kw_deangelis2.mplstyle')

    fig, axx = plt.subplots()
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
    axx.set_xticklabels((r"L$_v$dP/dT", "all-dLWC/dT", "all-dSWA/dT",
                         "dSH/dT", "clr-dLWC/dT", "clr-dSWA/dT"),
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
    fig.savefig(os.path.join(cfg[n.PLOT_DIR], 'exfig2a.' +
                             cfg[n.OUTPUT_FILE_TYPE]))
    plt.close()


def plot_rlnst_regression(cfg, dataset_name, data, variables, regs):
    """Plot linear regression used to calculate ECS."""
    if not (cfg[n.WRITE_PLOTS] and cfg.get('plot_ecs_regression')):
        return

    filepath = os.path.join(cfg[n.PLOT_DIR],
                            dataset_name + '.' + cfg[n.OUTPUT_FILE_TYPE])

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

    lab4 = 'LvP: slope (dLvP/dT) = {:.2f}, '.format(regs["lvp"].slope) + \
        'y−int  = {:.2f}, '.format(regs["lvp"].intercept) + \
        'r = {:.2f}'.format(regs["lvp"].rvalue)
    lab = 'LWC: slope (dLWC/dT) = {:.2f}, '.format(regs["rlnst"].slope) + \
        'y−int  = {:.2f}, '.format(regs["rlnst"].intercept) + \
        'r = {:.2f}'.format(regs["rlnst"].rvalue)
    lab2 = 'SWA: slope (dSWA/dT) = {:.2f}, '.format(regs["rsnst"].slope) + \
        'y−int  = {:.2f}, '.format(regs["rsnst"].intercept) + \
        'r = {:.2f}'.format(regs["rsnst"].rvalue)
    lab3 = 'SH: slope (dSH/dT) = {:.2f}, '.format(regs["hfss"].slope) + \
        'y−int  = {:.2f}, '.format(regs["hfss"].intercept) + \
        'r = {:.2f}'.format(regs["hfss"].rvalue)

    axhline_dict = {'y': 0, 'linestyle': 'dashed', 'color': 'k',
                    'linewidth': 2.0}

    e.plot.scatterplot(
        [data[0], np.linspace(0.0, 7.0, 2), data[0], np.linspace(0.0, 7.0, 2),
         data[0], np.linspace(0.0, 7.0, 2), data[0], np.linspace(0.0, 7.0, 2)],
        [data[4], y_reg4, data[1], y_reg, data[2], y_reg2, data[3], y_reg3],
        filepath,
        mpl_style_file='/work/bd0854/b380216/esmvaltool/v2private/' +
        'ESMValTool-private/esmvaltool/diag_scripts/testkw/' +
        'style_kw_deangelisex1.mplstyle',
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
                     'orientation': 'portrait'},
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
    model_nrs = dict([('ACCESS1-0', 1), ('ACCESS1-3', 2),
                      ('bcc-csm1-1', 3),
                      ('bcc-csm1-1-m', 4), ('CanESM2', 5),
                      ('CCSM4', 6), ('CNRM-CM5', 7),
                      ('CNRM-CM5-2', 8), ('GFDL-CM3', 9),
                      ('GFDL-ESM2G', 10),
                      ('GFDL-ESM2M', 11), ('GISS-E2-H', 12),
                      ('GISS-E2-R', 13),
                      ('HadGEM2-ES', 14), ('inmcm4', 15),
                      ('IPSL-CM5A-LR', 16),
                      ('IPSL-CM5A-MR', 17), ('IPSL-CM5B-LR', 18),
                      ('MIROC-ESM', 19), ('MIROC5', 20),
                      ('MPI-ESM-LR', 21),
                      ('MPI-ESM-MR', 22), ('MPI-ESM-P', 23),
                      ('MRI-CGCM3', 24), ('NorESM1-M', 25)])
    pathlist = data.get_path_list(short_name='tas', exp=PICONTROL)
    regressions = np.zeros((len(pathlist), 7))

    data_var = OrderedDict()
    reg_var = OrderedDict()
    varvar = var.short_names()

    for iii, dataset_path in enumerate(pathlist):

        # Substract piControl experiment from abrupt4xCO2 experiment
        dataset = data.get_info(n.DATASET, dataset_path)

        for jvar in varvar:
            data_var[jvar] = data.get_data(short_name=jvar, exp=ABRUPT4XCO2,
                                           dataset=dataset) - \
                data.get_data(short_name=jvar, exp=PICONTROL, dataset=dataset)

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
                            reg_var["rsnstcs"].slope,
                            model_nrs[dataset]]

    return regressions


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
    logging.debug("Found variables in recipe:\n%s", var)

    # Check for tas and rlnst
    if not var.vars_available('tas', 'lvp', 'rlnst', 'rsnst', 'hfss'):
        raise ValueError("This diagnostic needs 'tas', 'lvp', 'rlnst'," +
                         "'rsnst', and 'hfss' variables")

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

    data_model = substract_and_reg_deangelis2(cfg, data, var)

    plot_slope_regression(cfg, data_model)
    plot_slope_regression_all(cfg, data_model)


if __name__ == '__main__':
    with e.run_diagnostic() as config:
        main(config)
