#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###############################################################################
ipcc_ar5/ch09_fig09-42a.py
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


import esmvaltool.diag_scripts.shared as e

import iris

from datetime import datetime
from scipy import stats
import logging
import numpy as np
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    """This is the main routine of the diagnostic.

    Parameters
    ----------
    cfg : dict
        Configuration dictionary of the namelist.

    """
    ###########################################################################
    # Setup diagnostic
    ###########################################################################

    # Model data containers
    MODELS = e.Models(cfg)
    logging.info("Found models in namelist:\n{}".format(MODELS))

    # Variables
    VARS = e.Variables(cfg)
    logging.info("Found variables in namelist:\n{}".format(VARS))
    ECS = e.Variable('ecs', 'ecs', 'equilibrium climate sensitivity', 'K')
    VARS.add_var(ecs=ECS)

    # Experiments
    PICONTROL = 'piControl'
    HISTORICAL = 'historical'
    ABRUPT4XCO2 = 'abrupt4xCO2'
    DIFF = 'difference of abrupt4xCO2 and piControl'
    PICONTROL_TEMP_MEAN = 'total temporal mean of piControl'

    # Matplotlib instance
    style_file = e.plot.get_path_to_mpl_style('default.mplstyle')
    plt.style.use(style_file)
    fig, ax = plt.subplots()

    ###########################################################################
    # Read data
    ###########################################################################

    # Create iris cube for each model
    for model_path in MODELS:
        cube = iris.load(model_path, VARS.standard_names())[0]

        # Global mean
        for coord in [cube.coord(e.LAT), cube.coord(e.LON)]:
            if (not coord.has_bounds()):
                coord.guess_bounds()
        area_weights = iris.analysis.cartography.area_weights(cube)
        cube = cube.collapsed([e.LAT, e.LON], iris.analysis.MEAN,
                              weights=area_weights)

        # Historical: total temporal mean; else: annual mean
        if (MODELS.get_exp(model_path) == HISTORICAL):
            cube = cube.collapsed([e.TIME], iris.analysis.MEAN)
        else:
            cube = cube.aggregated_by(e.YEAR, iris.analysis.MEAN)

        MODELS.set_data(cube.data, model_path=model_path)

    ###########################################################################
    # Process data
    ###########################################################################

    # Substract piControl experiment from abrupt4xCO2 experiment and add total
    # temporal mean of piControl
    for model_path in MODELS.get_path_list(exp=PICONTROL):
        model = MODELS.get_model(model_path)
        short_name = MODELS.get_short_name(model_path)
        data = MODELS.get_data(model_path=model_path)
        data_diff = MODELS.get_data(short_name=short_name, exp=ABRUPT4XCO2,
                                    model=model) - data
        MODELS.add_model(short_name+DIFF+model, data_diff,
                         short_name=short_name, exp=DIFF, model=model)
        MODELS.add_model(short_name+PICONTROL_TEMP_MEAN+model, np.mean(data),
                         short_name=short_name, exp=PICONTROL_TEMP_MEAN,
                         model=model)

    # Calculate ECS (cf. Andrews et al. 2015)
    for model_path in MODELS.get_path_list(short_name=VARS.tas, exp=DIFF):
        model = MODELS.get_model(model_path)
        data_tas = MODELS.get_data(model_path=model_path)
        data_rtmt = MODELS.get_data(short_name=VARS.rtmt, exp=DIFF,
                                    model=model)
        reg_stats = stats.linregress(data_tas, data_rtmt)
        data_ecs = -reg_stats.intercept / (2*reg_stats.slope)
        MODELS.add_model(short_name+DIFF+model+VARS.ecs, data_ecs,
                         short_name=VARS.ecs, exp=DIFF, model=model)

        # Plot ECS regression if desired
        if (cfg['write_plots'] and cfg.get('plot_ecs_regression')):

            # Plot data
            ax.plot(data_tas, data_rtmt, linestyle='none',
                    markeredgecolor='b', markerfacecolor='none',
                    marker='s')

            # Plot regerssion line
            x_reg = np.linspace(-1.0, 8.0, 2)
            y_reg = reg_stats.slope*x_reg + reg_stats.intercept
            ax.plot(x_reg, y_reg, color='k', linestyle='-')

            # Options
            ax.set_title(model)
            ax.set_xlabel(VARS.TAS.standard_name + " / " + VARS.TAS.units)
            ax.set_ylabel(VARS.RTMT.standard_name + " / " + VARS.RTMT.units)
            ax.set_xlim(0.0, 7.0)
            ax.set_ylim(-2.0, 10.0)
            ax.axhline(linestyle='dotted', c='black')
            ax.text(0.05, 0.05,
                    "r = {:.2f}".format(reg_stats.rvalue),
                    transform=ax.transAxes)
            ax.text(0.05, 0.9,
                    r"$\alpha$ = {:.2f},  ".format(-reg_stats.slope) +
                    "F = {:.2f},  ".format(reg_stats.intercept) +
                    "ECS = {:.2f}".format(data_ecs),
                    transform=ax.transAxes)

            # Save plot
            filename = model + '.' + cfg['output_file_type']
            filepath = os.path.join(cfg['plot_dir'], filename)
            fig.savefig(filepath, bbox_inches='tight', orientation='landscape')
            logger.info("Writing {}".format(filepath))
            ax.cla()

    ###########################################################################
    # Plot data
    ###########################################################################

    if (cfg['write_plots']):
        for model_path in MODELS.get_path_list(short_name=VARS.ecs, exp=DIFF):
            model = MODELS.get_model(model_path)
            data_ecs = MODELS.get_data(model_path=model_path)
            data_tas_hist = MODELS.get_data(short_name=VARS.tas,
                                            exp=HISTORICAL, model=model)
            data_tas_piC = MODELS.get_data(short_name=VARS.tas,
                                           exp=PICONTROL_TEMP_MEAN,
                                           model=model)
            style = e.plot.get_model_style(model)

            # Plot
            ax.plot(data_ecs, data_tas_hist, linestyle='none',
                    markeredgecolor=style['color'],
                    markerfacecolor=style['facecolor'],
                    marker=style['mark'], markersize=10, label=model)
            ax.plot(data_ecs, data_tas_piC, linestyle='none',
                    markeredgecolor=style['color'],
                    markerfacecolor=style['facecolor'],
                    marker=style['mark'], markersize=6, label='_'+model)

        # Options
        ax.set_title("GMSAT vs. ECS for CMIP5 models")
        ax.set_xlabel(VARS.ECS.standard_name + " / " + VARS.ECS.units)
        ax.set_ylabel(VARS.TAS.standard_name + " / " + VARS.TAS.units)
        ax.set_xlim(1.5, 5.0)
        legend = ax.legend(loc='center left',
                           bbox_to_anchor=(1.05, 0.5), borderaxespad=0.0,
                           ncol=2)

        # Save plot
        filename = 'ch09_fig09-42a.' + cfg['output_file_type']
        filepath = os.path.join(cfg['plot_dir'], filename)
        fig.savefig(filepath, additional_artists=[legend],
                    bbox_inches='tight', orientation='landscape')
        logger.info("Writing {}".format(filepath))
        ax.cla()
    plt.close()

    ###########################################################################
    # Write nc file
    ###########################################################################

    if (cfg['write_netcdf']):
        data_ecs = MODELS.get_data_list(short_name=VARS.ecs, exp=DIFF)
        models = [model_info[e.MODEL_STR] for model_info in
                  MODELS.get_model_info_list(short_name=VARS.ecs, exp=DIFF)]
        model_coord = iris.coords.AuxCoord(models, long_name='models')
        attr = {'created_by': 'ESMValTool version {}'.format(cfg['version']) +
                              ', diagnostic {}'.format(cfg['script']),
                'creation_date': datetime.utcnow().isoformat(' ') + ' UTC'}
        cube = iris.cube.Cube(data_ecs, long_name=VARS.ECS.long_name,
                              var_name=VARS.ecs, units=VARS.ECS.units,
                              attributes=attr)
        cube.add_aux_coord(model_coord, 0)

        # Save file
        filename = VARS.ecs + '.nc'
        filepath = os.path.join(cfg['work_dir'], filename)
        iris.save(cube, filepath)
        logger.info("Writing {}".format(filepath))


if __name__ == '__main__':

    with e.run_diagnostic() as cfg:
        main(cfg)
